import altair as alt
import pandas as pd
import streamlit as st
import datetime
import numpy as np
import folium
from streamlit_folium import folium_static
import streamlit.components.v1 as components
from scipy.spatial import cKDTree
import numpy as np
from folium import Choropleth
import cdsapi
import xarray as xr
import requests
import json

# --------------- Page config ---------------
page_title = "Painel Seguro Paramétrico ☁️"
page_icon = "☁️"

st.set_page_config(page_title=page_title, page_icon=page_icon,layout='wide')
st.title(page_title)
# -------------------------------------------

@st.cache_data
def load_json(file_path):
    with open(file_path, 'r') as jsonFile:
        f = json.load(jsonFile)
        return f
    
@st.cache_data
def load_geo_json_data(file_path):
    with open(file_path) as f:
        geojson_data = json.load(f)
        return geojson_data

@st.cache_data(ttl='1d')
def load_csv_data(file_path,sep):
    df = pd.read_csv(file_path,sep=sep)
    return df

@st.cache_data(ttl='1d')  
def days_without_rain(city_option, actual_date, total_period_in_days):
    count = 0
    date = actual_date
    for i in range(0,total_period_in_days):
        day = meteorological_data[(meteorological_data["Data"] == date) & (meteorological_data["Cidade"] == city_option)]
        max_precipitation = day["Precipitação Acumulada"].sum()
        if(max_precipitation == 0):
            count += 1
        else:
            return count
        date -= datetime.timedelta(days=1)

@st.cache_data(ttl='1d')  
def days_without_rain_nc(df,last_date_update):
    df = df[df['Precipitação Acumulada'] >= 5]
    if df.size > 0:
        return (last_date_update - df['Data'].to_list()[0]).days, df['Data'].to_list()[0], df['Precipitação Acumulada'].to_list()[0]
    
    return None

@st.cache_data()
def get_lat_long(city):
    location = region[region['municipio'] == city]
    if not location.empty:
        return location["latitude"].values[0], location["longitude"].values[0]
    
    return None, None


def create_city_map(city_option):
    latitude, longitude = get_lat_long(city_option)

    city_map = folium.Map(location=[latitude, longitude], zoom_start=13, min_zoom=6, max_zoom=17)

    folium.Marker([latitude, longitude], tooltip=city_option).add_to(city_map)

    return city_map

def render_folium_map(map_obj,height=400):
    map_html = map_obj._repr_html_()

    map_html = f"""
    <html>
    <head>
    <style>
    .folium-map {{
        position: relative;
        width: 100%;
        height: 400vh;
    }}
    </style>
    </head>
    <body>
    <div class="folium-map">{map_html}</div>
    </body>
    </html>
    """
    components.html(map_html, height=height)  # Ajuste a altura conforme necessário

@st.cache_data(ttl='1d') 
def interpolateData(data=None):
    known_precipitation = data.dropna(subset=['Precipitação Acumulada'])
    unknown_precipitation = data[data['Precipitação Acumulada'].isna()]

    tree = cKDTree(known_precipitation[['latitude', 'longitude']])

    def inverse_distance_weighting(distances, values, num_nearest=3):
        inverse_distances = 1 / distances
        weights = inverse_distances / inverse_distances.sum()
        return np.dot(weights, values)

    def interpolate_missing_values(row):
        distances, indices = tree.query([row['latitude'], row['longitude']], k=3)
        nearest_values = known_precipitation.iloc[indices]['Precipitação Acumulada'].values
        return inverse_distance_weighting(distances, nearest_values)

    interpolated_values = unknown_precipitation.apply(interpolate_missing_values, axis=1)

    data.loc[data['Precipitação Acumulada'].isna(), 'Precipitação Acumulada'] = interpolated_values

    return data

def createDailyChart(data, metric):

    data["Data"] = data["Data"].dt.strftime('%d/%m')
    
    max_scale = 50

    if data[metric].max() > max_scale:
        max_scale = data[metric].max()

    chart_df = alt.Chart(data).mark_bar().encode(
        y=alt.Y(metric,scale=alt.Scale(domain=[0, max_scale]), axis=alt.Axis(title=f'{metric} {metrics_unit[metric]}')),
        x=alt.X('Data')
    )
    st.altair_chart(chart_df,use_container_width=True)

def createHourlyChart(df_chart,metric):
    
    chart = (
        alt.Chart(df_chart)
        .mark_line()
        .encode(
            x=alt.X("Horario:N", title="Horário"),
            y=alt.Y(f'{metric}:Q', title=f'{metric} {metrics_unit[metric]}'),
            color="Cidade:N"
        )
    )
    st.altair_chart(chart, use_container_width=True)


@st.cache_data(ttl='1d') 
def download_nc_data_from_source(url):
    url_split = url.split('/')
    file_name = url_split[-1]

    req = requests.get(url, allow_redirects=True)
    
    with open(file_name, 'wb') as file:
        file.write(req.content)

    data = xr.open_dataset(file_name)

    return data


@st.cache_data
def request_data_period_from_era5_api(period,lat,lon):
    era5_conn = cdsapi.Client()

    hash_table = {}
    for date in period:
        if not str(date.month) in hash_table:
            hash_table[str(date.month)] = []
        hash_table[str(date.month)].append(date.day)
    
    for key in hash_table.keys():
        era5_conn.retrieve(
            'reanalysis-era5-single-levels',
            {
            'product_type': 'reanalysis',
            'variable': 'total_precipitation',
            'year': '2024',
            'month': key,
            'day': hash_table[key], 
            'time': [
                    '00:00', '01:00', '02:00',
                    '03:00', '04:00', '05:00',
                    '06:00', '07:00', '08:00',
                    '09:00', '10:00', '11:00',
                    '12:00', '13:00', '14:00',
                    '15:00', '16:00', '17:00',
                    '18:00', '19:00', '20:00',
                    '21:00', '22:00', '23:00',
            ],
            'area': [
                lat,lon, 
                lat,lon,
            ],
            'format': 'netcdf',
        },
            'precipitation_data.nc'
        )

        data = xr.open_dataset('precipitation_data.nc')

        date_period = []
        day_precip = []


        for date in period:
            hours = [f'{date}T{x:02}:00:00.000000000' for x in range(24)]
            day = data.sel(time=hours)
            date_period.append(date)

            total_precip = 0
            for precip in day["tp"].data:
                total_precip += precip[0][0]

            total_precip = total_precip * 1000

            day_precip.append(total_precip)

        return pd.DataFrame(
            {
                "Data" : date_period,
                "Precipitação Acumulada" : day_precip
            }
        )

def num_to_human(num):
    if num >= 1_000_000_000:
        num_bilion = num / 1_000_000_000
        return f"{num_bilion:,.3f} bilhões"
    elif num >= 1_000_000:
        num_milion = num / 1_000_000
        return f"{num_milion:,.3f} milhões"
    elif num >= 1_000:
        num_thousand = num / 1_000
        return f"{num_thousand:,.3f} mil"
    else:
        return f"{num}"
    



metrics = ["Temperatura Média",
            "Umidade Relativa",
            "Velocidade do Vento",
            "Rajada do Vento",
            "Precipitação Acumulada",
            "Pressão Atmosférica Reduzida"]

metrics_unit = {
    "Temperatura Média": "(ºC)",
    "Umidade Relativa": "(%)",
    "Velocidade do Vento": "(Km/h)",
    "Rajada do Vento": "(Km/h)",
    "Precipitação Acumulada": "(mm)",
    "Pressão Atmosférica Reduzida": "(hPa)"
}


stations_coordinates = load_json("data/coordenadas_estacoes.json")
meteorological_data = load_csv_data("data/dados_meteorologicos_simepar_parana.csv",sep=',')

oeste_parana = load_csv_data("./data/lat_lon_oeste_parana.csv",sep=',')
sul_minas = load_csv_data("./data/lat_lon_sul_minas.csv",sep=',')




#processing simepar data
meteorological_data.fillna(0, inplace=True)
meteorological_data["Data"] = pd.to_datetime(meteorological_data["Data"])

cpc_data = download_nc_data_from_source('https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/precip.2024.nc')

chirps_data = download_nc_data_from_source('https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p05/chirps-v2.0.2024.days_p05.nc')

cpc_temperature_max = download_nc_data_from_source('https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/tmax.2024.nc')
cpc_temperature_min = download_nc_data_from_source('https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/tmin.2024.nc')

st.write(
    """
    Este painel visualiza dados meteorológicos de 4 fontes: CPC, ERA 5, SIMEPAR e CHIRPS
    """
)

secao1 = st.container()
col_map_city, col_analyse_data_from_city = secao1.columns([0.7, 1], gap='large')
col_analyse_data_from_city.header('Dados Meteorológicos')



with col_map_city:
    selected_region = st.radio('Selecione a Região:',['Oeste Paraná','Sul de Minas'],index=0, horizontal=True)

    if selected_region == 'Oeste Paraná':
        region = oeste_parana
        index = 1
    else:
        region = sul_minas
        index = 4

    selected_station = st.selectbox(
        "Cidade",
        region['municipio'].unique(),
        index=index
    )

    col_map_city.header(selected_station)


    pr_map = create_city_map(selected_station)
    render_folium_map(pr_map)

    population_cities = pd.read_csv('data/populacao_regioes.csv')

    selected_station_population = population_cities[population_cities['Município'] == selected_station]

    population = selected_station_population['População residente - pessoas [2022]'].astype(int).values[0]

    area_city = selected_station_population['Área Territorial - km² [2022]'].astype(int).values[0]

    population_density = selected_station_population['Densidade demográfica - hab/km² [2022]'].astype(int).values[0]

    pib = selected_station_population['PIB per capita - R$ [2021]'].astype(int).values[0]

    propriedades_habitacionais = pd.read_csv('./data/total_propriedades_habitacionais.csv',sep=',')

    total_propriedades_habitacionais = propriedades_habitacionais[propriedades_habitacionais['Município'] == selected_station]['Domicílios recenseados (Domicílios)'].values[0]

    df_agricola = pd.read_csv('./data/dados_agricola.csv')

    df_agricola['Município'] = df_agricola['Município'].str.strip()

    agricola_estacao = df_agricola[df_agricola['Município'] == selected_station]

    population_message = f"🌍 {num_to_human(population)}"
    area_message = f"🌾 {num_to_human(agricola_estacao['Área colhida (Hectares)'].astype(int).values[0] * 10000)} m²"
    valor_producao_message = f"💰 R${num_to_human(agricola_estacao['Valor da produção (Mil Reais)'].astype(int).values[0])}"


    st.title("Dados Agrícolas e População")

    

    st.markdown("---")

    # Informações de População
    col1, col2, col3 = st.columns(3,gap='small')
    with col1:
        st.subheader("População")
        st.subheader(population_message)
    with col2:
        st.subheader("Área Teritorial")
        st.subheader(str(area_city) + ' km²')
    with col3:
        st.subheader("Densidade Populacional")
        st.subheader(str(population_density) + ' (habitantes/km²)')

    col1, col2, col3 = st.columns(3,gap='small')  
    with col1:
        st.subheader("Total de Propriedades Habitacionais")
        st.subheader(num_to_human(total_propriedades_habitacionais))
    with col2:
        st.subheader("PIB")
        st.subheader(f'R${num_to_human(pib * population)}')

    st.markdown("---")

    # Informações de Agricultura
    st.header("Dados de Agricultura")
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Área Cultivada")
        st.subheader(area_message)


    st.markdown("---")



with col_analyse_data_from_city:

    actual_date = meteorological_data["Data"].max()
    total_period_in_days = (actual_date - meteorological_data["Data"].min()).days

    cpc_last_update = pd.to_datetime(cpc_data['time'][cpc_data["time"].size-1].data)
    cpc_date_list = [cpc_last_update - datetime.timedelta(days=x) for x in range(100)]

    lat, lon = get_lat_long(selected_station)

    precipitation_on_interval = cpc_data.sel(time=cpc_date_list,lat=lat, lon=360 + lon, method='nearest')

    df_cpc = pd.DataFrame({
        "Data": precipitation_on_interval['time'].data,
        "Precipitação Acumulada": precipitation_on_interval['precip'].data
    })

    count = days_without_rain_nc(df_cpc, cpc_last_update)

    tabs = st.tabs(['Chuva', 'Temperatura'])

    with tabs[0]:
        ct = st.container()
        col1, col2 = ct.columns([1, 1], gap='large')
        with col1:
            if count:
                st.write(f'# Dias sem chuva: {count[0]}')
                st.write(f'Último dia com chuva: {count[1]} : {count[2]:.2f}mm')
            else:
                st.write(f'# Dias sem chuva: sem dados!')
        with col2:
            st.subheader('Período')
            start_date = st.date_input("Início",datetime.date.today() - datetime.timedelta(days=15),key='all-sources',max_value=datetime.date.today())
            end_date = st.date_input("Fim",datetime.date.today() - datetime.timedelta(days=1),key='all-sources2',max_value=datetime.date.today())

            if start_date > end_date:
                st.write(':red[Período Inválido: Data de Início deve ser menor que a data de Fim]')
                period_all_sources = (start_date,start_date + datetime.timedelta(days=1))
            else:
                period_all_sources = (start_date,end_date)

            
        ct_grap = st.container()
        col_grap1, col_grap2 = ct.columns([1, 1], gap='large')

        with col_grap1:
            #CPC
            cpc_last_update = pd.to_datetime(cpc_data['time'][cpc_data["time"].size-1].data)
            st.markdown(f"**CPC**: última atualização: {pd.to_datetime(cpc_last_update).strftime('%d/%m/%Y')}")
            
            cpc_period = (cpc_last_update - datetime.timedelta(days=15),cpc_last_update)
            
            if pd.Timestamp(period_all_sources[1]) <= pd.Timestamp(cpc_last_update):
                cpc_period = period_all_sources
            else:
                if pd.Timestamp(cpc_period[1]) >= pd.Timestamp(period_all_sources[0]):
                    cpc_period = pd.to_datetime((period_all_sources[0],cpc_period[1]))

            if len(cpc_period) == 2:
                numdays = (cpc_period[1] - cpc_period[0]).days + 1 
                date_list = [cpc_period[1] - datetime.timedelta(days=x) for x in range(numdays)]
                lat, lon = get_lat_long(selected_station)

                precipitation_on_interval = cpc_data.sel(time=date_list,lat=lat, lon=360 + lon, method='nearest')
                df_cpc = pd.DataFrame({
                    "Data": precipitation_on_interval['time'].data,
                    "Precipitação Acumulada": precipitation_on_interval['precip'].data
                })

                createDailyChart(df_cpc ,metrics[4])
                
        with col_grap2:
            #CHIRPS
            chirps_last_update = pd.to_datetime(chirps_data['time'][chirps_data["time"].size-1].data)
            st.markdown(f"**CHIRPS**: última atualização: {pd.to_datetime(chirps_last_update).strftime('%d/%m/%Y')}")
            chirps_period = (chirps_last_update - datetime.timedelta(days=15),chirps_last_update)
            
            if pd.Timestamp(period_all_sources[1]) <= pd.Timestamp(chirps_last_update):
                chirps_period = period_all_sources
            else:
                if pd.Timestamp(chirps_period[1]) >= pd.Timestamp(period_all_sources[0]):
                    chirps_period = pd.to_datetime((period_all_sources[0],chirps_period[1]))

            if len(chirps_period) == 2:

                numdays = (chirps_period[1] - chirps_period[0]).days + 1

                date_list = [chirps_period[1] - datetime.timedelta(days=x) for x in range(numdays)]

                lat, lon = get_lat_long(selected_station)

                precipitation_on_interval = chirps_data.sel(time=date_list,latitude=lat, longitude=lon, method='nearest')

                df_cpc = pd.DataFrame({
                    "Data": precipitation_on_interval['time'].data,
                    "Precipitação Acumulada": precipitation_on_interval['precip'].data
                })

                createDailyChart(df_cpc, metrics[4]) 
        
        ct_grap2 = st.container()
        col_grap3, col_grap4 = ct.columns([1, 1], gap='large')
        with col_grap3:
            #SIMEPAR
            if selected_station.lower() in [x.lower() for x in meteorological_data['Cidade'].unique()]:
                simepar_last_update = meteorological_data["Data"].max().date()
                st.markdown(f"**SIMEPAR**: última atualização: {pd.to_datetime(simepar_last_update).strftime('%d/%m/%Y')}")

                simepar_period = (simepar_last_update - datetime.timedelta(days=15),simepar_last_update)
                
                if pd.Timestamp(period_all_sources[1]) <= pd.Timestamp(simepar_last_update):
                    simepar_period = period_all_sources
                else:
                    if pd.Timestamp(simepar_period[1]) >= pd.Timestamp(period_all_sources[0]):
                        simepar_period = pd.to_datetime((period_all_sources[0],simepar_period[1]))
                    
                
                if len(simepar_period) == 2:
                    df_filtered = meteorological_data[(meteorological_data["Cidade"].str.lower() == selected_station.lower()) & (meteorological_data["Data"] >= pd.to_datetime(simepar_period[0])) & (meteorological_data["Data"] <= pd.to_datetime(simepar_period[1]))]
                    df_simepar = df_filtered.groupby(['Data']).agg({metrics[4] : 'sum'}).reset_index()
                    
                    createDailyChart(df_simepar, metrics[4])
            else:
                st.write("Simepar: Sem dados para a cidade selecionada!")               
            
        with col_grap4:
            #ERA 5
            era5_req_date = requests.get("https://cds.climate.copernicus.eu/api/v2.ui/resources/reanalysis-era5-single-levels")
            date_string = era5_req_date.json()['update_date']
            format = "%Y-%m-%d"
            era5_last_update = datetime.datetime.strptime(date_string, format) - datetime.timedelta(days=6) 
            st.markdown(f"**ERA 5**: última atualização: {pd.to_datetime(era5_last_update).strftime('%d/%m/%Y')}")

            era5_period = (era5_last_update - datetime.timedelta(days=15),era5_last_update)
            
            if pd.Timestamp(period_all_sources[1]) <= pd.Timestamp(era5_last_update):
                era5_period = period_all_sources
            else:
                if pd.Timestamp(era5_period[1]) >= pd.Timestamp(period_all_sources[0]):
                    era5_period = pd.to_datetime((period_all_sources[0],era5_period[1]))
            
            if len(era5_period) == 2:
                numdays = (era5_period[1] - era5_period[0]).days + 1

                date_list = [era5_period[1] - datetime.timedelta(days=x) for x in range(numdays)]    

                lat, lon = get_lat_long(selected_station)

                df_era5 = request_data_period_from_era5_api(date_list,lat,lon)
                df_era5['Data'] = pd.to_datetime(df_era5['Data'])

                createDailyChart(df_era5, metrics[4])
    with tabs[1]:

        #CPC TEMP
        st.subheader("CPC")
        numdays = (cpc_period[1] - cpc_period[0]).days + 1 
        date_list = [cpc_period[1] - datetime.timedelta(days=x) for x in range(numdays)]
        lat, lon = get_lat_long(selected_station)

        cpc_temperature_max = cpc_temperature_max.sel(time=date_list,lat=lat, lon=360 + lon, method='nearest')
        cpc_temperature_min = cpc_temperature_min.sel(time=date_list,lat=lat, lon=360 + lon, method='nearest')
        temp_mean = []


        temp_mean = (cpc_temperature_max['tmax'].data + cpc_temperature_min['tmin'].data)/2

        df_cpc = pd.DataFrame({
            "Data": cpc_temperature_max['time'].data,
            'Temperatura Média' : temp_mean,
        })

        createDailyChart(df_cpc,metrics[0])

        #SIMEPAR
        st.subheader("Simepar")
        if selected_station.lower() in [x.lower() for x in meteorological_data['Cidade'].unique()]:
            df_simepar = df_filtered.groupby(['Data']).agg({'Temperatura Média':'mean'}).reset_index()
            createDailyChart(df_simepar,metrics[0])
        else:
            st.write("Simepar: Sem dados para a cidade selecionada!")   




secao2 = st.container()
col_precipitation_map, col_2 = secao2.columns([1, 1], gap='large')
col_precipitation_map.header('Precipitação acumulada no Período - (Paraná)')
col_2.header('Precipitação acumulada no Período - (Minas Gerais)')

meteorological_data = load_csv_data("data/dados_meteorologicos_simepar_parana.csv",sep=',')
meteorological_data.fillna(0, inplace=True)
meteorological_data["Data"] = pd.to_datetime(meteorological_data["Data"])

with col_precipitation_map:
    selected_period_days = st.selectbox(
            "Período",
            [7,15,30],
            format_func=lambda x: "últimos " + str(x) + " dias",
            index=0
        )

    interval_period = meteorological_data["Data"].max() - datetime.timedelta(days=selected_period_days), meteorological_data["Data"].max()

    df_period_selected = meteorological_data[(meteorological_data["Data"] >= interval_period[0]) & (meteorological_data["Data"] <= interval_period[1])]

    # Estimate precipitation on all cities 
    cities_geojson_data = load_geo_json_data("./data/geojs-41-mun.json")

    cities_lat_long = pd.read_csv("./data/lat_long_cidades_pr.csv",sep=",")

    df_period_total =df_period_selected.groupby(['Cidade']).agg({'Precipitação Acumulada': 'sum'}).reset_index()

    accumulated_precipitation_by_station = df_period_total[["Cidade","Precipitação Acumulada"]]

    accumulated_precipitation_by_city = cities_lat_long.merge(accumulated_precipitation_by_station,on="Cidade",how="outer")

    accumulated_precipitation_estimation = interpolateData(accumulated_precipitation_by_city)


    #Precipitation map in the period
    pr_map = folium.Map(location=[-24.5452, -51.5652], zoom_start=7)

    style = lambda x: {"color" : "white",
                    "fillOpacity": 1,
                    "weight": 1}

    folium.GeoJson(
        cities_geojson_data,
        name='geojson',
        style_function=style
    ).add_to(pr_map)


    Choropleth(
        geo_data=cities_geojson_data,
        name='precipitação acumulada no período',
        data=accumulated_precipitation_estimation,
        columns=['Cidade', 'Precipitação Acumulada'],
        key_on='feature.properties.name',
        fill_color='YlOrRd',
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name='Precipitação Acumulada (mm)',
        nan_fill_opacity=0.0,
        nan_fill_color="white"
    ).add_to(pr_map)

    render_folium_map(pr_map,height=600)

with col_2:
    selected_period_days = st.selectbox(
            "Período",
            [7,15,30],
            format_func=lambda x: "últimos " + str(x) + " dias",
            index=0,
            key='precip_MG'
        )

    cities_geojson_data = load_geo_json_data("./data/geojs-31-mun.json")

    cities_lat_long = pd.read_csv("./data/lat_long_cidades_mg.csv",sep=",")

    date_list = [cpc_period[1] - datetime.timedelta(days=x) for x in range(selected_period_days)]

    df_cpc = pd.DataFrame({
        "Data": [],
        "Precipitação Acumulada": []
    })

    for x in cities_lat_long['Cidade'].unique():
        city =  cities_lat_long[cities_lat_long['Cidade'] == x]

        precipitation_on_interval = cpc_data.sel(time=date_list,lat=city['latitude'].values[0], lon=360 + city['longitude'].values[0], method='nearest')
        
        city_data = {
            "Cidade": x,
            "Precipitação Acumulada": np.sum(precipitation_on_interval['precip'].data)
        }
        df_cpc = pd.concat([df_cpc, pd.DataFrame([city_data])], ignore_index=True)
   
    accumulated_precipitation_by_city = cities_lat_long.merge(df_cpc,on="Cidade",how="outer")

    #Precipitation map in the period
    mg_map = folium.Map(location=[-18.64288124580271, -44.686599834296274], zoom_start=6)

    style = lambda x: {"color" : "white",
                    "fillOpacity": 1,
                    "weight": 1}

    folium.GeoJson(
        cities_geojson_data,
        name='geojson',
        style_function=style
    ).add_to(mg_map)


    Choropleth(
        geo_data=cities_geojson_data,
        name='precipitação acumulada no período',
        data=accumulated_precipitation_by_city,
        columns=['Cidade', 'Precipitação Acumulada'],
        key_on='feature.properties.name',
        fill_color='YlOrRd',
        fill_opacity=0.7,
        line_opacity=0.2,
        legend_name='Precipitação Acumulada (mm)',
        nan_fill_opacity=0.0,
        nan_fill_color="white"
    ).add_to(mg_map)

    render_folium_map(mg_map,height=600)