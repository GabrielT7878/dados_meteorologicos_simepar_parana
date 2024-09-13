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
import xarray as xr
import requests
import json
from bs4 import BeautifulSoup

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
def load_csv_data(file_path):
    df = pd.read_csv(file_path)
    return df

@st.cache_data(ttl=datetime.timedelta(days=1))  
def days_without_rain(df,last_date_update):
    df = df[df['Precipitação Acumulada'] >= 5]
    if df.size > 0:
        return (last_date_update - df['Data'].to_list()[0]).days, df['Data'].to_list()[0], df['Precipitação Acumulada'].to_list()[0]
    
    return None

@st.cache_data()
def get_lat_lon(city):
    location = region[region['Cidade'] == city]
    if not location.empty:
        return location["latitude"].values[0], location["longitude"].values[0]
    
    return None, None


def create_city_map(city_option):
    latitude, longitude = get_lat_lon(city_option)

    city_map = folium.Map(location=[latitude, longitude], zoom_start=13, min_zoom=6, max_zoom=17)

    folium.Marker([latitude, longitude], tooltip=city_option).add_to(city_map)

    return city_map

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


@st.cache_data(ttl=datetime.timedelta(days=1))
def download_nc_data_from_source(url):
    url_split = url.split('/')
    file_name = url_split[-1]

    req = requests.get(url, allow_redirects=True)
    
    with open(file_name, 'wb') as file:
        file.write(req.content)

    data = xr.open_dataset(file_name)

    return data

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
    

@st.cache_data(ttl=datetime.timedelta(days=1))
def get_current_coffe_price():

    url = 'https://www.cepea.esalq.usp.br/br/indicador/cafe.aspx'
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 11.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }

    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        df_coffe_value = pd.read_html(str(soup))
        coffe_value = df_coffe_value[0].loc[0]['Valor R$']
        return coffe_value

    return None

@st.cache_data
def calculate_days_without_rain_for_region(region,last_update):

    cities_lat_long = region
    cidade = []
    dias_sem_chuva = []

    date_list = [pd.to_datetime(last_update) - datetime.timedelta(days=x) for x in range(250)]

    for x in cities_lat_long['Cidade'].unique():
        city =  cities_lat_long[cities_lat_long['Cidade'] == x]

        precipitation_on_interval = cpc_data.sel(time=date_list,lat=city['latitude'].values[0], lon=360 + city['longitude'].values[0], method='nearest')

        precip_bigger_than_5mm = precipitation_on_interval.where(precipitation_on_interval['precip']>5)

        precip_bigger_than_5mm = precip_bigger_than_5mm.dropna(dim='time')

        cidade.append(x)
        dias_sem_chuva.append(int((pd.to_datetime(last_update) - pd.to_datetime(precip_bigger_than_5mm['time'].data[0])).days))


    df_cpc = pd.DataFrame({
        "Cidade": cidade,
        "Dias Sem Chuva": dias_sem_chuva
    })

    return df_cpc

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

oeste_parana = load_csv_data("./data/lat_lon_oeste_parana.csv")
sul_minas = load_csv_data("./data/lat_lon_sul_minas.csv")
metropolitana_sao_paulo = load_csv_data("./data/lat_lon_rmsp.csv")

cpc_data = download_nc_data_from_source('https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/precip.2024.nc')
cpc_temperature_max = download_nc_data_from_source('https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/tmax.2024.nc')
cpc_temperature_min = download_nc_data_from_source('https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/tmin.2024.nc')

lavouras_totais = pd.read_csv('data/agricola_cidades_MG_PR_SP.csv')

secao1 = st.container()
col_map_city, col_analyse_data_from_city = secao1.columns([1, 1], gap='large')
col_analyse_data_from_city.header('Dados Meteorológicos')

with col_map_city:
    selected_region = st.radio('Selecione a Região:',['Oeste Paraná','Sul de Minas','Metropolitana de São Paulo '],index=0, horizontal=True)

    if selected_region == 'Oeste Paraná':
        region = oeste_parana
        index = 1
    elif selected_region == 'Sul de Minas':
        region = sul_minas
        index = 4
    else:
        region = metropolitana_sao_paulo
        index = 1

    selected_station = st.selectbox(
        "Cidade",
        region['Cidade'].unique(),
        index=index
    )

    col_map_city.header(selected_station)


    pr_map = create_city_map(selected_station)
    render_folium_map(pr_map)

    df_agricola = pd.read_csv('./data/dados_agricola.csv')

    df_agricola['Município'] = df_agricola['Município'].str.strip()

    agricola_estacao = df_agricola[df_agricola['Município'] == selected_station]

    culturas_comums = lavouras_totais[(lavouras_totais['Nome_Município'] == selected_station) & (lavouras_totais['produto'] != 'Café (em grão) Arábica')].sort_values(by='area_plantada',ascending=False).reset_index()['produto'].values

    area_total = lavouras_totais[(lavouras_totais['Nome_Município'] == selected_station) & (lavouras_totais['produto'] != 'Café (em grão) Arábica')].agg({'area_plantada':'sum'}).values[0]

    if area_total > 0:
        area_message = f"🌾 {num_to_human(area_total)} Hectares"
    else:
        area_message = "Dados não disponíveis"

    
    st.title("Dados de Agricultura")

    st.markdown("---")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Área Cultivada")
        st.subheader(area_message)
    with col2:
        if len(culturas_comums):
            st.subheader("Culturas Cultivadas Mais Comuns")
            if len(culturas_comums) < 3:
                qtd_culturas = len(culturas_comums)
            else:
                qtd_culturas = 3
            for i in range(qtd_culturas):
                st.subheader(f'{i+1}º {culturas_comums[i]}')
        else:
            st.subheader(f'Sem Dados')   

    st.markdown("---")


with col_analyse_data_from_city:

    cpc_last_update = pd.to_datetime(cpc_data['time'][cpc_data["time"].size-1].data)
    cpc_date_list = [cpc_last_update - datetime.timedelta(days=x) for x in range(100)]

    lat, lon = get_lat_lon(selected_station)

    precipitation_on_interval = cpc_data.sel(time=cpc_date_list,lat=lat, lon=360 + lon, method='nearest')

    df_cpc = pd.DataFrame({
        "Data": precipitation_on_interval['time'].data,
        "Precipitação Acumulada": precipitation_on_interval['precip'].data
    })

    count_days_without_rain = days_without_rain(df_cpc, cpc_last_update)

    tabs = st.tabs(['Chuva', 'Temperatura'])

    with tabs[0]:
        ct = st.container()
        col1, col2 = ct.columns([1, 1], gap='large')
        with col1:
            if count_days_without_rain:
                st.write(f'# Dias sem chuva: {count_days_without_rain[0]}')
                st.write(f'Último dia com chuva: {count_days_without_rain[1]} : {count_days_without_rain[2]:.2f}mm')
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

        with ct_grap:
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
                lat, lon = get_lat_lon(selected_station)

                precipitation_on_interval = cpc_data.sel(time=date_list,lat=lat, lon=360 + lon, method='nearest')
                df_cpc = pd.DataFrame({
                    "Data": precipitation_on_interval['time'].data,
                    "Precipitação Acumulada": precipitation_on_interval['precip'].data
                })

                createDailyChart(df_cpc ,metrics[4])
                
                #MAP
                secao2 = st.container()
                secao2.header(f'Mapa Dias sem Chuva - Região: {selected_region}')

                df_cpc = calculate_days_without_rain_for_region(region,cpc_last_update)

                st.subheader(f"Cidades com {count_days_without_rain[0]} dias sem chuva: {round(len(df_cpc[df_cpc['Dias Sem Chuva'] == count_days_without_rain[0]]) /  len(df_cpc) * 100)} %")            
                coords_map = None

                if selected_region == 'Oeste Paraná':
                    region_code = 41
                    coords_map = [-24.292611466381025, -53.31514183892323]
                elif selected_region == 'Sul de Minas':
                    region_code = 31
                    coords_map = [-21.544771753269067, -45.452079858011935]
                else:
                    region_code = 35
                    coords_map = [-23.51712720453375, -46.58213425292997]
                
                cities_geojson_data = load_json(f"./data/{region_code}_geo.json")
                
                cities_lat_long = pd.read_csv("./data/lat_lon_sul_minas.csv",sep=",")
            
                accumulated_precipitation_by_city = cities_lat_long.merge(df_cpc,on="Cidade",how="outer")

                mg_map = folium.Map(location=coords_map, zoom_start=8)

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
                    name='Dias Sem Chuva',
                    data=df_cpc,
                    columns=['Cidade', 'Dias Sem Chuva'],
                    key_on='feature.properties.name',
                    fill_color='YlOrRd',
                    fill_opacity=0.7,
                    line_opacity=0.2,
                    legend_name='Dias Sem Chuva',
                    nan_fill_opacity=0.0,
                    nan_fill_color="white",
                    highlight=True
                ).add_to(mg_map)

                accumulated_precipitation_by_city = accumulated_precipitation_by_city.dropna()

                render_folium_map(mg_map,height=600)

                #st.dataframe(df_cpc,use_container_width=True,height=200)
             

    with tabs[1]:

        #CPC TEMP
        st.subheader("CPC")
        numdays = (cpc_period[1] - cpc_period[0]).days + 1 
        date_list = [cpc_period[1] - datetime.timedelta(days=x) for x in range(numdays)]
        lat, lon = get_lat_lon(selected_station)

        cpc_temperature_max = cpc_temperature_max.sel(time=date_list,lat=lat, lon=360 + lon, method='nearest')
        cpc_temperature_min = cpc_temperature_min.sel(time=date_list,lat=lat, lon=360 + lon, method='nearest')
        temp_mean = []

        temp_mean = (cpc_temperature_max['tmax'].data + cpc_temperature_min['tmin'].data)/2

        df_cpc = pd.DataFrame({
            "Data": cpc_temperature_max['time'].data,
            'Temperatura Média' : temp_mean,
        })

        df_cpc_min_max_temp = pd.DataFrame({
            "Data": cpc_temperature_max['time'].data,
            "Temp_Max": cpc_temperature_max['tmax'].data,
            "Temp_Min": cpc_temperature_min['tmin'].data
        })

        df_cpc_min_max_temp["Data"] = df_cpc_min_max_temp["Data"].dt.strftime('%d/%m')

        max_scale = 50

        chart_df_temp_max = alt.Chart(df_cpc_min_max_temp).mark_bar(opacity=1,color='darksalmon').encode(
            y=alt.Y('Temp_Max',scale=alt.Scale(domain=[0, max_scale]), axis=alt.Axis()),
            x=alt.X('Data')
        )

        chart_df_temp_min = alt.Chart(df_cpc_min_max_temp).mark_bar(opacity=0.6,color='aliceblue').encode(
            y=alt.Y('Temp_Min',scale=alt.Scale(domain=[0, max_scale]), axis=alt.Axis(title=f'Temperatura Máxima e Miníma ºC')),
            x=alt.X('Data')
        )

        c = alt.layer(chart_df_temp_max, chart_df_temp_min)

        st.altair_chart(c, use_container_width=True)


simulacao_seguro_container = st.container()

dias_de_corbertura_seguro = 151

with simulacao_seguro_container:
    st.title('📊 Simulação do Seguro')
    

    tabs = st.tabs(['☕ Café', '🌱 Soja'])
    
    with tabs[0]:
        st.header('Seguro Paramétrico com cobertura contra a Seca (Café)')
        st.write('')
        st.markdown("<br>", unsafe_allow_html=True)

        st.header(f'🗺️ Região: {selected_station}')
        st.subheader(f' ')

        col1, col2, col3= st.columns([0.7,0.1,0.3])

        with col1:
            sub_col1, sub_col2 = st.columns(2)

            with sub_col1:
                st.subheader('📅 Período de Cobertura')
                st.metric(label="", value="121 Dias", delta="(Setembro a Dezembro)")
            with sub_col2:
                st.subheader('🌧️ Gatilho de Indenização')
                st.metric(label="", value="Dias sem chuva")

            sub_col3, sub_col4 = st.columns(2)

            with sub_col3:
                st.subheader('🌍 Área Total Segurada')
                st.metric(label="", value=f"{num_to_human(area_total)} hectares")

            with sub_col4:
                st.subheader('💰 Valor Atual da Saca de Café')
                st.metric(label="", value=f"R$ {get_current_coffe_price()}")


        with col2:
        
            st.html(
            '''
                <div class="divider-vertical-line"></div>
                <style>
                    .divider-vertical-line {
                        border-left: 2px solid rgba(49, 51, 63, 0.5);
                        height: 320px;
                        margin: auto;
                    }
                </style>
            '''
            )

        with col3:
            seca_count = 6
            if count_days_without_rain[0] > 12:
                seca_count = 12
            else:
                seca_count = count_days_without_rain[0]       
         
            valor_saca_cafe_arabica = float(get_current_coffe_price().replace('.', '').replace(',', '.'))
            lmi = area_total * valor_saca_cafe_arabica
            valor_indenizacao = round((seca_count / dias_de_corbertura_seguro) * lmi)

            st.subheader(f'🌧️ Dias sem chuva')

            cpc_date_list = [cpc_last_update - datetime.timedelta(days=x) for x in range((pd.to_datetime(cpc_last_update) - pd.to_datetime(datetime.date(2024,8,1))).days)]
            
            lat, lon = get_lat_lon(selected_station)

            precipitation_on_interval = cpc_data.sel(time=cpc_date_list,lat=lat, lon=360 + lon, method='nearest')

            df_cpc = pd.DataFrame({
                "Data": precipitation_on_interval['time'].data,
                "Precipitação Acumulada": precipitation_on_interval['precip'].data
            })


            st.metric(label="", value=f"{seca_count}")
            st.subheader('💵 Valor da Indenização')
            st.write(' ')
            st.header(f":green[R$ {num_to_human(valor_indenizacao)}]")

st.write(' ')
st.write(' ')
st.write(' ')