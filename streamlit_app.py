import altair as alt
import pandas as pd
import streamlit as st
import datetime
import numpy as np
import folium
from streamlit_folium import folium_static
from scipy.spatial import cKDTree
import numpy as np
from folium import Choropleth
import json

# --------------- Page config ---------------
page_title = "Dados Meteorológicos Simepar Paraná ☁️"
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
    
def load_csv_data(file_path):
    df = pd.read_csv(file_path)
    return df

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

def get_lat_long(city):
    location = stations_coordinates[city]
    if location:
        return location["latitude"], location["longitude"]
    else:
        return None, None


def create_city_map(city_option):
    latitude, longitude = get_lat_long(city_option)

    city_map = folium.Map(location=[latitude, longitude], zoom_start=12, min_zoom=6, max_zoom=13)

    folium.Marker([latitude, longitude], tooltip=city_option).add_to(city_map)

    return city_map

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

def createDailyChart(station,metric):
    period = meteorological_data["Data"].max()-datetime.timedelta(days=7), meteorological_data["Data"].max()

    df_filtered = meteorological_data[(meteorological_data["Cidade"] == station) & (meteorological_data["Data"] >= period[0]) & (meteorological_data["Data"] <= period[1])]
    df_filtered = df_filtered.groupby(['Data']).agg({metric : 'sum'}).reset_index()

    df_filtered["Data"] = df_filtered["Data"].dt.strftime('%d/%m')

    st.bar_chart(df_filtered, x="Data", y=metric)

def createHourlyChart(station,metric):
    ActualDate = meteorological_data["Data"].max()

    df_filtered = meteorological_data[(meteorological_data["Cidade"] == station) & (meteorological_data["Data"] == pd.to_datetime(ActualDate))]
    
    df_reshaped = df_filtered.pivot_table(
        index="Horario", columns="Cidade", values=metric, fill_value=0
    )

    df_chart = pd.melt(
        df_reshaped.reset_index(), id_vars="Horario", var_name="Cidade", value_name=metric
    )


    chart = (
        alt.Chart(df_chart)
        .mark_point()
        .encode(
            x=alt.X("Horario:N", title="Horário"),
            y=alt.Y(f'{metric}:Q', title=f'{metric} {metrics_unit[metric]}'),
            color="Cidade:N"
        )
        .properties(height=320,title=f"Atualização {ActualDate.strftime('%d/%m/%Y')}")
    )
    st.altair_chart(chart, use_container_width=True)

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
meteorological_data = load_csv_data("data/dados_meteorologicos_simepar_parana.csv")

#processing data
meteorological_data.fillna(0, inplace=True)
meteorological_data["Data"] = pd.to_datetime(meteorological_data["Data"])

st.write(
    """
    Este aplicativo visualiza dados do [SIMEPAR](http://www.simepar.br/prognozweb/simepar).
    Ele mostra métricas sobre o clima coletadas em estações meteorológicas no Paraná-RS.
    """
)

st.write(f'Registros na base de dados: {meteorological_data["Data"].min().date()} até {meteorological_data["Data"].max().date()}')

secao1 = st.container()
col_map_city, col_analyse_data_from_city = secao1.columns([1, 1], gap='large')
col_analyse_data_from_city.header('Dados das Estações')
col_map_city.header("Mapa da Cidade")

with col_analyse_data_from_city:
    selected_station = st.selectbox(
        "Estação Metereológica",
        meteorological_data["Cidade"].unique(),
        3
    )

    actual_date = meteorological_data["Data"].max()
    total_period_in_days = (actual_date - meteorological_data["Data"].min()).days

    count = days_without_rain(selected_station ,actual_date,total_period_in_days)

    tabs = st.tabs(['Chuva', 'Temperatura'])

    with tabs[0]:
        ct = st.container()
        col1, col2 = ct.columns([1, 1], gap='large')
        with col1:
            st.write(f'# Dias sem chuva: {count}')
        with col2:
            chart_view = st.selectbox(
                "Dados",
                ["Diário","Hora em Hora"],
                0
            )

        if(count < total_period_in_days):

            last_day_with_rain = meteorological_data["Data"].max() - datetime.timedelta(days=count)
            precipitation_last_day_with_rain = (meteorological_data[(meteorological_data["Data"] == last_day_with_rain) & (meteorological_data["Cidade"] == selected_station)])["Precipitação Acumulada"].sum()

            rain_classification = ""

            if precipitation_last_day_with_rain <= 10:
                rain_classification = "Fraca <= 10mm"
            elif precipitation_last_day_with_rain > 10 and precipitation_last_day_with_rain <= 50:
                rain_classification = "Moderada - entre 10mm e 50mm"
            elif precipitation_last_day_with_rain > 50:
                rain_classification = "Forte > 50mm"

            st.write(f'Último dia com chuva: {last_day_with_rain.date()}: {precipitation_last_day_with_rain:.2f}mm - Classificação: {rain_classification}')
        else:
            st.write(f'Último dia com chuva: Sem dados!')

        if chart_view == "Diário":
            createDailyChart(selected_station, metrics[4])
        elif chart_view == "Hora em Hora":
            createHourlyChart(selected_station, metrics[4])

    with tabs[1]:
        createHourlyChart(selected_station, metrics[4])



with col_map_city:
    pr_map = create_city_map(selected_station)
    folium_static(pr_map)

secao2 = st.container()
col_precipitation_map, col_2 = secao2.columns([1, 1], gap='large')
col_precipitation_map.header('Precipitação acumulada no Período')

meteorological_data = load_csv_data("data/dados_meteorologicos_simepar_parana.csv")
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

    folium_static(pr_map)