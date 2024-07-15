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
page_title = "Dados Meteorológicos Simepar Paraná"
page_icon = "☁️"

st.set_page_config(page_title=page_title, page_icon=page_icon)
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

selected_station = st.selectbox(
    "Estação Metereológica",
    meteorological_data["Cidade"].unique(),
    3
)

# Select the metrics to analyse
selected_metric = st.selectbox(
    "Métrica",
    metrics,
    index=4
)

st.write(f'Registros na base de dados: {meteorological_data["Data"].min().date()} até {meteorological_data["Data"].max().date()}')
date = st.date_input("Selecione uma data",meteorological_data["Data"].max()-datetime.timedelta(days=1),min_value=meteorological_data["Data"].min(),max_value=meteorological_data["Data"].max())


#Filter the dataframe based on the widget input and reshape it.
df_filtered = meteorological_data[(meteorological_data["Cidade"] == selected_station) & (meteorological_data["Data"] == pd.to_datetime(date))]


df_reshaped = df_filtered.pivot_table(
    index="Horario", columns="Cidade", values=selected_metric, fill_value=0
)

# Display the data as an Altair chart using `st.altair_chart`.
df_chart = pd.melt(
    df_reshaped.reset_index(), id_vars="Horario", var_name="Cidade", value_name=selected_metric
)


chart = (
    alt.Chart(df_chart)
    .mark_line()
    .encode(
        x=alt.X("Horario:N", title="Horário"),
        y=alt.Y(f'{selected_metric}:Q', title=f'{selected_metric} {metrics_unit[selected_metric]}'),
        color="Cidade:N"
    )
    .properties(height=320)
)
st.altair_chart(chart, use_container_width=True)


#Days without rain Section
city_option = st.selectbox(
    "Selecione uma Estação",
    meteorological_data.Cidade.unique(),
    index=3
)

#calculating days without rain
def days_without_rain(start_date,total_period_in_days):

    count = 0
    date = start_date
    for i in range(0,total_period_in_days):
        day = meteorological_data[(meteorological_data["Data"] == date) & (meteorological_data["Cidade"] == city_option)]
        max_precipitation = day["Precipitação Acumulada"].sum()
        if(max_precipitation <= 0):
            count += 1
        else:
            return count
        date -= datetime.timedelta(days=1)

start_date = meteorological_data["Data"].max()
total_period_in_days = (start_date - meteorological_data["Data"].min()).days

count = days_without_rain(start_date,total_period_in_days)

st.write(f'# Dias sem chuva: {count}')


if(count < total_period_in_days):

    last_day_with_rain = meteorological_data["Data"].max() - datetime.timedelta(days=count)
    precipitation_last_day_with_rain = (meteorological_data[(meteorological_data["Data"] == last_day_with_rain) & (meteorological_data["Cidade"] == city_option)])["Precipitação Acumulada"].sum()

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


#Show the City Map
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

pr_map = create_city_map(city_option)
folium_static(pr_map)


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

# Estimate precipitation on all cities 
cities_geojson_data = load_geo_json_data("./data/geojs-41-mun.json")

cities_lat_long = pd.read_csv("./data/lat_long_cidades_pr.csv",sep=",")

df_period_total = meteorological_data.groupby(['Cidade']).agg({'Precipitação Acumulada': 'sum'}).reset_index()

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


st.title('Mapa de Precipitação Acumulada no Período')
folium_static(pr_map)