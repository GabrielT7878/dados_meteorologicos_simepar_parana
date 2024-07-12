import altair as alt
import pandas as pd
import streamlit as st
import datetime
import numpy as np

import folium
from streamlit_folium import folium_static

from scipy.spatial import cKDTree

import json

# Show the page title and description.
st.set_page_config(page_title="Dados Meteorológicos Simepar Paraná", page_icon="☁️")
st.title("Dados Meteorológicos Simepar Paraná")
st.write(
    """
    Este aplicativo visualiza dados do [SIMEPAR](http://www.simepar.br/prognozweb/simepar).
    Ele mostra métricas sobre o clima coletadas em estações meteorológicas no Paraná-RS.
    """
)


# Load the data from a CSV. We're caching this so it doesn't reload every time the app
# reruns (e.g. if the user interacts with the widgets).
#@st.cache_data
def load_data():
    df = pd.read_csv("data/dados_meteorologicos_simepar_parana.csv")

    return df

#carregando coordenadas das estações
with open('data/coordenadas_cidades.json', 'r') as arquivo:
    estacoes = json.load(arquivo)   



metricas = ["Temperatura Média","Umidade Relativa","Velocidade do Vento","Rajada do Vento","Precipitação Acumulada","Pressão Atmosférica Reduzida"]


df = load_data()
df.fillna(0, inplace=True)
df["Data"] = pd.to_datetime(df["Data"])

# Show a multiselect widget with the genres using `st.multiselect`.
# cidade = st.multiselect(
#     "Estações Metereológicas",
#     df.Cidade.unique(),
#     ["Apucarana","Palmas","Guaratuba"]
# )

cidade = st.selectbox(
    "Estação Metereológica",
    df["Cidade"].unique(),
    3
)



# Select the metrics to analyse
metrica = st.selectbox(
    "Métrica",
    metricas,
    index=4
)


st.write(f'Registros na base de dados: {df["Data"].min().date()} até {df["Data"].max().date()}')
#dates = st.date_input("Selecione um período",(datetime.date(2024,7,1),datetime.date(2024,7,3)),min_value=datetime.date(2024,6,30),max_value=datetime.date(2024,7,3))
data = st.date_input("Selecione uma data",df["Data"].max()-datetime.timedelta(days=1),min_value=df["Data"].min(),max_value=df["Data"].max())


data = pd.to_datetime(data)


#Filter the dataframe based on the widget input and reshape it.
df_filtered = df[(df["Cidade"] == cidade) & (df["Data"] == data)]


df_reshaped = df_filtered.pivot_table(
    index="Horario", columns="Cidade", values=metrica, fill_value=0
)

print(df_reshaped)

#
#Display the data as a table using `st.dataframe`.
# st.dataframe(
#     df_reshaped,
#     use_container_width=True,
#     column_config={"Horario": st.column_config.TextColumn("Horário")},
# )

# Display the data as an Altair chart using `st.altair_chart`.
df_chart = pd.melt(
    df_reshaped.reset_index(), id_vars="Horario", var_name="Cidade", value_name=metrica
)

unidadeMetrica = {
    "Temperatura Média": "(ºC)",
    "Umidade Relativa": "(%)",
    "Velocidade do Vento": "(Km/h)",
    "Rajada do Vento": "(Km/h)",
    "Precipitação Acumulada": "(mm)",
    "Pressão Atmosférica Reduzida": "(hPa)"
}

chart = (
    alt.Chart(df_chart)
    .mark_line()
    .encode(
        x=alt.X("Horario:N", title="Horário"),
        y=alt.Y(f'{metrica}:Q', title=f'{metrica} {unidadeMetrica[metrica]}'),
        color="Cidade:N"
    )
    .properties(height=320)
)
st.altair_chart(chart, use_container_width=True)


optionChuva = st.selectbox(
    "Selecione uma Estação",
    df.Cidade.unique(),
    index=3
)

#Calculando os dias sem chuva
data = df["Data"].max()
periodo = data - df["Data"].min()
count = 0

print(data)
print()


for i in range(0,periodo.days):
    dia = df[(df["Data"] == data) & (df["Cidade"] == optionChuva)]
    precipitacaoMax = dia["Precipitação Acumulada"].max()
    if(precipitacaoMax <= 0):
        count += 1
    else:
        break
    data += datetime.timedelta(days=-1)

st.write(f'# Dias sem chuva: {count}')


classificacaoChuva = ""

if precipitacaoMax < 2.5:
    classificacaoChuva = "Fraca < 2.5mm"
elif precipitacaoMax < 10:
    classificacaoChuva = "Moderada < 10mm"
elif precipitacaoMax > 50:
    classificacaoChuva = "Forte > 50mm"

if(precipitacaoMax != 0 and precipitacaoMax != np.nan):
    st.write(f'Último dia com chuva: {data.date()}: {precipitacaoMax}mm - Classificação: {classificacaoChuva}')
else:
    st.write(f'Último dia com chuva: Sem dados!')

def get_lat_long(city):
    location = estacoes[city]
    if location:
        return location["latitude"], location["longitude"]
    else:
        return None, None

latitude, longitude = get_lat_long(optionChuva)


m = folium.Map(location=[latitude, longitude], zoom_start=12, min_zoom=6, max_zoom=13)

# Adicionando um marcador para a estação selecionada
folium.Marker([latitude, longitude], tooltip=optionChuva).add_to(m)

# Exibindo o mapa no Streamlit
folium_static(m)


#Teste Mapa de calor

import folium
import branca.colormap
import numpy as np
from folium.plugins import HeatMap
from folium import Choropleth

data = "./data/geojs-41-mun.json"
df_merged2 = pd.read_csv("conversaoDados/heat_map.csv",sep=",")
# df_merged = pd.read_csv("./data/Interpolated_Precipitation_Data.csv",sep=",")

df_base = pd.read_csv("./data/intepolation_base_cidades.csv",sep=",")

df_period_total = df.groupby(['Cidade']).agg({'Precipitação Acumulada': 'sum'}).reset_index()

df_period_total = df_period_total[["Cidade","Precipitação Acumulada"]]

df_merged = df_base.merge(df_period_total,on="Cidade",how="outer")

print(len(df_merged))

def interpolateData(data=None):
    # Separar os dados com precipitação (não NaN) e sem precipitação (NaN)
    known_precipitation = data.dropna(subset=['Precipitação Acumulada'])
    unknown_precipitation = data[data['Precipitação Acumulada'].isna()]

    # Construir uma árvore KD para os pontos conhecidos
    tree = cKDTree(known_precipitation[['latitude', 'longitude']])

    # Definir uma função de inverso da distância
    def inverse_distance_weighting(distances, values, num_nearest=3):
        inverse_distances = 1 / distances
        weights = inverse_distances / inverse_distances.sum()
        return np.dot(weights, values)

    # Encontrar os pontos mais próximos e calcular a interpolação
    def interpolate_missing_values(row):
        distances, indices = tree.query([row['latitude'], row['longitude']], k=3)
        nearest_values = known_precipitation.iloc[indices]['Precipitação Acumulada'].values
        return inverse_distance_weighting(distances, nearest_values)

    # Aplicar a interpolação para cada linha com precipitação desconhecida
    interpolated_values = unknown_precipitation.apply(interpolate_missing_values, axis=1)

    # Atualizar o DataFrame original com os valores interpolados
    data.loc[data['Precipitação Acumulada'].isna(), 'Precipitação Acumulada'] = interpolated_values

    return data


df_merged = interpolateData(df_merged)


with open(data) as f:
    geojson_data = json.load(f)

m = folium.Map(location=[-24.5452, -51.5652], zoom_start=7)

estilo = lambda x: {"color" : "white",
                   "fillOpacity": 1,
                   "weight": 1}

# Adiciona os limites municipais ao mapa
folium.GeoJson(
    geojson_data,
    name='geojson',
    style_function=estilo
).add_to(m)


Choropleth(
    geo_data=geojson_data,
    name='choropleth',
    data=df_merged,
    columns=['Cidade', 'Precipitação Acumulada'],
    key_on='feature.properties.name',
    fill_color='YlOrRd',
    fill_opacity=0.7,
    line_opacity=0.2,
    legend_name='Precipitação Acumulada (mm)',
    nan_fill_opacity=0.0,
    nan_fill_color="white"
).add_to(m)


#Adicionar as cidades ao mapa com a precipitação acumulada
# for _, row in df_merged2.iterrows():
#     folium.CircleMarker(
#         location=[row['latitude'], row['longitude']],
#         radius=0,
#         popup=folium.Popup(f"{row['Cidade']}: {row['Precipitação Acumulada']} mm", parse_html=True),
#         tooltip=folium.Tooltip(f"{row['Cidade']}: {row['Precipitação Acumulada']} mm"),
#         color='orange',
#         fill=True,
#         fill_color='orange',
#     ).add_to(m)

# Exibir o mapa no Streamlit
st.title('Mapa de Precipitação Acumulada no Período')
folium_static(m)