import altair as alt
import pandas as pd
import streamlit as st
import datetime
import numpy as np

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



metricas = ["Temperatura Média","Umidade Relativa","Velocidade do Vento","Rajada do Vento","Precipitação Acumulada","Pressão Atmosférica Reduzida"]


df = load_data()
df.fillna(0, inplace=True)
df["Data"] = pd.to_datetime(df["Data"])

# Show a multiselect widget with the genres using `st.multiselect`.
cidades = st.multiselect(
    "Estações Metereológicas",
    df.Cidade.unique(),
    ["Apucarana","Palmas","Guaratuba"]
)

# Select the metrics to analyse
metrica = st.selectbox(
    "Métrica",
    metricas,
    index=0
)


st.write(f'Registros na base de dados: {df["Data"].min().date()} até {df["Data"].max().date()}')
#dates = st.date_input("Selecione um período",(datetime.date(2024,7,1),datetime.date(2024,7,3)),min_value=datetime.date(2024,6,30),max_value=datetime.date(2024,7,3))
data = st.date_input("Selecione uma data",df["Data"].max()-datetime.timedelta(days=1),min_value=df["Data"].min(),max_value=df["Data"].max())


data = pd.to_datetime(data)


#Filter the dataframe based on the widget input and reshape it.
df_filtered = df[(df["Cidade"].isin(cidades)) & (df["Data"] == data)]


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
    index=0
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
