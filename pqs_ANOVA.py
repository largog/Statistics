import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from statannotations.Annotator import Annotator

# Leer el archivo
df = pd.read_excel("../data/ow_anova.xlsx")

# Extraer la fila con los valores
values = df.iloc[0]

# Construir el long a partir de los nombres de columnas
df_long = pd.DataFrame({
    "strain": values.index,
    "PQS_conc": values.values
})

# Limpieza
df_long["strain"] = df_long["strain"].str.replace(r"\.\d+$", "", regex=True)

#Eliminar la primera fila de títulos

df_long = df_long.drop(index=0)
df_long["PQS_conc"] = pd.to_numeric(df_long["PQS_conc"], errors="coerce")
# Agrupar valores por strain
groups = [g["PQS_conc"].astype(float).values
          for name, g in df_long.groupby("strain")]

# ANOVA
F, p = f_oneway(*groups)

print(f"El valor de F es: {F}")
print(f"El valor de p es: {p}")

print(df_long)

# Tukey HSD
tukey = pairwise_tukeyhsd(
    endog=df_long["PQS_conc"],     # dependent variable
    groups=df_long["strain"],      # strain categories
    alpha=0.05
)

# Crear lista de comparaciones significativas
pairs = []
pvalues = []

for res in tukey.summary().data[1:]:
    g1, g2, meandiff, p_adj, lower, upper, reject = res
    pairs.append((g1, g2))
    pvalues.append(p_adj)


print("\nTukey HSD results:")
print(tukey)

plt.figure(figsize=(5,5))

# Paleta con un color por categoría
unique_strains = df_long["strain"].unique()
palette = ["#f04442", "#1afe1d", "#3f44d7"]

ax = sns.barplot(
    data=df_long,
    x="strain",
    y="PQS_conc",
    edgecolor="none",
    linewidth=2
)
# Barra sin relleno
for patch, color in zip(ax.patches, palette):
    patch.set_facecolor('none')
    patch.set_edgecolor(color)
    patch.set_linewidth(2)
# Puntos individuales
sns.swarmplot(
    data=df_long,
    x="strain",
    y="PQS_conc",
    size=6,
    palette=palette,
)

# Añadir anotaciones de significancia
annot = Annotator(
    ax,
    pairs,
    data=df_long,
    x="strain",
    y="PQS_conc"
)

annot.configure(
    test=None,
    text_format="star",
    loc="outside",
    verbose=False
)

annot.set_pvalues(pvalues)
annot.annotate()

sns.despine()
plt.xlabel("Strain")
plt.ylabel("PQS concentration (μg/mL)")

plt.rcParams["font.family"] = "Helvetica"
plt.tight_layout()
#Save the plot
plt.savefig("ANOVA_tukey.png",dpi=600)
plt.show()
