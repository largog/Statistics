import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#  definir 4-parameter log-logistic model
def loglogistic4(x, E0, Emax, EC50, h):
    return E0 + Emax * ((x**h) / ( x**h + (EC50**h)))

# Leer el archivo
df = pd.read_excel("../data/model.xlsx")

# Extraer la mediana para cada concentración

x = df["PQS"].to_numpy()
y = df.iloc[:, 1:].median(axis=1).to_numpy()

#Ajustar al modelo
# condiciones iniciales para el ajuste
p0 = [np.min(y), np.max(y)-np.min(y), np.median(x), 1.0]

params, cov = curve_fit(loglogistic4, x, y,
                        p0=p0)

E0, Emax, EC50, h = params
pred_4p = loglogistic4(x, E0, Emax, EC50, h)


#  Gráfico
sns.scatterplot(x=x, y=y, s=60, color="black")
sns.lineplot(x=x, y=pred_4p)
sns.despine()
plt.xlabel("Concentration")
plt.ylabel("Viability (%)")
plt.title("CCK8-assay 4-parameter log-logistic fit")
plt.savefig("Model1.png",dpi=600)
plt.show()


# Parámetros
print("Fit parameters:")
print(f"  E0     = {E0:.3f}")
print(f"  Emax   = {Emax:.3f}")
print(f"  EC50   = {EC50:.3f}")
print(f"  h      = {h:.3f}")


# Generar el rango de graficación en escala logarítmica
x_smooth = np.logspace(np.log10(min(x[x>0])), np.log10(max(x)), 300)

# Pasar los valores del modelo al rango logarítmico
y_fit = loglogistic4(x_smooth, *params)
# Datos reales
plt.scatter(x, y, label="Data", color="red")
# Curva del modelo
plt.plot(x_smooth, y_fit, label="4p LL fit", color="red")
# Log scale on x-axis
plt.xscale("log")


plt.axvline(EC50, color="gray", linestyle="--", linewidth=1)
plt.text(
    EC50,
    np.interp(EC50, x_smooth, y_fit),     # y-position on the curve
    f" EC₅₀ = {EC50:.2f} µg/mL",
    fontsize=10,
    verticalalignment="bottom",
    horizontalalignment="left"
)

# Optional formatting
plt.xlabel("PQS concentration (µg/mL, log scale)")
plt.ylabel("Cell viability (%)")
plt.legend()
sns.despine()
plt.savefig("Model2_log.png",dpi=600)
plt.show()