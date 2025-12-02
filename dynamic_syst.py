import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Definir el sistema de ecuaciones diferenciales
def odes(y, t):
    X_H, X_T, S1, S2, E, P, V = y
    X = X_H + X_T  # Concentración total de biomasa [g/L]

    Q = 0 #Caudal inicial [L/h]
    if S1 <= 0.01:
      Q = 0 #(600*0.021)*np.exp(0.021*t) #Caudal [L/h]
    dVdt = Q # Volumen [L]
    #Balance de glucosa
    S1_in = 0 # Concentracion de sustrato de entrada [g/L]
    Y_S1 = 0.46188 # c []
    K_S1 = 0.1 # Constante de saturación (glucosa) [g/L]
    I1 = 20 # Constante de represión de diferenciación morfológica (glucosa) [g/L]
    u_maxS1 = 0.042 # Tasa específica de crecimiento máxima (glucosa) []
    u_S1 = u_maxS1*S1/(K_S1+S1) # Tasa específica de crecimiento (glucosa) []
    dS1dt = -(u_S1/Y_S1)*X + (S1_in - S1)*(Q/V) #Balance de sustrato (glucosa)
    #Balance de Sacarosa
    S2_in = 30 # Concentracion de sustrato de entrada [g/L]
    m_S2 = 0.02267 # Mantenimiento (sacarosa) [1/h]
    Y_S2 = 0.4 # Factor de rendimiento (sacarosa)
    K_S2 = 10.199 # Constante de saturación (sacarosa)
    I2 = 300 # Constante de represión de uso de sacarosa por glucosa
    u_maxS2 = 0.021 # Tasa específica de crecimiento máxima (sacarosa)
    phi = 0
    if S1 <= 0.01:
      phi = 1
    u_S2 = phi*(u_maxS2*S2/(K_S2+S2)) # Tasa específica de crecimiento (sacarosa) []
    dS2dt = -((u_S2/Y_S2)+(m_S2/(K_S2+S2+I2*S1)))*X + (S2_in - S2)*(Q/V) #Balance de sustrato (sacarosa)
    #Balance para las células
    X_H_in = 0 # Concentracion de hifas de entrada [g/L]
    X_T_in = 0 # Concentracion de artrosporas de entrada [g/L]
    d_Hmax = 0.00668 # Tasa específica de muerte máxima para hifas []
    d_H = d_Hmax/(1+I1*S1) # Tasa específica de muerte para hifas []
    u_Tmax = 0.04526 # Tasa específica de crecimiento máxima para artrosporas []
    u_T = u_Tmax/(1+I1*S1) # Tasa específica de crecimiento para artrosporas []
    d_T = 0.00441 # Tasa específica de muerte para artrosporas []
    u_H = u_S1 # Tasa específica de crecimiento para hifas []
    if S1 <= 0.01:
      u_H = u_S2 # Tasa específica de crecimiento para hifas []
    dX_Hdt = (u_H*X - u_T*X_H - d_H*X_H) + (X_H_in - X_H)*(Q/V) #Balance para hifas
    dX_Tdt = (u_T*X_H - d_T*X_T) + (X_T_in - X_T)*(Q/V) #Balance para artrosporas

    a = 9 # Tasa de formación de enzima ligada a crecimiento []
    b = 2.18787 # Tasa de degradacíon de enzima ligada a crecimiento []
    dEdt = (a*u_T*X_H - b*E) - E*(Q/V) #Balance de enzima

    gamma = 0.01076 # Tasa de degradación de CPC []
    dPdt = (E*X_T - gamma*P) - P*(Q/V) #Balance de CPC

    return [dX_Hdt, dX_Tdt, dS1dt, dS2dt, dEdt, dPdt, dVdt]

# Condiciones iniciales
y0 = [5, 0, 40, 5, 0, 0, 12000] # Concentraciones iniciales (hifas, artrosporas, glucosa, sacarosa, enzima, CPC) [g/L]

# Definir el span de tiempo
t = np.linspace(0, 100, 100) # [horas]

# Solucionar el sistema
sol = odeint(odes, y0, t)

# Extraer las soluciones de las EDOs
X_H = sol[:, 0]
X_T = sol[:, 1]
S1 = sol[:, 2]
S2 = sol[:, 3]
E = sol[:, 4]
P = sol[:, 5]
V = sol[:, 6]

# Plots
fig, ax = plt.subplots()
fig.subplots_adjust(right=1)

twin1 = ax.twinx()
twin2 = ax.twinx()
twin3 = ax.twinx()

# Para colocar a la derecha
twin2.spines.right.set_position(("axes", 1.12))
twin3.spines.right.set_position(("axes", 1.24))
#Colores de cada subset
p1, = ax.plot(t, X_H + X_T, color='#548235', label='$X$')
p2, = twin1.plot(t, S1, color='#ff8000', label='$S_1$')
p3, = twin1.plot(t, S2, color='#fbd415', label='$S_2$')
p4, = twin2.plot(t, P, 'r', label='$P$')
p5, = twin3.plot(t, E, label='$E$')
#límites de cada eje
ax.set_xlim(0, t[-1])
ax.set_ylim(0, 40)
twin1.set_ylim(0, 40)
twin2.set_ylim(0, max(P)+100)
twin3.set_ylim(0, 10)
#Labels por eje
ax.set_xlabel("Tiempo $t$ [h]")
ax.set_ylabel("Biomasa [g/L]")
twin1.set_ylabel("Sustrato [g/L]")
twin2.set_ylabel("Producto [mg/L]")
twin3.set_ylabel("Enzima [mg/g.hr]")

ax.grid(True)
plt.minorticks_on()
ax.grid(True, which='minor', axis='both', linestyle='-', alpha=0.3)

ax.legend(handles=[p1, p2, p3, p4, p5], loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=5)

# Apply tight layout
plt.tight_layout()
plt.savefig("batch.png",dpi=600)
plt.show()


print('\nEl valor final de CPC es:',round(P[-1],3),'mg/L')
print('\nEl valor maximo de CPC es:',round(max(P),3),'mg/L')
print('\nEl valor total de CPC es:',round((P[-1]*V[-1])/1000000,3),'kg \n')



