import numpy as np
import math
import matplotlib.pyplot as plt

# Standaard Variabelen
alpha_l = 19  # Dimensieloos
kin_vis = 10 ** -6  # m^2 / s
lengte_model = 31.5 / 19  # m
lengte_waterlijn = 31.5  # m
g = 9.81  # m / s^2
rho_m = 998.8  # kg / m^3
rho_s = rho_m * 1.0282  # kg / m^3
oppervlak_m = 0.755124654  # m^2
oppervlak_s = 272.6  # m^2

# Data van het modelschip
snelheid = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
weerstand = np.array([0.110, 0.220, 0.410, 0.620, 0.870, 1.180, 1.570, 2.060, 2.620, 3.480, 4.790,
                      6.220, 8.060, 11.340, 16.650])

# Data van de kotter
snelheid_s = np.arange(0.868125, 7.3790625, 0.4340625)
weerstand_s = np.zeros(len(snelheid_s))

# Arrays
"""Arrays van het model"""
froude_array = np.zeros(len(snelheid))
re_array = np.zeros(len(snelheid))
cf_array = np.zeros(len(snelheid))
ct_array = np.zeros(len(snelheid))
fr4_cf_array = np.zeros(len(snelheid))
ct_cf_array = np.zeros(len(snelheid))
cr_array = np.zeros(len(snelheid))

"""Arrays van het schip"""
re_s_array = np.zeros(len(snelheid_s))
cf_s_array = np.zeros(len(snelheid_s))
ct_s_array = np.zeros(len(snelheid_s))


# BEREKENINGEN
for n in range(len(snelheid)):
    # Bereken het Froudegetal
    froude_array[n] = snelheid[n] / (g * lengte_waterlijn / alpha_l) ** 0.5

    # Bereken het Reynold's getal
    re_array[n] = snelheid[n] * (lengte_waterlijn / alpha_l) / kin_vis

    # Bereken cf
    if re_array[n] == 0:
        cf_array[n] = 0
    else:
        cf_array[n] = 0.075 / ((math.log(re_array[n], 10) - 2) ** 2)

    # Bereken ct
    if snelheid[n] == 0:
        ct_array[n] = 0
    else:
        ct_array[n] = weerstand[n] / (0.5 * rho_m * snelheid[n] * snelheid[n] * oppervlak_m)

    # Bereken fr4_cf
    if cf_array[n] == 0:
        fr4_cf_array[n] = 0
    else:
        fr4_cf_array[n] = froude_array[n] ** 4 / cf_array[n]

    # Bereken ct_cf
    if cf_array[n] == 0:
        ct_cf_array[n] = 0
    else:
        ct_cf_array[n] = ct_array[n] / cf_array[n]

    # Bereken cr
    cr_array[n] = ct_array[n] - cf_array[n]

    # Bereken het Reynold's getal van de kotter zelve
    re_s_array[n] = snelheid_s[n] * lengte_waterlijn / kin_vis

    # Bereken cf van de kotter
    if re_s_array[n] == 0:
        cf_s_array[n] = 0
    else:
        cf_s_array[n] = 0.075 / ((math.log(re_s_array[n], 10) - 2) ** 2)

    # Bereken ct van de kotter
    ct_s_array[n] = cf_s_array[n] + cr_array[n]

    # Bereken de algehele weerstand van de kotter
    weerstand_s[n] = ct_s_array[n] * 0.5 * rho_s * (snelheid_s[n] ** 2) * oppervlak_s


# plot the data itself
plt.plot(fr4_cf_array, ct_cf_array, 'o', label='plot punten')

# calc the trendline
z = np.polyfit(fr4_cf_array, ct_cf_array, 1)
p = np.poly1d(z)
plt.plot(fr4_cf_array, p(fr4_cf_array), "r--", label='regressielijn')
plt.legend(loc = 'lower center', shadow = True, fontsize = 'large')
plt.xlim(0, 7)
plt.ylim(0, 5)
plt.title("Prohaska Plot, MT05")
plt.xlabel("Fr^4 / Cf")
plt.ylabel("Ct / Cf")
plt.grid()
plt.show()







