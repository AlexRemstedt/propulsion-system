import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import scipy.interpolate as interp

# Variabelen
diameter = 0.1333  # [m]
spoedverhouding = 0.66  # Dimensieloos
rho_water = 999  # [kg/m^3]
toerental = 15   # [1/s]
coef_a = -0.3821  # Dimensieloos
coef_b = 0.2885  # Dimensieloos
coef_c = -0.03346  # Dimensieloos
coef_d = 0.0308  # Dimensieloos

# Basisdata arrays
snelheden = np.arange(0.0, 1.8, 0.2)
stuwkracht = np.array([21, 19, 16.6, 13.9, 10.9, 7.83, 4.57, 0.88, -3.10])
koppel = np.array([0.282, 0.262, 0.240, 0.211, 0.178, 0.144, 0.105, 0.063, 0.012])
snelheidsgraden = np.zeros(len(snelheden))
oude_snelheidsgraden = np.arange(0, 1.0, 0.01)
k_t_oud = np.zeros(len(oude_snelheidsgraden))
k_q_oud = np.zeros(len(oude_snelheidsgraden))

# Bepaal een functie waarmee data wordt geinterpoleerd


def interpolator(interpolating_array, reference_length):
    """Interpoleert data over een gegeven reference array.
        Interpolating_array: numpy array
            Basisdata die over de reference gestretched moet worden.
        Reference_length: getal
            Tot welke waarde de basisdata gestretched moet worden.
    """
    reference = np.linspace(0, reference_length, 1000)

    array_interp = interp.interp1d(np.arange(interpolating_array.size), interpolating_array, kind='quadratic')
    interpolated_data = array_interp(np.linspace(0, interpolating_array.size - 1, reference.size))

    return interpolated_data


# Bereken de nieuwe geinterpoleerde waardes
nieuwe_stuwkracht = interpolator(stuwkracht, np.amax(stuwkracht))
nieuwe_koppel = interpolator(koppel, np.amax(koppel))


# Bereken de snelheidsgraad J
for n in range(len(snelheden)):
    snelheidsgraden[n] = snelheden[n] / (toerental * diameter)

nieuwe_snelheidsgraden = interpolator(snelheidsgraden, np.amax(snelheidsgraden))


# Bereken de waardes van K_T
k_t = np.zeros(len(nieuwe_stuwkracht))

for n in range(len(nieuwe_stuwkracht)):
    k_t[n] = nieuwe_stuwkracht[n] / (rho_water * (diameter ** 4) * (toerental ** 2))


# Bereken de waardes van K_Q
k_q = np.zeros(len(nieuwe_koppel))

for n in range(len(nieuwe_koppel)):
    k_q[n] = nieuwe_koppel[n] / (rho_water * (diameter ** 5) * (toerental ** 2))


# Bereken de waardes van Eta_0
eta_0 = np.zeros(len(k_t))

for n in range(len(nieuwe_snelheidsgraden)):
    eta_0[n] = k_t[n] * nieuwe_snelheidsgraden[n] / (k_q[n] * 2 * np.pi)


# Fit de data arrays om de coefficienten te bepalen
print("Voor K_T bestaan de waardes c, b, en a als", poly.polyfit(nieuwe_snelheidsgraden, k_t, 2))
print("Voor K_Q bestaan de waardes f, e, en d als", poly.polyfit(nieuwe_snelheidsgraden, k_q, 2))


# Genereer de waardes voor de oude variabelen
for n in range(len(oude_snelheidsgraden)):
    k_t_oud[n] = coef_a * oude_snelheidsgraden[n] + coef_b
    k_q_oud[n] = coef_c * oude_snelheidsgraden[n] + coef_d


# Plot de gegevens
plt.plot(nieuwe_snelheidsgraden, k_t, 'r', label="K_T, nieuw")
plt.plot(nieuwe_snelheidsgraden, k_q * 10, 'm', label="10 * K_Q, nieuw")
plt.plot(nieuwe_snelheidsgraden, eta_0, 'c', label="Open water rendement, nieuw")
plt.plot(oude_snelheidsgraden, k_t_oud, 'k--', label="K_T, oud")
plt.plot(oude_snelheidsgraden, k_q_oud * 10, 'k:', label="10 * K_Q, oud")
plt.legend(loc="best", shadow=True)
plt.ylabel("Waardes KT, KQ, Eta_0")
plt.xlabel("J")
plt.xlim(0, 0.8)
plt.ylim(0, 0.8)
plt.grid()
plt.show()
