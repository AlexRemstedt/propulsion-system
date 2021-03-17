"""
Alex Remstedt
17/03/2021

De opdrachten uit de handleiding op pagina 14 worden in dit bestand berekent.
Er wordt gerefereerd naar een figuur, dit figuur vind men terug als figuur 3 uit de handleiding of figuur 56 uit het
dictaat van dr.ir. P. de Jong (2014)
"""
# imports
import numpy as np

# variables
dichtheid = 1024  # kg/m^3
snelheid_kn = np.array([3, 11.5])  # kn
ship_speed = snelheid_kn * .514  # m/s
volgstroomgetal = .25  # (geen eenheid)
toerental = np.array([8/6, 14/6])  # s^-1
schroefdiameter = 2.53  # m
K_T = np.array([.32, .16])  # Stuwkracht coëfficient
K_Q = np.array([.043, .027])  # Askoppel coëfficient


# functions
def instroomsnelheid(v_s):
    """
    Bereken de instroomsnelheid rondom de schroef

    :param v_s: ship_speed
    :type v_s: Union(np.ndarray, float)
    :return: v_0
    """
    v_0 = v_s * (1 - volgstroomgetal)
    return v_0


def snelheidsgraad(v_s):
    """
    Bereken de snelheidsgraad

    :param v_s: ship speed
    :type v_s: Union(np.ndarray, float)
    :return: J
    """
    j = np.zeros(len(v_s))
    for k in range(len(j)):
        j[k] = instroomsnelheid(v_s)[k] / toerental[k] / schroefdiameter
    return j


def open_water_rendement(v_s):
    """
    Bereken open water rendement

    :param v_s:
    :return:
    """
    eta_0 = np.zeros(len(v_s))
    for k in range(len(eta_0)):
        eta_0[k] = K_T[k] * snelheidsgraad(v_s)[k] / K_Q[k] / 2 / np.pi
    return eta_0


def thrust():
    """
    Bereken de voortstuwende kracht

    :return: thrust
    """
    t = K_T * dichtheid * toerental ** 2 * schroefdiameter ** 4
    return t


def as_koppel():
    """
    Bereken het koppel in de as

    :return: askoppel
    """
    q = K_Q * dichtheid * toerental ** 2 * schroefdiameter ** 5
    return q


# prints & plots
if __name__ == '__main__':
    print(f"""
    a)  Bepaal de instroomsnelheid ter plekke van de schroef v_0 voor beide condities: {instroomsnelheid(ship_speed)}
    b)  Bepaal de snelheidsgraad J en vervolgens het openwater rendement en de K_T en de K_Q 
        van de schroef in beide condities:
            J = {snelheidsgraad(ship_speed)}
            K_T = {K_T}
            K_Q = {K_Q}
            open water rendement = {open_water_rendement(ship_speed)}
    c)  Bepaal de stuwkracht T en het benodigde askopppel Q in beide condities:
            T = {thrust()}
            Q = {as_koppel()}
    d)  Hoeveel kan het open water rendement verbeterd worden in de vissende conditie door gebruik te maken 
        van een verstelbare schroef (variabele spoed) en hoe groot is dan de optimale spoedverhouding voor de 
        vissende conditie? J, K_T en K_Q blijven gelijk.
        
        Uit de tabel blijkt dat een P/D van {.6} het optimaalste is voor vissende condities.
    """)
