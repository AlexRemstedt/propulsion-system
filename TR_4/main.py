"""
Alex Remstedt

Een aanvulling op "TR_3/main.py" van Kevin en JRM en EU
De verliezen van de motor
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import time
import scipy.interpolate as interp
from scipy import integrate

# ----------- parameters for simulation --------------------------------------

tmax = 36000  # simulation time [s]
dt = 1  # timestep [s]

# fuel properties
LHV = 42700  # Lower Heating Value [kJ/kg]
print('fuel properties loaded')

# water properties
rho_sw = 1025  # density of seawater [kg/m3]
kin_vis = 10 ** -6  # kinematische viscositeitscoefficient [m^2 / s]
print('water properties loaded')

# ship data
m_ship = 358000  # ship mass [kg]
c1 = 1500  # resistance coefficient c1 in R = c1*vs^2
lengte_s = 31.5  # ship length [m]
vormfactor = 0.22  # vormfactor
oppervlakte_s = 272.6  # vrije wateroppervlak van het schip
v_s0 = 6.5430  # ship design speed [m/s]
snelheden = np.linspace(0, v_s0, tmax)  # snelheden in een array
t = 0.1600  # thrust deduction factor[-]
w = 0.2000  # wake factor [-]
print('ship data loaded')

# propellor data
D_p = 3  # diameter of propellor [m]
k_t_a = -0.1709  # value within the formula: a * J^2 + b * J + c
k_t_b = -0.2882  # value within the formula: a * J^2 + b * J + c
k_t_c = 0.2977  # value within the formula: a * J^2 + b * J + c
k_q_d = -0.0204  # value within the formula: d * J^2 + e * J + f
k_q_e = -0.0190  # value within the formula: d * J^2 + e * J + f
k_q_f = 0.0298  # value within the formula: d * J^2 + e * J + f
eta_R = 1.0100  # relative rotative efficiency [-]
print('propellor data loaded')

# engine data
m_f_nom = 1.314762  # nominal fuel injection [g]
n_eng_nom = 900 / 60  # Hz
# eta_e = 0.3800  # nominal engine efficiency [-]
i = 6  # number of cylinders [-]
k_es = 2  # k-factor for engines based on nr.of strokes per cycle
P_b = np.zeros(tmax)  # engine power [kW]
P_b[0] = 960  # Nominal engine power [kW]
M_b = np.zeros(tmax)  # engine torque [Nm]
M_b[0] = P_b[0] * 1000 / 2 / math.pi / n_eng_nom  # ([P_b*1000/2/math.pi/n_eng_nom])
print('engine data loaded')

# gearbox data
eta_TRM = 0.9500  # transmission efficiency [-]
i_gb = 4.2100  # gearbox ratio [-]
I_tot = 200  # total mass of inertia of propulsion system [kg*m^2]
print('gearbox data loaded')

# efficiency parameters
Q_comb_loss = 0
eta_td = np.full(tmax, .52)
Q_comb = 56140  # m_f * h ** L

# initial values
in_p = 3.2830  # initial rpm
iv_t_control = np.array([0, tmax])
X_parms = np.array([.2, 1])  # % maximum fuelrack
Y_parms = np.array([1, 1])  # disturbance factor

# simulation control parameters
xvals = np.linspace(0, tmax - 1, tmax)
ov_X_set = np.interp(xvals, iv_t_control, X_parms)
ov_Y_set = np.interp(xvals, iv_t_control, Y_parms)


# --------- Start van de funtie definities
def cws_berekening():
    cwm = np.array([-0.00418599, -0.00336616, 0.00016877, 0.00100681, 0.00140657, 0.00259207,
                    0.00174022, 0.00312651, 0.00185042, 0.00349679, 0.00328257, 0.00535085,
                    0.00555112, 0.00736961, 0.00933834, 0.01421556])
    reference = np.linspace(0, tmax - 1, tmax)

    array_interp = interp.interp1d(np.arange(cwm.size), cwm, kind='cubic')
    cws_f = array_interp(np.linspace(0, cwm.size - 1, reference.size))

    return cws_f


wrijvingscoef = cws_berekening()


def r_schip(snelheid_schip, wrijvingcoef):
    reynolds = snelheid_schip * lengte_s / kin_vis

    if reynolds == 0:
        cfs = 0
    else:
        cfs = 0.075 / ((math.log(reynolds, 10) - 2) ** 2)

    cts = (1 + vormfactor) * cfs + wrijvingcoef
    global Y
    weerstand = abs(Y * 0.5 * cts * rho_sw * oppervlakte_s * snelheid_schip ** 2)
    return weerstand


def vormweerstand(snelheid_schip):
    reynolds = snelheid_schip * lengte_s / kin_vis

    if reynolds == 0:
        cfs = 0
    else:
        cfs = 0.075 / ((math.log(reynolds, 10) - 2) ** 2)

    cts = vormfactor * cfs
    global Y
    vw = abs(Y * 0.5 * cts * rho_sw * oppervlakte_s * snelheid_schip ** 2)
    return vw  # vormweerstand


def wrijfweerstand(snelheid_schip):
    reynolds = snelheid_schip * lengte_s / kin_vis

    if reynolds == 0:
        cfs = 0
    else:
        cfs = 0.075 / ((math.log(reynolds, 10) - 2) ** 2)

    cts = cfs
    global Y
    weerstand = abs(Y * 0.5 * cts * rho_sw * oppervlakte_s * snelheid_schip ** 2)
    return weerstand


def golfweerstand(snelheid_schip, wrijvingcoef):
    cts = wrijvingcoef
    global Y
    vw = abs(Y * 0.5 * cts * rho_sw * oppervlakte_s * snelheid_schip ** 2)
    return vw  # vormweerstand


def q_cooling_loss(x):
    """
    Calculate energy loss due to cooling

    Parameters
    ----------
    x
        Fuel ratio

    Returns
    -------
        lost cooling warmth
    """
    q = 1908.8 + 7635.2 * x
    return q


def w_mechanical_loss(n):
    """
    calculate lost mechanical work

    Parameters
    ----------
    n
        ratio of rpm

    Returns
    -------
        work
    """
    work = 711.1 + 1659.3 * n / n_eng_nom
    return work


def q_in(x):
    """
    Calculate heat input

    P_i = W_i * n_e / 2
    Q_in = W_i / eta_td

    Parameters
    ----------
    x : float
        brandstofheugelstand

    Returns
    -------
    float
        input heat

    """
    return Q_comb - q_cooling_loss(x)


def combustion_eff():
    return 1


def mech_eff(n, x):
    w_i = q_in(x) * eta_td[0]
    eta = (w_i - w_mechanical_loss(n)) / w_i
    return eta


def heat_eff(x):
    eta = q_in(x) / Q_comb
    return eta


def effective_engine_efficiency(n, x):
    """
    Calculate the efficiency

    Parameters
    ----------
    p : float
        Power
    n : float
        Rpm
    x : float
        Fuel rack

    Returns
    -------
    Effective efficiency
    """
    eta = mech_eff(n, x) * combustion_eff() * heat_eff(x) * eta_td[0]
    return eta


# -------- Make arrays ---------

# Time
mytime = np.linspace(0, tmax - 1, tmax)
# Velocity of the ship [m/s]
v_s = np.zeros(tmax)
v_s[0] = 0
# Distance traveled [m]
s = np.zeros(tmax)
# Advance velocity [m/s]
v_a = np.zeros(tmax)
v_a[0] = 0
# Rpm propellor [Hz]
n_p = np.zeros(tmax)
n_p[0] = in_p
# Rpm diesel engine [Hz]
n_e = np.zeros(tmax)
n_e[0] = 900 / 60  # Nominal engine speed in rotations per second [Hz]
# Resistance [N]
R = np.zeros(tmax)
Y = ov_Y_set[0]
R[0] = r_schip(v_s0, 0)
R_vorm = np.zeros(tmax)
R_vorm[0] = vormweerstand(v_s0)
R_wrijving = np.zeros(tmax)
R_wrijving[0] = 0
R_golf = np.zeros(tmax)
R_golf[0] = 0
# Acceleration ship [m/s^2]
sum_a = np.zeros(tmax)
# Acceleration propellor[1/s^2]
sum_dnpdt = np.zeros(tmax)
m_flux_f = np.zeros(tmax)
out_fc = np.zeros(tmax)

M_Trm = np.zeros(tmax)  # M_B * i_gb * eta_TRM
KT = np.zeros(tmax)  # Thrust coefficient [-]
KQ = np.zeros(tmax)  # Torque coefficient [-]
Rsp = np.zeros(tmax)  # Resistance propelled situation [N]
F_prop = np.zeros(tmax)  # Thrust power [N]
M_prop = np.zeros(tmax)  # Torque [Nm]
P_O = np.zeros(tmax)  # Open water propellor power
P_p = np.zeros(tmax)  # Propellor power [kW]
P_b = np.zeros(tmax)  # Engine brake power [kW]
P_T = np.zeros(tmax)  # Thrust power [kW]
P_E = np.zeros(tmax)  # Engine power [kW]
J = np.zeros(tmax)  # Advance ratio [-]

eta_e = np.zeros(tmax)
eta_q = np.zeros(tmax)
eta_m = np.zeros(tmax)
W_e = np.zeros(tmax)

# ------------- Run simulation ---------------
start = time.perf_counter()

for k in range(tmax - 1):
    # advance ratio
    J[k + 1] = ((v_a[k] / n_p[k]) / D_p)
    # Thrust and torque
    F_prop[k] = (((k_t_a * J[k + 1] ** 2) + (k_t_b * J[k + 1]) + k_t_c) * n_p[k] ** 2) * rho_sw * D_p ** 4
    M_prop[k] = ((((k_q_d * J[k + 1] ** 2) + (k_q_e * J[k + 1]) + k_q_f) * n_p[k] ** 2) * rho_sw * D_p ** 5) / eta_R
    KT[k + 1] = ((k_t_a * J[k + 1] ** 2) + (k_t_b * J[k + 1]) + k_t_c)
    KQ[k + 1] = ((k_q_d * J[k + 1] ** 2) + (k_q_e * J[k + 1]) + k_q_f)
    P_O[k + 1] = ((((k_q_d * J[k + 1] ** 2) + (k_q_e * J[k + 1]) + k_q_f) * n_p[k] ** 2) * rho_sw * D_p ** 5) * \
        n_p[k] * 2 * math.pi
    P_p[k + 1] = M_prop[k] * n_p[k] * 2 * math.pi
    # Calculate acceleration from resulting force --> ship speed & tr.distance
    sum_a[k + 1] = ((F_prop[k] - (R[k] / (1 - t))) / m_ship)
    # v_s_new = (np.trapz(sum_a[k:k+2], dx=0.01)) + v_s[k]
    v_s[k + 1] = integrate.simps(sum_a[:k + 2], dx=0.01)
    # v_s[k+1] = v_s_new
    Rsp[k + 1] = R[k] / (1 - t)
    # Traveled distance
    s[k + 1] = s[k] + v_s[k + 1] * dt
    # Advance velocity
    v_a[k + 1] = v_s[k + 1] * (1 - w)
    P_T[k + 1] = F_prop[k] * v_a[k + 1]
    # Resistance
    cws = np.interp(v_s[k + 1], snelheden, wrijvingscoef)
    Y = ov_Y_set[k]
    R[k + 1] = r_schip(v_s[k + 1], cws)
    P_E[k + 1] = v_s[k + 1] * R[k + 1]
    R_vorm[k + 1] = vormweerstand(v_s[k + 1])
    R_wrijving[k + 1] = wrijfweerstand(v_s[k + 1])
    R_golf[k + 1] = golfweerstand(v_s[k + 1], cws)
    # Calculate acceleration from resulting force --> propellor np
    sum_dnpdt[k + 1] = ((M_b[k] * i_gb * eta_TRM) - M_prop[k]) / (2 * math.pi * I_tot)
    n_p[k + 1] = integrate.simps(sum_dnpdt[:k + 2], dx=0.01) + n_p[0]
    # Engine speed
    n_e[k + 1] = n_p[k + 1] * i_gb
    # Fuel rack
    X = ov_X_set[k]
    m_flux_f[k + 1] = (X * m_f_nom * n_e[k + 1]) * i / k_es
    # Fuel consumption
    out_fc[k + 1] = integrate.simps(m_flux_f[:k + 2], dx=0.01) + out_fc[0]
    Q_f = X * m_f_nom * LHV
    eta_e[k + 1] = effective_engine_efficiency(n_e[k + 1], X)
    W_e[k + 1] = Q_f * eta_e[k + 1]
    # Brake power
    P_b[k + 1] = (W_e[k + 1] * n_e[k + 1] * i) / k_es
    # Engine torque
    M_b[k + 1] = P_b[k + 1] / (2 * math.pi * n_e[k + 1])
    M_Trm[k + 1] = M_b[k + 1] * i_gb * eta_TRM
    eta_q[k + 1] = heat_eff(X)
    eta_m[k + 1] = mech_eff(n_e[k + 1], X)

# EU just to be sure
v_s[0] = 0
v_s[1] = 0
# -------------- Plot Figure -------------------------------------------------
fig = plt.figure()

# Figuur 1
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(P_b[15:], eta_e[15:], label='Effectief')
ax1.plot(P_b[15:], eta_td[15:], label='Thermodynamisch')
ax1.plot(P_b[15:], eta_q[15:], label='Warmte')
ax1.plot(P_b[15:], eta_m[15:], label='Mechanisch')
ax1.set(xlabel='Vermogen',
        ylabel='Rendement [-]')
ax1.grid()
ax1.legend()

fig.savefig('fig/TR4.png')
