# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:02:26 2017

    Viskotter simulation file
    Version 1.0H
    J. Rodrigues Monteiro, based on matlab code from P. de Vos
    Delft University of Technology
    3ME / MTT / SDPO / ME

Edited for assignment (changes denoted by '*')
    Alex Remstedt

History:
    20171108I: initial python version    JRM
    20200108J: simps integration         EU
              no graphing endpoits...   EU
    20200228H: simps correctie           EU
    *20210215: line 27, and 74-77 changed for assignment AR
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time

from scipy import integrate

# ----------- parameters for simulation --------------------------------------

tmax = 1500  # simulation time [s]
dt = 1  # timestep [s] *

# fuel properties
LHV = 42700  # Lower Heating Value [kJ/kg]
print('fuel properties loaded')

# water properties
rho_sw = 1025  # density of seawater [kg/m3]
print('water properties loaded')

# ship data
m_ship = 358000  # ship mass [kg]
c1 = 1500  # resistance coefficient c1 in R = c1*vs^2
v_s0 = 6.5430  # ship design speed [m/s]
t = 0.1600  # thrust deduction factor[-]
w = 0.2000  # wake factor [-]
print('ship data loaded')

# propellor data
D_p = 3  # diameter of propellor [m]
K_T_a = -0.3821  # factor a in K_T = a*J + b [-]
K_T_b = 0.2885  # factor b in K_T = a*J + b [-]
K_Q_a = -0.03346  # factor a in K_Q = a*J + b [-]
K_Q_b = 0.0308  # factor b in K_Q = a*J + b [-]
eta_R = 1.0100  # relative rotative efficiency [-]
print('propellor data loaded')

# engine data
m_f_nom = 1.314762  # nominal fuel injection [g]
eta_e = 0.3800  # nominal engine efficiency [-]
i = 6  # number of cylinders [-]
k_es = 2  # k-factor for engines based on nr.of strokes per cycle
P_b = np.zeros(tmax)  # engine power [kW]
P_b[0] = 960  # Nominal engine power [kW]
M_b = np.zeros(tmax)  # engine torque [Nm]
M_b[0] = P_b[0] * 1000 / 2 / math.pi / (900 / 60)  # ([P_b*1000/2/math.pi/n_eng_nom])
print('engine data loaded')

# gearbox data
eta_TRM = 0.9500  # transmission efficiency [-]
i_gb = 4.2100  # gearbox ratio [-]
I_tot = 200  # total mass of inertia of propulsion system [kg*m^2]
print('gearbox data loaded')

# initial values
in_p = 3.2830  # initial rpm
iv_t_control = np.array([0, 0.2 * tmax, 0.2 * tmax, 0.4 * tmax,  0.4 * tmax, 0.6 * tmax, 0.6 * tmax, 0.8 * tmax,
                         0.8 * tmax,  tmax])  # *
X_parms = np.array([1.0, 1.0, 0.75, 0.75, 0.5, 0.5, 0.25, 0.25, 0.01, 0.01])  # % maximum fuelrack *
Y_parms = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])  # disturbance factor *

# simulation control parameters
xvals = np.linspace(0, tmax - 1, tmax)
ov_X_set = np.interp(xvals, iv_t_control, X_parms)
ov_Y_set = np.interp(xvals, iv_t_control, Y_parms)


# --------- Start van de funtie definities

def R_schip(snelheid_schip):
    global Y, c1
    weerstand = Y * c1 * snelheid_schip ** 2
    return weerstand


# -------- Make arrays -------------------------------------------------------

# Time
mytime = np.linspace(0, tmax - 1, tmax)
# Velocity of the ship [m/s]
v_s = np.zeros(tmax)
v_s[0] = v_s0
# Distance traveled [m]
s = np.zeros(tmax)
# Advance velocity [m/s]
v_a = np.zeros(tmax)
v_a[0] = (1 - w) * v_s0
# Rpm propellor [Hz]
n_p = np.zeros(tmax)
n_p[0] = in_p
# Rpm diesel engine [Hz]
n_e = np.zeros(tmax)
n_e[0] = 900 / 60  # Nominal engine speed in rotations per second [Hz]
# Resistance [N]
R = np.zeros(tmax)
Y = ov_Y_set[0]
R[0] = R_schip(v_s0)
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

Q_f = np.zeros(tmax)
eta_TRM = np.full(tmax, 0.95)
eta_e = np.zeros(tmax)

# ------------- Run simulation -----------------------------------------------
start = time.perf_counter()

for k in range(tmax - 1):
    # advance ratio
    J[k + 1] = ((v_a[k] / n_p[k]) / D_p)
    # Thrust and torque
    F_prop[k] = ((((J[k + 1] * K_T_a) + K_T_b) *
                  n_p[k] ** 2) * rho_sw * D_p ** 4)
    M_prop[k] = (((((J[k + 1] * K_Q_a) + K_Q_b) *
                   n_p[k] ** 2) * rho_sw * D_p ** 5) / eta_R)
    KT[k + 1] = J[k + 1] * K_T_a + K_T_b
    KQ[k + 1] = J[k + 1] * K_Q_a + K_Q_b
    P_O[k + 1] = ((((J[k + 1] * K_Q_a) + K_Q_b) *
                   n_p[k] ** 2) * rho_sw * D_p ** 5) * n_p[k] * 2 * math.pi
    P_p[k + 1] = M_prop[k] * n_p[k] * 2 * math.pi
    # Calculate acceleration from resulting force --> ship speed & tr.distance
    sum_a[k + 1] = ((F_prop[k] - (R[k] / (1 - t))) / m_ship)
    # v_s_new = (np.trapz(sum_a[k:k+2], dx=0.01)) + v_s[k]
    v_s[k + 1] = integrate.simps(sum_a[:k + 2], dx=0.01) + v_s0
    # v_s[k+1] = v_s_new
    Rsp[k + 1] = R[k] / (1 - t)
    # Traveled distance
    s[k + 1] = s[k] + v_s[k + 1] * dt
    # Advance velocity
    v_a[k + 1] = v_s[k + 1] * (1 - w)
    P_T[k + 1] = F_prop[k] * v_a[k + 1]
    # Resistance
    Y = ov_Y_set[k]
    R[k + 1] = R_schip(v_s[k + 1])
    P_E[k + 1] = v_s[k + 1] * R[k + 1]
    # Calculate acceleration from resulting force --> propellor np
    sum_dnpdt[k + 1] = ((M_b[k] * i_gb * eta_TRM[k]) - M_prop[k]) / (2 * math.pi * I_tot)
    n_p[k + 1] = integrate.simps(sum_dnpdt[:k + 2], dx=0.01) + n_p[0]
    # Engine speed
    n_e[k + 1] = n_p[k + 1] * i_gb
    # Fuel rack
    X = ov_X_set[k]
    m_flux_f[k + 1] = (X * m_f_nom * n_e[k + 1]) * i / k_es
    # Fuel consumption
    out_fc[k + 1] = integrate.simps(m_flux_f[:k + 2], dx=0.01) + out_fc[0]
    Q_f[k] = X * m_f_nom * LHV
    eta_e[k] = 0.3800
    W_e = Q_f[k] * eta_e[k]
    # Brake power
    P_b[k + 1] = (W_e * n_e[k + 1] * i) / k_es
    # Engine torque
    M_b[k + 1] = P_b[k + 1] / (2 * math.pi * n_e[k + 1])
    eta_TRM[k] = 0.9500
    M_Trm[k + 1] = M_b[k + 1] * i_gb * eta_TRM[k]

# EU just to be sure
v_s[0] = v_s0
v_s[1] = v_s0
# -------------- Plot Figure -------------------------------------------------
# create figure with four subplots
fig = plt.figure(figsize=(10, 7.5))
fig2 = plt.figure(figsize=(10, 7.5))

# Figuur 1
ax1 = fig.add_subplot(3, 1, 1)  # fig.add_subplot(#rows, #cols, #plot)
ax1.plot(mytime[1:tmax-2], v_s[1:tmax-2])
ax1.set(ylabel='Scheeps snelheid [m/s]',
        xlabel='Tijd [s]')
ax1.grid()


# ax2 = fig.add_subplot(20, 2, 2)
# ax2.plot(mytime[1:tmax-2], s[1:tmax-2])
# ax2.set(ylabel='Afgelegde weg [m]',
#         xlabel='Tijd [s]')
# ax2.grid()

ax3 = fig.add_subplot(3, 1, 2)
ax3.plot(mytime[1:tmax-2], out_fc[1:tmax-2])
ax3.set(ylabel='Verbruikte Brandstof [g]',
        xlabel='Tijd [s]')
ax3.grid()

ax4 = fig.add_subplot(3, 1, 3)
ax4.plot(mytime[1:tmax-2], ov_X_set[1:tmax-2])
ax4.set(ylabel='Fuel rack [%]',
        xlabel='Tijd [s]')
ax4.grid()

# Figuur 2
# ax6 = fig2.add_subplot(20, 2, 5)
# ax6.plot(mytime[1 : tmax], 60 * n_p[1 : tmax], label = "Propeller")
# ax6.plot(mytime[1 : tmax], 60 * n_e[1 : tmax], label = "Diesel Motor")
# ax6.set(ylabel='Toerental [Hz]',
#         xlabel='Tijd [s]')
# ax6.grid()
# ax6.legend()

# figuur 3
# ax5 = fig.add_subplot(20, 2, 6)
# ax5.plot(v_s[1:tmax - 2], R[1:tmax - 2])
# ax5.set(xlabel='Scheeps snelheid [m/s]',
#         ylabel='Scheeps weerstand [N]')
# ax5.grid()


#figuur 4
# ax7 = fig.add_subplot(20, 2, 7)
# ax7.plot(v_a[1: tmax - 1], F_prop[1 : tmax - 1])
# ax7.set(ylabel = 'Voortstuwingskracht [N]',
#         xlabel = 'Voortgangssnelheid [m/s]')
# ax7.grid()

#figuur 5
# ax9 = fig.add_subplot(20, 2, 8)
# ax9.plot(60 * n_p[1:tmax - 2], M_prop[1:tmax - 2])
# ax9.set(xlabel='Propeller Toerental [Hz]',
#         ylabel='Propeller Koppel [Nm]')
# ax9.grid()

#figuur 6
# ax8 = fig.add_subplot(20, 2, 9)
# ax8.plot(60 * n_e[1 : tmax-2], M_b[1 : tmax-2])
# ax8.set(ylabel = 'Motor Koppel [Nm]',
#         xlabel = 'Diesel Motor Toerental [Hz]')
# ax8.grid()

#figuur 7

# ax10 = fig.add_subplot(20, 2, 10)
# ax10.plot(mytime[1 : tmax-2], P_E[1:tmax-2])
# ax10.set(xlabel='Tijd [s]',
#          ylabel='Motor Power [W]')
# ax10.grid()

# ax11 = fig.add_subplot(20, 2, 11)
# ax11.plot(mytime[1 : tmax-2], P_p[1 : tmax-2])
# ax11.set(xlabel = "Tijd [s]",
#          ylabel = "Propeller Power [W]")
# ax11.grid()

# ax12 = fig.add_subplot(20, 2, 12)
# ax12.plot(mytime[1 : tmax-2], P_b[1 : tmax-2])
# ax12.set(xlabel = "Tijd [s]",
#          ylabel = "Brake Power [W]")
# ax12.grid()

# ax13 = fig.add_subplot(20, 2, 13)
# ax13.plot(mytime[1: tmax - 2], Q_f[1 : tmax - 2])
# ax13.set(xlabel = "Tijd [s]",
#          ylabel = "Potentiele warmte [-]")   
# ax13.grid()

ax14 = fig2.add_subplot(3, 1, 1)
ax14.plot(mytime[1 : tmax-2], eta_TRM[1 : tmax-2])
ax14.set(xlabel = "Tijd [s]",
         ylabel = "Transmissie Rendement [-]")
ax14.grid()

ax15 = fig2.add_subplot(3, 1, 2)
ax15.plot(mytime[1 : tmax-2], eta_e[1 : tmax-2])
ax15.set(xlabel = "Tijd [s]",
         ylabel = "Nominaal Motor Rendement [-]")
ax15.grid()

ax16 = fig2.add_subplot(3, 1, 3)
ax16.plot(eta_e[1 : tmax -2], P_b[1 : tmax-2])
ax16.set(xlabel = "Nominaal Motor Rendement [-]",
         ylabel = "Arbeid van de motor [W]")
ax16.grid()
fig.tight_layout()
fig2.tight_layout()

fig.savefig('images/opdracht_1_1.png', bbox_inches='tight')
fig2.savefig('images/opdracht_1_2.png', bbox_inches='tight')
print(time.perf_counter()-start)
