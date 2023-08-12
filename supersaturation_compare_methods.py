# Suppress warnings
import warnings 
warnings.simplefilter('ignore')

import math
import matplotlib.pyplot as plt
import numpy as np
import pyrcel as pm
import time

from pyrcel import binned_activation

# Sea level: 101325 Pa
P0 = 95000. # Pressure, Pa
T0 = 285.   # Temperature, K
S0 = -0.05  # Supersaturation, 1-RH (95% here)

# Mean and standard deviation constants for sea salt
# and large and small sulfate particles
mean_sea, stdev_sea = (0.307, 2.)
mean_lg_sulf, stdev_lg_sulf = (0.075, 1.6)
mean_sm_sulf, stdev_sm_sulf = (0.03, 1.6)

# Concentration values for sea salt and large/small sulfate
Nt_sea = 10.
Nt_lg_sulf= 50.
Nt_sm_sulf = 50.

# Hygroscopicity constant (κ) values for sea salt and sulfate
kappa_sea = 1.2
kappa_sulf = 0.54


# ---- CONSTANTS ----

# vertical wind speed
w = 0.5 # meters per second

# gravitational constant
g = 9.8 # m / s^2

# specific heat of dry air at constant pressure
c_p = 1004.0 # J / kg

# latent heat of water vaporization
L_v = 2.5e6

# diffusivity of water vapor in air
D_v = 3e-5 # m^2 / s

# density of water at STP
rho_liquid = 1000. # kg / m^3

# gas constant for dry air
R_d = 287.0 # J/kg/K

# gas constant for water vapor
R_v = 461.50 # J/kg/K 

# Constants
beta_val = 0.5
Mw = 18.0 / 1e3
R = 8.314


# ---- FUNCTION CALCULATIONS ----

def get_e_sat(T):
    # saturation vapor pressure over liquid in Pascal
    # Following IFS documentation Cy47r3, p. 211
    # https://www.ecmwf.int/en/elibrary/81271-ifs-documentation-cy47r3-part-iv-physical-processes
    T0 = 273.16
    return 611.21 * np.exp(17.502 * (T - T0) / (T - 32.19))


def get_q_sat(T, P): # q_sat does not depend on RH
    e_sat = get_e_sat(T)  
    return 0.622 * (e_sat / (P - e_sat))


def get_rho_air(T, P, q_v):
    # Tv is the "virtual temperature"
    Tv = T * (1.0 + 0.61 * q_v)
    rho_a = P / R_d / Tv  # air density

    return rho_a


def get_desat_dT(T):
    esat = get_e_sat(T)
    return L_v * esat / ( R_v * T * T )


def get_dq_sat_dT(T, P):
    # partial derivative of saturation vapor pressure over liquid
    #    with respect to temperature, units of kg/kg/K
    # Following IFS documentation Cy47r3, p. 212, eqn 12.20
    # https://www.ecmwf.int/en/elibrary/81271-ifs-documentation-cy47r3-part-iv-physical-processes
    qsat = get_q_sat(T, P)
    return L_v * qsat / ( R_v * T * T )


# ---- UNKNOWNS ----

N_c = 1e8 # TODO: numerically determine the appropriate value to use for this


# ---- PRIMARY DIFFERENTIAL EQUATION ----

def evolution(y, t):
    #   delta is defined as q_v - q_sat(p,T) ~ S*q_sat(p,T)
    z, p, T, q_v, q_c, q_i, delta = y
    e_sat = get_e_sat(T)
    rho_air = get_rho_air(T, p, q_v)
    dq_sat_dT = get_dq_sat_dT(T, p)

    if (q_c > 0):
        N_c = 1e8
        r_v = (3 * rho_air * q_c / (4 * np.pi * N_c * rho_liquid)) ** (1/3)
        tau_c = 1 / (4 * np.pi * D_v * N_c * r_v)
    elif (delta > 0):
        # estimate r_v using q_c = delta
        N_c = 1e8
        r_v = (3 * rho_air * delta / (4 * np.pi * N_c * rho_liquid)) ** (1/3)
        tau_c = 1 / (4 * np.pi * D_v * N_c * r_v)
    else:   
        tau_c = 1e8

    q_s = get_q_sat(T,p)
    A_delta = (g / c_p) * dq_sat_dT - q_s * rho_air * g / (p - e_sat)
    Gamma = 1 + (L_v / c_p) * dq_sat_dT
    C = delta / (tau_c * Gamma)

    dydt = [
        w,
        -rho_air * g * w,
        -g * 2/c_p + (L_v/c_p) * C,
        -C,
        C,
        0,
        A_delta * w - delta/tau_c,
        N0
    ]

    return dydt


prev_N_c = None

def evolution_increment(y, delta_t):
    # This function returns an increment of y over a step delta_t

    # delta is defined as q_v - q_sat(p,T) ~ S*q_sat(p,T)
    z, p, T, q_v, q_c, q_i, delta = y
    e_sat = get_e_sat(T)
    rho_air = get_rho_air(T, p, q_v)
    dq_sat_dT = get_dq_sat_dT(T, p)

    if (q_c > 0):
        N_c = 1e8
        r_v = (3 * rho_air * q_c / (4 * np.pi * N_c * rho_liquid)) ** (1/3)
        tau_c = 1 / (4 * np.pi * D_v * N_c * r_v)
    elif (delta > 0):
        N_c = 1e8
        r_v = (3 * rho_air * delta / (4 * np.pi * N_c * rho_liquid)) ** (1/3)
        tau_c = 1 / (4 * np.pi * D_v * N_c * r_v)
    else:   
        tau_c = 1e8

    q_s = get_q_sat(T,p)

    A_delta = (g / c_p) * dq_sat_dT - q_s * rho_air * g / (p - e_sat)
    Gamma = 1 + (L_v / c_p) * dq_sat_dT
    C = delta / (tau_c * Gamma)

    delta_final = A_delta*w*tau_c + (delta - A_delta*w*tau_c) * np.exp(- delta_t/tau_c)
    C_integ = 1/2 * (delta + delta_final) / Gamma / tau_c / 2

    C_integ = np.max([-q_c, np.min([q_v, C_integ])])

    delta_y = [
        w*delta_t,
        -rho_air * g * w * delta_t,
        -(g/c_p) * w * delta_t + (L_v/c_p) * C_integ,
        -C_integ,
        C_integ,
        0,
        delta_final - delta
    ]

    return delta_y


# ---- SOLUTION CALCULATION - PARCEL MODEL ----

lg_sulfate =  pm.AerosolSpecies('large sulfate', 
                             pm.Lognorm(mu=mean_lg_sulf, sigma=stdev_lg_sulf, N=Nt_lg_sulf),
                             kappa=kappa_sulf, bins=200)
sm_sulfate =  pm.AerosolSpecies('small sulfate', 
                             pm.Lognorm(mu=mean_sm_sulf, sigma=stdev_sm_sulf, N=Nt_sm_sulf),
                             kappa=kappa_sulf, bins=200)
sea_salt = pm.AerosolSpecies('sea salt',
                             pm.Lognorm(mu=mean_sea, sigma=stdev_sea, N=Nt_sea),
                             kappa=kappa_sea, bins=200)

# Initialize the parcel model and find the 
start_time = time.time()
initial_aerosols = [lg_sulfate, sm_sulfate, sea_salt]
V = 0.5 # updraft speed, m/s
dt = 1.0 # timestep, seconds
height = 300.
t_end = height/V # end time, seconds

model = pm.ParcelModel(initial_aerosols, V, T0, S0, P0, console=False, accom=0.3)
parcel_trace, aerosol_traces = model.run(t_end, dt, solver='cvode')

print('Parcel model execution time:', time.time() - start_time)

Smax = parcel_trace['S'].max()*100
z_at_smax = parcel_trace['z'].iloc[parcel_trace['S'].argmax()]


# ---- SATURATION PLOT ----

lg_sulf_array = aerosol_traces['large sulfate'].values
sm_sulf_array = aerosol_traces['small sulfate'].values
sea_array = aerosol_traces['sea salt'].values

fig, [axS, axA] = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

axS.plot(parcel_trace['S']*100., parcel_trace['z'], color='k', lw=2)
axT = axS.twiny()
axT.plot(parcel_trace['T'], parcel_trace['z'], color='r', lw=1.5)

axS.annotate("max S = %0.2f%%" % Smax, 
             xy=(Smax, z_at_smax), 
             xytext=(Smax-0.3, z_at_smax+50.),
             arrowprops=dict(arrowstyle="->", color='k',
                             connectionstyle='angle3,angleA=0,angleB=90'),
             zorder=10)

axS.set_xlim(0, 0.7)
axS.set_ylim(0, height)

axT.xaxis.label.set_color('red')
axT.tick_params(axis='x', colors='red')

axS.set_xlabel("Supersaturation, %")
axT.set_xlabel("Temperature, K")
axS.set_ylabel("Height, m")

# Color definitions for graphing
lg_sul_c = "#CC0066"
sm_sul_c = "#DF88AA"
sea_c = "#0099FF"

lg_ss = axA.plot(lg_sulf_array[:, ::10]*1e6, parcel_trace['z'], color=lg_sul_c, 
         label="large sulfate")
sm_ss = axA.plot(sm_sulf_array[:, ::10]*1e6, parcel_trace['z'], color=sm_sul_c, 
         label="small sulfate")
sa = axA.plot(sea_array[:, ::10]*1e6, parcel_trace['z'], color=sea_c, label="sea salt")
axA.semilogx()
axA.set_xlim(1e-2, 10.)
axA.set_xticks([1e-2, 1e-1, 1e0, 1e1], [0.01, 0.1, 1.0, 10.0])
axA.legend([lg_ss[0], sm_ss[0], sa[0]], ['large sulfate', 'small sulfate', 'sea salt'], loc='upper right')
axA.set_xlabel("Droplet radius, micron")

for ax in [axS, axA, axT]:
    ax.grid()

plt.show()


# ---- SOLUTION CALCULATION - EULER'S METHOD ----

Dv = 3.0e-5
N0 = 1e8
H = 300

def num_CCN_activated(s, T, mean, stdev, Nt, kappa):
    # https://journals.ametsoc.org/view/journals/atsc/65/3/2007jas2374.1.xml
    # Equation (6)
    sigma_w =  0.0761 - 1.55e-4 * (T - 273.15)
    A =  (2.0 * Mw * sigma_w) / (R * T * rho_liquid * 1e-6)
    # print('A:', A)

    s0 = mean ** -(1 + beta_val) * np.sqrt(4 * (A ** 3) / (27 * kappa))
    sig = stdev ** (1 + beta_val)

    if s <= 0:
        return 0
    u = np.log(s0 / s) / (np.sqrt(2) * np.log(sig))
    return (Nt / 2) * math.erfc(u)


def droplet_size_check(population, aerosol, threshold):
    return np.sum(aerosol.Nis[population > threshold])


def get_tau(N, mean):
    return 1 / (4 * np.pi * Dv * N * mean)


def get_A(S_0, S_max, tau, t_max):
    return (S_max - S_0 * np.exp(-t_max / tau)) / (tau * (1 - np.exp(-t_max / tau)))


# Isolate altitude, supersaturation, and temperature values
z_vals = parcel_trace['z']
S_vals = parcel_trace['S']
T_vals = parcel_trace['T']
S_maximizer = np.argmax(S_vals)

# Isolate aerosol data for different aerosol types
lg_sulf_trace = aerosol_traces['large sulfate']
sm_sulf_trace = aerosol_traces['small sulfate']
sea_trace = aerosol_traces['sea salt']

ind_final = int(t_end/dt) - 1

lg_sulf_vals = []
sm_sulf_vals = []
sea_vals = []

index = final_time_index = len(T_vals) - 1
lg_sulf_wet_radii = sm_sulf_wet_radii = sea_wet_radii = None


# Solve the differential equation for three different aerosol types
# Uses a basic Euler's method approach and prints the time elapsed
start_time = time.time()

lg_sulf_wet_radii = lg_sulf_trace.iloc[index]
sm_sulf_wet_radii = sm_sulf_trace.iloc[index]
sea_wet_radii = sea_trace.iloc[index]
T = parcel_trace['T'].iloc[index]

eq_lg_sulf, kn_lg_sulf, alpha_lg_sulf, phi_lg_sulf = \
    binned_activation(Smax/100, T, lg_sulf_trace.iloc[index], lg_sulfate)
eq_lg_sulf *= lg_sulfate.total_N
kn_lg_sulf *= lg_sulfate.total_N

eq_sm_sulf, kn_sm_sulf, alpha_sm_sulf, phi_sm_sulf = \
    binned_activation(Smax/100, T, sm_sulf_trace.iloc[index], sm_sulfate)
eq_sm_sulf *= sm_sulfate.total_N
kn_sm_sulf *= sm_sulfate.total_N

eq_sea, kn_sea, alpha_sea, phi_sea = \
    binned_activation(Smax/100, T, sea_trace.iloc[index], sea_salt)
eq_sea *= sea_salt.total_N
kn_sea *= sea_salt.total_N

lg_sulf_vals.append(kn_lg_sulf)
sm_sulf_vals.append(kn_sm_sulf)
sea_vals.append(kn_sea)

y0 = np.array([parcel_trace[axis][0] for axis in ['z', 'P', 'T', 'wv', 'wc', 'wi']]
              + [parcel_trace['S'][0] * get_q_sat(parcel_trace['T'][0], parcel_trace['P'][0])])

Tfinal = H / w
Nt = int(2*Tfinal+1)
t = np.linspace(0, Tfinal, Nt)
delta_t = t[2] - t[1]

sol = np.zeros((Nt,7))
sol[0,:] = y0
for it in range(1,Nt):
    y = sol[it-1,:]
    sol[it,:] = y + evolution_increment(y,delta_t)

sol_sea, sol_lg_sulf, sol_sm_sulf = [np.zeros((Nt,7))] * 3
sol_sea[0,:], sol_lg_sulf[0,:], sol_sm_sulf[0,:] = [y0] * 3
sol_sea[0,-1], sol_lg_sulf[0,-1], sol_sm_sulf[0,-1] = (Nt_sea, Nt_lg_sulf, Nt_sm_sulf)
for it in range(1, Nt):
    y_sea = sol_sea[it-1,:]
    y_lg_sulf = sol_lg_sulf[it-1,:]
    y_sm_sulf = sol_sm_sulf[it-1,:]

    sol_sea[it,:] = y_sea + evolution_increment(y_sea, delta_t)
    sol_lg_sulf[it,:] = y_lg_sulf + evolution_increment(y_lg_sulf, delta_t)
    sol_sm_sulf[it,:] = y_sm_sulf + evolution_increment(y_sm_sulf, delta_t)

# Show time elapsed for calculating using this model
print('Euler\'s method execution time:', time.time() - start_time)


# ---- COMPARATIVE PLOTS ----

plt.plot(z_vals, S_vals, label='Pyrcel supersaturation')

zmg = sol[:,0]
temps = sol[:,2]
supersat = sol[:,6] / get_q_sat(sol[:,2], sol[:,1])

plt.title('Supersaturation over time')
plt.xlabel('Time')
plt.ylabel('Supersaturation')
plt.plot(zmg, supersat, label='Mor&Grab supersaturation')
plt.legend(loc='best')
plt.show()

plt.title('q_v value over time - two models')
plt.xlabel('q_v')
plt.ylabel('Supersaturation')
plt.plot(z_vals, parcel_trace['wv'], label='Pyrcel q_v')
plt.plot(zmg, sol[:, 3], label='Mor&Grab q_v')
plt.legend(loc='best')
plt.show()

plt.title('Temperature over time - two models')
plt.xlabel('Temperature')
plt.ylabel('Supersaturation')
plt.plot(z_vals, parcel_trace['T'], label='Pyrcel Temperature')
plt.plot(zmg, sol[:, 2], label='Mor&Grab Temperature')
plt.legend(loc='best')
plt.show()

thresholds = np.logspace(-9, -4, 500)
total_Nis_sea = sum(sea_salt.Nis)
total_Nis_lg_sulf = sum(lg_sulfate.Nis)
total_Nis_sm_sulf = sum(sm_sulfate.Nis)

plt.title('Sea salt activation model comparison')
sea_supersat = np.maximum.accumulate([num_CCN_activated(supersat[i], temps[i], mean_sea, stdev_sea, Nt_sea, kappa_sea) for i in range(1200)])
plt.plot(np.array(range(1200)) / 2, sea_supersat, label='Mor&Grab Eq6 salt')
plt.plot(range(600), np.maximum.accumulate([num_CCN_activated(S_vals[i], T_vals[i], mean_sea, stdev_sea, Nt_sea, kappa_sea) for i in range(600)]), label='Mor&Grab with Pyrcel supersat')
plt.plot(range(600), Nt_sea * np.array([droplet_size_check(sea_trace.iloc[i], sea_salt, 1e-6) for i in range(600)]) / total_Nis_sea, label='Pyrcel model salt')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Activated N')
plt.show()

plt.title('Large sulfate activation model comparison')
lg_sulf_supersat = np.maximum.accumulate([num_CCN_activated(supersat[i], temps[i], mean_lg_sulf, stdev_lg_sulf, Nt_lg_sulf, kappa_sulf) for i in range(1200)])
plt.plot(np.array(range(1200)) / 2, lg_sulf_supersat, label='Mor&Grab Eq6 large sulfate')
plt.plot(range(600), np.maximum.accumulate([num_CCN_activated(S_vals[i], T_vals[i], mean_lg_sulf, stdev_lg_sulf, Nt_lg_sulf, kappa_sulf) for i in range(600)]), label='Mor&Grab with Pyrcel supersat')
plt.plot(range(600), Nt_lg_sulf * np.array([droplet_size_check(lg_sulf_trace.iloc[i], lg_sulfate, 1e-6) for i in range(600)]) / total_Nis_lg_sulf, label='Pyrcel model large sulfate')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Activated N')
plt.show()

plt.title('Small sulfate activation model comparison')
sm_sulf_supersat = np.maximum.accumulate([num_CCN_activated(supersat[i], temps[i], mean_sm_sulf, stdev_sm_sulf, Nt_sm_sulf, kappa_sulf) for i in range(1200)])
plt.plot(np.array(range(1200)) / 2, sm_sulf_supersat, label='Mor&Grab Eq6 small sulfate')
plt.plot(range(600), np.maximum.accumulate([num_CCN_activated(S_vals[i], T_vals[i], mean_sm_sulf, stdev_sm_sulf, Nt_sm_sulf, kappa_sulf) for i in range(600)]), label='Mor&Grab with Pyrcel supersat')
plt.plot(range(600), Nt_sm_sulf * np.array([droplet_size_check(sm_sulf_trace.iloc[i], sm_sulfate, 1e-6) for i in range(600)]) / total_Nis_sm_sulf, label='Pyrcel model small sulfate')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Activated N')
plt.show()


# ---- FINAL CONSOLE OUTPUT ----

print('\n------------------------\nParcel model concentration results:')
print("  CDNC(large sulfate) = {:3.1f}".format(eq_lg_sulf))
print("  CDNC(small sulfate) = {:3.1f}".format(eq_sm_sulf))
print("  CDNC(sea salt) = {:3.1f}".format(eq_sea))
print("------------------------")

# Display the total number of particles activated, based on the 
print("          total = {:3.1f} / {:3.0f} ~ act frac = {:1.2f}".format(
      kn_lg_sulf + kn_sm_sulf + kn_sea, 
      sea_salt.total_N + lg_sulfate.total_N + sm_sulfate.total_N,
      (kn_lg_sulf + kn_sm_sulf + kn_sea)/(sea_salt.total_N + lg_sulfate.total_N + sm_sulfate.total_N)
))