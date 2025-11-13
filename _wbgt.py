"""
Created by Riddhima Puri on 2/12/2024.
Adapted using PyWBGT from Kong et al. 2022. 

"""
import math
import numpy as np
import xarray as xr
import pandas as pd
import netCDF4

# import xclim
import cftime

from xclim.indices.helpers import (
    _gather_lat,
    _gather_lon,
    cosine_of_solar_zenith_angle,
    distance_from_sun,
    solar_declination,
    wind_speed_height_conversion,
)

from scipy import optimize


# physical constants
mair = 28.97 # molecular weight of dry air (grams per mole)
mh2o = 18.015 # molecular weight of water vapor (grams per mole)
rgas = 8314.34 # ideal gas constant (J/kg mol · K)
cp = 1003.5 # Specific heat capacity of air at constant pressure (J·kg-1·K-1)
stefanb = 0.000000056696  # stefan-boltzmann constant
ratio = cp * mair * (mh2o**(-1))
rair = rgas * (mair**(-1))
Pr = cp * ((cp + 1.25 * rair)**(-1)) # Prandtl number 

#globe constants
diamglobe = 0.0508 # diameter of globe (m)
emisglobe = 0.95 # emissivity of globe
albglobe = 0.05 # albedo of globe

#wick constants
emiswick = 0.95 # emissivity of the wick
albwick = 0.4 # albedo of the wick
diamwick = 0.007 # diameter of the wick
lenwick = 0.0254 # length of the wick

#surface constant
albsfc=0.45

PI=3.1415926535897932


# Arguments for equations of Tg and Tnwb.
class Tg_params:
    def __init__(self):
        self.C0 = 0.0
        self.C1 = 0.0
        self.C2 = 0.0
        self.C3 = 0.0

class Tnwb_params:
    def __init__(self):
        self.D0 = 0.0
        self.D1 = 0.0
        self.D2 = 0.0
        self.D3 = 0.0
        self.D4 = 0.0


# Calculation of physical quantities.
def viscosity(tas):
    """    
    Parameters:
    tas: air temperature (K)

    Returns: air viscosity (kg/(m s))
    """
    omega = 1.2945 - tas / 1141.176470588
    visc = 0.0000026693 * (math.sqrt(28.97 * tas)) * ((13.082689 * omega) ** -1)
    return visc


def thermcond(tas):
    """
    Parameters:
    tas: air temperature (K)

    Returns: thermal conductivity of air (W/(m K))
    """
    tc = (cp + 1.25 * rair) * viscosity(tas)
    return tc

def esat(tas, ps):
    """
    Parameters:
    tas: air temperature (K)
    ps: surface pressure (Pa)
    
    Returns: saturation vapor pressure (Pa)
    """
    if tas > 273.15:
        es = 611.21 * math.exp(17.502 * (tas - 273.15) / (tas - 32.18))
        es = (1.0007 + (3.46e-6 * ps / 100)) * es
    else:
        es = 611.15 * math.exp(22.452 * (tas - 273.15) / (tas - 0.6))
        es = (1.0003 + (4.18e-6 * ps / 100)) * es
    return es

def emisatm(tas, hurs, ps):
    """
    Parameters:
    tas: air temperature (K)
    hurs: relative humidity (%)
    ps: surface pressure (Pa)
    
    Returns: atmospheric emissivity
    """
    e = hurs * 0.01 * (esat(tas, ps) * 0.01)
    emis_atm = 0.575 * (e ** 0.143)
    return emis_atm

def diffusivity(tas, ps):
    """
    Parameters:
    tas: air temperature (K)
    ps: surface pressure (Pa)
    
    Returns: diffusivity of water vapor in air (m^2/s)
    """
    return 2.471773765165648e-05 * ((tas * 0.0034210563748421257) ** 2.334) * ((ps / 101325) ** -1)

def h_evap(tas):
    """
    Parameters:
    tas: air temperature (K)

    Returns: heat of evaporation (J/(kg K))
    """
    return 1665134.5+2370.0*tas

def h_sphere_in_air(tas, ps, sfcwind):
    """
    Parameters:
    tas: air temperature (K)
    ps: surface pressure (Pa)
    sfcwind: 2 meter wind (m/s)
    
    Returns: convective heat tranfer coefficient for flow around a sphere (W/(m2 K))
    """
    thermcon = thermcond(tas)
    density = ps * ((rair * tas)**(-1))
    Re = sfcwind * density * diamglobe * ((viscosity(tas))**(-1))
    Nu = 2 + 0.6 * math.sqrt(Re) * math.pow(Pr,0.3333)
    h = Nu * thermcon * (diamglobe**(-1))
    return h

def h_cylinder_in_air(tas, ps, sfcwind):
    """
    Parameters:
    tas: air temperature (K)
    ps: surface pressure (Pa)
    sfcwind: 2 meter wind (m/s)

    Returns: convective heat transfer coefficient for a long cylinder (W/(m2 K))
    """
    thermcon = thermcond(tas)
    density = ps * ((rair * tas)**(-1))
    Re = sfcwind * density * diamwick * ((viscosity(tas))**(-1))
    Nu = 0.281 * (Re ** 0.6) * (Pr ** 0.44)
    h = Nu * thermcon * (diamwick**(-1))
    return h


# Calculate Tg.
def fTg(x, args):
    """
    
    Equation of Tg that needs to be solved by iteration.

    Parameters:
    x: variable to solve for
    args: arguments of the equation, from Tg_params

    Returns: the equation
    """

    myargs = args
    h = h_sphere_in_air(0.5 * (myargs.C1 + x), myargs.C2, myargs.C3)
    return (myargs.C0 - ((emisglobe * stefanb) ** -1) * h * (x - myargs.C1)) - x ** 4

def Tg_brentq_wrapper(args, xa, xb, xtol, rtol, mitr):
    """
    Use scipy.optimize.brentq algorithm to solve Tg iteratively.
    
    Parameters:
    args (Tg_Params): parameters for the equation
    xa: lower bound of the interval
    xb: upper bound of the interval
    xtol: absolute error tolerance
    rtol: relative error tolerance
    mitr: maximum number of iterations
    
    Returns: root of the equation
    """
    if rtol < 1e-5:
        rtol = 1e-5  # Set a minimum rtol value to avoid the error
    return optimize.brentq(fTg, xa, xb, args=(args,), xtol=xtol, rtol=rtol, maxiter=mitr)

def Tg_10mwind(tas, ps, sfcwind, rsds, rsus, rlds, rlus, fdir, cosz, xtol=0.01, rtol=0.0, mitr=100000):
    """
    Build and solve the Tg equation for 10m wind speed.
    
    Parameters:
    tas: near-surface air temperature (K)
    sfcwind: near-surface meter wind speed (m/s) (2m in Lilgrejen's code?)
    ps: near-surface pressure (Pa)
    rsds: surface downward solar radiation (W/m2)
    rsus: surface reflected solar radiation (W/m2)
    rlds: surface downward long-wave radiation (W/m2)
    rlus: surface upwelling long-wave radiation (W/m2)
    fdir: the ratio of direct solar radiation 
    cosz: average cosine zenith angle during only the sunlit period of each interval

    Returns: outdoor black globe temperature (K)

    Modified from WBGT.pyx under https://github.com/QINQINKONG/PyWBGT/tree/v1.0.0 """

    x_max = tas.shape[0]
    y_max = tas.shape[1]
    z_max = tas.shape[2]
    args = Tg_params()
    result = np.zeros([x_max, y_max, z_max])

    for i in range(x_max):
        for j in range(y_max):
            for k in range(z_max):
                    if (np.isnan(tas[i,j,k]) or np.isnan(ps[i,j,k]) or np.isnan(sfcwind[i,j,k]) or np.isnan(rsds[i,j,k]) or np.isnan(rsus[i,j,k]) or np.isnan(rlds[i,j,k]) or np.isnan(rlus[i,j,k]) or np.isnan(fdir[i,j,k]) or np.isnan(cosz[i,j,k])):
                        result[i,j,k] = math.nan
                    else:
                        args.C0 = 0.5*(stefanb**(-1))*(rlds[i,j,k]+rlus[i,j,k])+rsds[i,j,k]*((2*emisglobe*stefanb)**(-1))*(1-albglobe)*(1-fdir[i,j,k]+0.5*fdir[i,j,k]*(cosz[i,j,k]**(-1)))+(1-albglobe)*((2*emisglobe*stefanb)**(-1))*rsus[i,j,k]
                        args.C1 = tas[i,j,k]
                        args.C2 = ps[i,j,k]
                        args.C3 = wind_speed_height_conversion(sfcwind[i,j,k], h_source="10 m", h_target="2 m", method = 'log')
                        xa=tas[i,j,k]-50
                        xb=tas[i,j,k]+90
                        # result[i,j,k]=Tg_brentq_wrapper(args, xa, xb, xtol, rtol, mitr)
                        try:
                            result[i, j, k] = Tg_brentq_wrapper(args, xa, xb, xtol, rtol, mitr)
                        except ValueError as e:
                            print(f"Fallback for i={i}, j={j}, k={k}: {e}")
                            result[i, j, k] = np.nan  # Assign a fallback value
    return result

# Calculate Tnwb.
def fTnwb(x, args):
    """
    Equation of Tnwb that needs to be solved by iteration.
    
    Parameters:
    x: variable to solve for
    args: arguments of the equation, from Tnwb_params

    Returns: the equation
    """

    myargs = args
    evap = h_evap(0.5 * (x + myargs.D0))
    es = esat(x, myargs.D1)
    density = myargs.D1 / (0.5 * (myargs.D0 + x) * rair)
    Sc = viscosity(0.5 * (myargs.D0 + x)) / (density * diffusivity(0.5 * (myargs.D0 + x), myargs.D1))
    h = h_cylinder_in_air(0.5 * (myargs.D0 + x), myargs.D1, myargs.D3)
    Fatm = myargs.D4 - emiswick * stefanb * (x ** 4)
    return myargs.D0 - evap / ratio * (es - myargs.D2) / (myargs.D1 - es) * (Pr / Sc ** 0.56) + Fatm / h - x

def Tnwb_brentq_wrapper(args, xa, xb, xtol, rtol, mitr):
    """
    Use scipy.optimize.brentq algorithm to solve Tnwb iteratively.
    
    Parameters:
    args (TnwbParams): Parameters for the equation
    xa: lower bound of the interval
    xb: upper bound of the interval
    xtol: absolute error tolerance
    rtol: relative error tolerance
    mitr: maximum number of iterations
    
    Returns: root of the equation
    """
    if rtol < 1e-5:
        rtol = 1e-5  # Set a minimum rtol value to avoid the error
    
    return optimize.brentq(fTnwb, xa, xb, args=(args,), xtol=xtol, rtol=rtol, maxiter=mitr)


def Tnwb_10mwind(tas, hurs, ps, sfcwind, rsds, rsus, rlds, rlus, fdir, cosz, xtol=0.01, rtol=0.0, mitr=100000):
    """
    Build and solve the Tg equation for 10m wind speed.
    
    Parameters:
    tas: near-surface air temperature (K)
    hurs: relative humidity (%)
    sfcwind: near-surface meter wind speed (m/s) (2m in Lilgrejen's code?)
    ps: near-surface pressure (Pa)
    rsds: surface downward solar radiation (W/m2)
    rsus: surface reflected solar radiation (W/m2)
    rlds: surface downward long-wave radiation (W/m2)
    rlus: surface upwelling long-wave radiation (W/m2)
    fdir: the ratio of direct solar radiation 
    cosz: average cosine zenith angle during only the sunlit period of each interval

    Returns: outdoor natural wet bulb temperature (K)
    """
    x_max = tas.shape[0]
    y_max = tas.shape[1]
    z_max = tas.shape[2]
    args = Tnwb_params()
    result = np.zeros([x_max, y_max, z_max])

    for i in range(x_max):
        for j in range(y_max):
            for k in range(z_max):
                    if (np.isnan(tas[i,j,k]) or np.isnan(hurs[i,j,k]) or np.isnan(ps[i,j,k]) or np.isnan(sfcwind[i,j,k]) or np.isnan(rsds[i,j,k]) or np.isnan(rsus[i,j,k]) or np.isnan(rlds[i,j,k]) or np.isnan(rlus[i,j,k]) or np.isnan(fdir[i,j,k]) or np.isnan(cosz[i,j,k])):
                        result[i,j,k] = math.nan
                    else:
                        args.D0 = tas[i,j,k]
                        args.D1 = ps[i,j,k]
                        args.D2 = hurs[i,j,k]*0.01*esat(tas[i,j,k],ps[i,j,k])
                        args.D3 = wind_speed_height_conversion(sfcwind[i,j,k], h_source="10 m", h_target="2 m", method = 'log')
                        args.D4 = emiswick*0.5*(rlds[i,j,k]+rlus[i,j,k])+(1+diamwick*((4*lenwick)**(-1)))*(1-albwick)*(1-fdir[i,j,k])*rsds[i,j,k]+(math.tan((math.acos(cosz[i,j,k])))*(PI**(-1))+diamwick*((4*lenwick)**(-1)))*(1-albwick)*fdir[i,j,k]*rsds[i,j,k]+(1-albwick)*rsus[i,j,k]
                        xa = tas[i,j,k]-((100-hurs[i,j,k])/5.0)-50
                        xb = min(tas[i,j,k]+70,340.0)
                        # result[i,j,k]=Tnwb_brentq_wrapper(args, xa, xb, xtol, rtol, mitr)
                        try:
                            result[i, j, k] = Tnwb_brentq_wrapper(args, xa, xb, xtol, rtol, mitr)
                        except ValueError as e:
                            print(f"Fallback for i={i}, j={j}, k={k}: {e}")
                            result[i, j, k] = np.nan  # Assign a fallback value
    return result


# Calculate WBGT. 
def WBGT(tas,hurs,ps,sfcwind,rsds,rsus,rlds,rlus,fdir,cosz):
    """ 
    
    Parameters:
    tas: air temperature (K)
    hurs: relative humidity (%)
    sfcwind: 10 meter wind speed (m/s) (to be converted to 2m)
    ps: surface pressure (Pa)
    rsds: surface downward solar radiation (W/m2)
    rsus: surface reflected solar radiation (W/m2)
    rlds: surface downward long-wave radiation (W/m2)
    rlus: surface upwelling long-wave radiation (W/m2)
    fdir: the ratio of direct solar radiation 
    cosz: cosine solar zenith angle

    Returns: outdoor wet bulb globe temperature (K)
    """
    tg=Tg_10mwind(tas,ps,sfcwind,rsds,rsus,rlds,rlus,fdir,cosz)
    tnwb=Tnwb_10mwind(tas,hurs,ps,sfcwind,rsds,rsus,rlds,rlus,fdir,cosz)
    wbgt=0.7*tnwb+0.2*tg+0.1*tas
    return wbgt