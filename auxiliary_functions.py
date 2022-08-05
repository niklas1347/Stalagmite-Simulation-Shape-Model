#!/usr/bin/env python
# coding: utf-8

# In[49]:


import numpy as np
import sys
import pickle


# In[50]:


def Konstanten(TC, pCO2):
    #
    #EQUILIBRIUM code Dreybrodt1988/S15:
    #
    #function KONSTANTEN(TC, pCO2) calculates the calcium concentration [mol/m³] comprised in the
    #CO2-H2O-CaCO3 system.
    #
    #input Varibales:
    #TC = Temperature [°C]
    #pCO2 = Cave CO2  [atm]
    #
    #output Variables:
    #ca = calcium concentration [mol/m³]
    
    #fixed Parameters:
    Ca = 1e-4 #any starting concentration of calcium
    gammaCO2 = 1

    #temperature in Kelvin:
    TK = TC + 273.16
    
    #temperature (only) dependent variables   
    #mass action constants
    K1 = 10**(-356.3094 - 0.06091964 * TK + 21834.37 / TK + 126.8339 * np.log10(TK) - 1684915 / (TK**2)) #KAUFMANN2003
    K2 = 10**(-107.8871 - 0.03252849 * TK + 5151.79 / TK + 38.92561 * np.log10(TK) - 563713.9 / (TK**2))
    KC = 10**(-171.9065 - 0.077993 * TK + 2839.319 / TK + 71.595 * np.log10(TK))
    KH = 10**(108.3865 + 0.01985076 * TK - 6919.53 / TK - 40.45154 * np.log10(TK) + 669365 / (TK**2))

    #constants for activity coefficients
    A = 0.48809 + 8.074e-4 * TC
    B = 0.3241 + 1.600e-4 * TC
    
    #calculation of the concentrations in equilibrium at a given pCO2
    #Ca = 1e-3; %any starting concentration of calcium, hab ich oben definiert
    
    #calculation of the activity coefficients acoording to the boundary
    #condition IS = 3 * Ca. Equilibrium values establish after a few runs
    #(less than 10)
    
    for i in range(1,11):
        #ionic strength
        IS = 3 * Ca
        
        #activity coefficients
        gammaCa = 10**(-A * 4 * np.sqrt(IS) / (1 + B * 5 * np.sqrt(IS)) + 0.165*IS) #Nach Wateq, Truesdell1974
        gammaHCO = 10**(-A * np.sqrt(IS) / (1 + B * 5.4 * np.sqrt(IS))) #Nach Wateq, Truesdell1974 
        
        #calcium
        Ca = (pCO2 * K1 * KC * KH / (4 * K2 * gammaCa * gammaHCO**2))**(1/3)
    
    Ca = Ca*10**3 #to convert [mol/l] in [mol/m³]
    return Ca


# In[51]:


def alpha_poly(temp, delta):
# rate constant alpha for calcite precipitation as a function of
# temperature and film thickness, Note that alpha does not depend 
# on the CO2 Concentration of the Cave, Formula after "Testing 
# theoretically predicted stalgmite growth...; Baker; 1997" 

# Input Variables:
# temp: Temperature of the Cave [°C]
# delta: film thickness of the drop [m]

# Output Variables:
# alpha: rate constant [m/s]

    if delta == 1e-4:
        alpha = 1e-7*(0.51549+0.04015*temp+0.00418*temp*temp)
        
    elif delta == 7.5e-5:
        alpha = 1e-7*(0.4615+0.03192*temp+0.00408*temp*temp)
        
    elif delta == 5e-5:
        alpha = 1e-7*(0.43182+0.02103*temp+0.00381*temp*temp)
        
    else:
        sys.exit('filmthickness too small')
    
    return alpha


# In[52]:


def get_pickled_data(filename):
    #These Codes provide a routine to load data from CaveCalc and 
    #to get atmospheric CO2 data from an icecore
    #
    #Input Variables:
    #filename = pathname to the CaveCalc file
    #
    #Output Variables:
    #Age = Alter [y]
    #Ca = calcium concentration of the water drop [mol/m³]
    #CO2 = CO2 concentration of the water drop [mol/m³]
    
    infile = open(filename,'rb')
    age_and_models = pickle.load(infile)
    infile.close()

    Age = age_and_models[0]
    models = age_and_models[1]
    Ca = []
    CO2 = []


    for model in models:

        data = model[0]
        try:
            Ca.append(data['Ca(mol/kgw)'][1])
        except IndexError:
            Ca.append(np.nan)
        try:
            CO2.append(data['CO2(mol/kgw)'][1])
        except IndexError:
            CO2.append(np.nan)


    Age = np.array(Age)
    Ca = np.array(Ca)*999.975 #[mol/m³]; 1 m^3=999.975 kg
    CO2 = np.array(CO2)*999.975 #[mol/m³]; 1 m^3=999.975 kg

    return Age, Ca, CO2 


# In[53]:


def same_timescale(age,Ca,CO2,age_CO2,scale,rescale,T=15):
    #this Code takes a timeseries for the age and the Ca and calculates mean values in timesteps of the scale
    #
    #Input:
    #age = age of the Data of the Sofular Spelothem which gives the Ca Data, array should start with old age and goes to young [y]
    #Ca = Calcium of the drip water calculated by CaveCalc [mol/m^3]
    #CO2 = atmospehric CO2 Data from an Icechore [µmol/mol = ppm]
    #age_CO2 = age of the Data from the CO2 Data, array should start with young age and end with old [y]
    #scale = timesteps of the age data, here we have 10 y and 5 y timesteps, if Input is 10 then we have konstant 10 y timesteps []
    #rescale = to make even lager timesteps, timestep in the end is scale*rescale []
    #T = Temperature, if no value is given 10°C is set as standart [°C]
    #
    #Output:
    #Ca = calcium concentration [mol/m³]
    #CO2 = CO2 concentration [atm]
    #timestep = time intervall between the entries of the time series of Ca and CO2 [y]
    #age = age distribution of the time series of Ca and CO2 [y]
    
    #determine the timesteps of the age Input (Sofular Speleothem)
    t=[]
    for i in range(0,len(age)-1):
        diff = age[i]-age[i+1]
        t.append(diff)
    timestep = np.array(t)
    
    #loop if all timesteps are already at the wanted sale
    if np.all(timestep==scale*rescale):
        age_new = np.arange(min(age),max(age)+scale*rescale,scale*rescale)
        age_new = age_new[::-1]
        co2_atmo_averaged_list = []
        for i in age_new:
            co2_mean = np.mean(CO2[(age_CO2>i-50)&(age_CO2<i+50)]) #average over 100 y
            co2_atmo_averaged_list.append(co2_mean)
        co2_atmo_averaged = np.array(co2_atmo_averaged_list)
        CO2_mean = np.array(co2_atmo_averaged)*10**-6*(1 - 0.61094*np.exp(17.625 * T / ( T + 243.04))/ 1013.25) #due to [ppm] in [atm], pCO2 [atm] = pCO2 [ppm] * (P_barometric pressure - P_waterpressure) / 1013.25, August-Roche-Magnus Equation and conversion from hPa in ATM          
        return Ca,CO2_mean,rescale*scale,age_new
        sys.exit()
        
    #get konstant Timesteps for Ca-Data
    i=0
    Ca_konstant = []
    while i in range(0,len(age)-2):
        if timestep[i] == scale:
            Ca_konstant.append(Ca[i])
            i = i+1
        elif (timestep[i]+timestep[i+1]) == scale:
            Ca_konstant.append((Ca[i] + Ca[i+1])/2)
            i = i+2
        else:
            sys.exit('wrong scale Parameter')
            
    #get same Timesteps for CO2 and average over 100 y
    age_new = np.arange(min(age),max(age)+scale,scale)
    age_new = age_new[::-1]
    co2_atmo_averaged_list = []
    for i in age_new:
        co2_mean = np.mean(CO2[(age_CO2>i-50)&(age_CO2<i+50)])
        co2_atmo_averaged_list.append(co2_mean)
    co2_atmo_averaged = np.array(co2_atmo_averaged_list)
    
    #rescale the timesteps to scale*rescale
    CO2_mean = []
    Ca_mean = []
    start = 0
    stop = rescale
    n = len(Ca_konstant)
    for i in range(0,int(n/rescale)):
        Ca_mean.append(np.mean(Ca_konstant[start:stop]))
        CO2_mean.append(np.mean(co2_atmo_averaged_list[start:stop]))
        start = stop+1
        stop = stop + rescale
    CO2_mean = np.array(CO2_mean)*10**-6*(1 - 0.61094*np.exp(17.625 * T / ( T + 243.04)) / 1013.25)  #due to [ppm] in [atm], pCO2 [atm] = pCO2 [ppm] * (P_barometric pressure - P_waterpressure) / 1013.25, August-Roche-Magnus Equation and conversion from hPa in ATM
    age_new_new = np.arange(min(age),max(age)-scale*rescale,scale*rescale)
    age_new_new = age_new_new[::-1]
    final_scale = rescale*scale
    return Ca_mean,CO2_mean,final_scale,age_new_new


# In[ ]:




