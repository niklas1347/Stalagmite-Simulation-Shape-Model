#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from auxiliary_functions import * #to read in my self definded functions


# ### define Functions

# In[2]:


def Simulator_comparsion_eq_flow_equalaxis(TC, tau, pCO2cave, Ca_drip, timestep=1):
    # This Code shall provide an Algorithm to simulate the shape of a speleothem 
    # using the FLOW-Modell as statet in: Modeling stalagmite growth 
    # by first principles of chemistry and physics of calcite precipitation, Kaufmann, 2008"
    # and starts growing in equilibrium and has equal stepsize at the x- and y-axis, 
    # no output values 
    #
    # Input Variables:
    # TC: Temperatur [°C]
    # tau: Drip Rate [s]
    # pCO2cave: CO2 concentration in the Cave [atm]
    # Ca_drip: Ca concentration of the water drop [mol/m^3]
    # timestep: defines how many years pass between a layer, default is 1 [y]
    
    # Expansion of all time series by 250 times the first value
    TC_before = []
    tau_before = []
    pCO2cave_before = []
    Ca_drip_before = []
    for i in range(0,250):
        TC_before.append(TC[0])
        tau_before.append(tau[0])
        pCO2cave_before.append(pCO2cave[0])
        Ca_drip_before.append(Ca_drip[0])

    TC = np.concatenate((TC_before,TC),axis=0)
    tau = np.concatenate((tau_before,tau),axis=0)
    pCO2cave = np.concatenate((pCO2cave_before,pCO2cave),axis=0)
    Ca_drip = np.concatenate((Ca_drip_before,Ca_drip),axis=0)

    TK =  TC + 273.15
    
    #fixed Parameters
    V = 1e-7       # water volume of the drop [m³], after "Processes in Karst Systems: Physics, Chem- istry, and Geology (Band 4); Dreybrodth; 1988"
    delta = 1e-4   # water film thickness [m]
    
    
    alpha = alpha_poly(TC,delta)
    
    # Calculation of the excess calzit concentration, if the CO2 concentration in the Cave is higher than in the drop there is no growth,
    # calculations after "Regular Satalgmites: The theory behind their shape; Dreybrodt; 2008"
    CaEx = np.zeros(len(TC))
    for i in range(0,len(TC)):
        tmp_cave = Konstanten(TC[i], pCO2cave[i])
 
        if Ca_drip[i] > (tmp_cave / np.sqrt(0.8)):
            CaEx[i] = ( Ca_drip[i] - tmp_cave  / np.sqrt(0.8) ) #excess calcium, apparent after "Chemical kinetics, speleothem growth and climate, Dreybrodt, 1999" [mol/m³]
        else:
            CaEx[i] = 0
    
    Z=delta/alpha

    # Equilibrium radius
    R=np.sqrt(V/(np.pi*alpha*tau))

    #Growth
    W_0 = 1.168*10**3 * alpha * CaEx #[m/y]
    Growth = W_0
    
    #inital grid points
    N=5000 #defines the resolution of the stalagmite and the numbers of the polygon
    l=3*R/N   #defines the length of the Polygon, has just visiual influences
    
    X=np.zeros(shape=(len(tau)+1,N))
    Y=np.zeros(shape=(len(tau)+1,N))

    X[0,] = np.linspace(0,max(l)*N,N)
    Y[0,] = 0
    
    
    # Growing the Stalgmite
    for i in range(0,len(tau)):
        #Calcualtion of the growth of the surface
        g = np.zeros(N)
        g[0] = W_0[i]*timestep

        beta = np.zeros(N)
        L = np.zeros(N)

        # Calculation of the next layer
        X[i+1,0]=0
        Y[i+1,0]=Y[i,0]+g[0]

        for j in range(1,N):
            L[j]=np.sqrt((X[i,j]-X[i,j-1])**2+(Y[i,j]-Y[i,j-1])**2)
            g[j]=g[j-1]*(1-2*L[j]*X[i,j-1]/(R[i]**2)) #nach "Regular Satalgmites: The theory behind their shape; Dreybrodt; 2008"
            beta=np.arctan((Y[i,j-1]-Y[i,j])/(X[i,j]-X[i,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i,j]+g[j]*np.sin(beta)
            Y[i+1,j]=Y[i,j]+g[j]*np.cos(beta)

        #Putting the grid points at the same distance again
        for j in range(1,N):
            beta=np.arctan((Y[i+1,j-1]-Y[i+1,j])/(X[i+1,j]-X[i+1,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i+1,j-1]+max(l)*np.cos(beta)
            Y[i+1,j]=Y[i+1,j-1]-max(l)*np.sin(beta)

        Y[Y<0]=0

    # Reversal of the expansion
    X=X[250:]
    Y=Y[250:]-max(Y[250])
    #Calculation of the depth:
    
    #Calculation of the depth:
    Depth = [0]
    for i in range(1,len(W_0[250:])+1):
        Depth.append(Depth[i-1] + W_0[250:][i-1]*timestep)
    Depth = Depth[1:]
    Depth = Depth[-1]-Depth
    
    fig = plt.figure(num=1,figsize=(9,6),dpi=300)
    ax = fig.add_subplot(111)
    plt.scatter(X,Y,marker=".",color='royalblue',s=1, linewidths=0)
    plt.scatter(-X,Y,marker=".",color='royalblue',s=1, linewidths=0)
    
    #print every 1000 years a red layer
    i=10
    while i in range(0,len(Depth)):
        plt.scatter(X[i],Y[i],marker=".",color='indigo',s=1, linewidths=0)
        plt.scatter(-X[i],Y[i],marker=".",color='indigo',s=1, linewidths=0)
        i=i+10
    
    plt.yticks(np.arange(0, 1.42+.1, .1))
    plt.ylim(0,1.44)
    ax.yaxis.set_minor_locator(MultipleLocator(0.02))# defines setting of the small ticks
    loc = plticker.MultipleLocator(base=.05) # this locator puts ticks at regular intervals
    plt.xticks(np.arange(-.4, .45, .1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.02))# defines setting of the small ticks
    plt.xlim(-0.15,0.15)
    ax.set_aspect('equal', 'box')
 
    plt.xlabel(r'Radius [m]',size="large")
    plt.ylabel(r'Heigth [m]',size="large")
    
    return


# In[3]:


def Simulator_comparsion_flow_upper(TC, tau, pCO2cave, Ca_drip, timestep=1):
    # This Code shall provide an Algorithm to simulate the shape of a speleothem 
    # using the FLOW-Modell as statet in: Modeling stalagmite growth 
    # by first principles of chemistry and physics of calcite precipitation, Kaufmann, 2008"
    # and starts growing in equilibrium and just displays the upper layer, 
    # no output values 
    #
    # Input Variables:
    # TC: Temperatur [°C]
    # tau: Drip Rate [s]
    # pCO2cave: CO2 concentration in the Cave [atm]
    # Ca_drip: Ca concentration of the water drop [mol/m^3]
    # timestep: defines how many years pass between a layer, default is 1 [y]
    
    # Expansion of all time series by 250 times the first value
    TC_before = []
    tau_before = []
    pCO2cave_before = []
    Ca_drip_before = []
    for i in range(0,250):
        TC_before.append(TC[0])
        tau_before.append(tau[0])
        pCO2cave_before.append(pCO2cave[0])
        Ca_drip_before.append(Ca_drip[0])

    TC = np.concatenate((TC_before,TC),axis=0)
    tau = np.concatenate((tau_before,tau),axis=0)
    pCO2cave = np.concatenate((pCO2cave_before,pCO2cave),axis=0)
    Ca_drip = np.concatenate((Ca_drip_before,Ca_drip),axis=0)

    TK =  TC + 273.15
    
    #fixed Parameters
    V = 1e-7       # water volume of the drop [m³], after "Processes in Karst Systems: Physics, Chem- istry, and Geology (Band 4); Dreybrodth; 1988"
    delta = 1e-4   # water film thickness [m]
    
    
    alpha = alpha_poly(TC,delta)
    
    # Calculation of the excess calzit concentration, if the CO2 concentration in the Cave is higher than in the drop there is no growth,
    # calculations after "Regular Satalgmites: The theory behind their shape; Dreybrodt; 2008"
    CaEx = np.zeros(len(TC))
    for i in range(0,len(TC)):
        tmp_cave = Konstanten(TC[i], pCO2cave[i])
 
        if Ca_drip[i] > (tmp_cave / np.sqrt(0.8)):
            CaEx[i] = ( Ca_drip[i] - tmp_cave  / np.sqrt(0.8) ) #excess calcium, apparent after "Chemical kinetics, speleothem growth and climate, Dreybrodt, 1999" [mol/m³]
        else:
            CaEx[i] = 0
    
    Z=delta/alpha

    # Equilibrium radius
    R=np.sqrt(V/(np.pi*alpha*tau))

    #Growth
    W_0 = 1.168*10**3 * alpha * CaEx #[m/y]
    Growth = W_0
    
    #inital grid points
    N=2000 #defines the resolution of the stalagmite and the numbers of the polygon
    l=5*R/N   #defines the length of the Polygon, has just visiual influences
    
    X=np.zeros(shape=(len(tau)+1,N))
    Y=np.zeros(shape=(len(tau)+1,N))

    X[0,] = np.linspace(0,max(l)*N,N)
    Y[0,] = 0
    
    
    # Growing the Stalgmite
    for i in range(0,len(tau)):
        #Calcualtion of the growth of the surface
        g = np.zeros(N)
        g[0] = W_0[i]*timestep

        beta = np.zeros(N)
        L = np.zeros(N)

        # Calculation of the next layer
        X[i+1,0]=0
        Y[i+1,0]=Y[i,0]+g[0]

        for j in range(1,N):
            L[j]=np.sqrt((X[i,j]-X[i,j-1])**2+(Y[i,j]-Y[i,j-1])**2)
            g[j]=g[j-1]*(1-2*L[j]*X[i,j-1]/(R[i]**2)) #nach "Regular Satalgmites: The theory behind their shape; Dreybrodt; 2008"
            beta=np.arctan((Y[i,j-1]-Y[i,j])/(X[i,j]-X[i,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i,j]+g[j]*np.sin(beta)
            Y[i+1,j]=Y[i,j]+g[j]*np.cos(beta)

        #Putting the grid points at the same distance again
        for j in range(1,N):
            beta=np.arctan((Y[i+1,j-1]-Y[i+1,j])/(X[i+1,j]-X[i+1,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i+1,j-1]+max(l)*np.cos(beta)
            Y[i+1,j]=Y[i+1,j-1]-max(l)*np.sin(beta)

        Y[Y<0]=0

    # Reversal of the expansion
    X=X[250:]
    Y=Y[250:]-max(Y[250])
    
    
    fig = plt.figure(num=1,figsize=(3,3),dpi=300)
     
    ax = fig.add_subplot(111)
    plt.scatter(X[-1],Y[-1],marker=".",color='royalblue',s=1, linewidths=0)
    plt.scatter(-X[-1],Y[-1],marker=".",color='royalblue',s=1, linewidths=0)
    plt.ylim(1.25, 1.40)
    plt.xlim(-0.1,0.1)
    ax.set_aspect('equal', 'box')
    plt.axis('off')
    plt.text(0.80, 0.97, 'Flow',transform=plt.gca().transAxes,fontsize=8,verticalalignment = 'top',horizontalalignment = 'right',color='royalblue')
    
    return


# In[4]:


def Simulator_comparsion_gauss_upper(TC, tau, pCO2cave, Ca_drip,phi=1, timestep=1):
    # This Code shall provide an Algorithm to simulate the shape of a speleothem 
    # using the GAUSS-Modell as statet in: Modeling stalagmite growth 
    # by first principles of chemistry and physics of calcite precipitation, Kaufmann, 2008"
    # and starts growing in equilibrium and just displays the upper layer, 
    # no output values 
    #
    # Input Variables:
    # TC: Temperatur [°C]
    # tau: Drip Rate [s]
    # pCO2cave: CO2 concentration in the Cave [atm]
    # Ca_drip: Ca concentration of the water drop [mol/m^3]
    # timestep: defines how many years pass between a layer, default is 1 [y]
    
    # Expansion of all time series by 250 times the first value
    TC_before = []
    tau_before = []
    pCO2cave_before = []
    Ca_drip_before = []
    for i in range(0,250):
        TC_before.append(TC[0])
        tau_before.append(tau[0])
        pCO2cave_before.append(pCO2cave[0])
        Ca_drip_before.append(Ca_drip[0])

    TC = np.concatenate((TC_before,TC),axis=0)
    tau = np.concatenate((tau_before,tau),axis=0)
    pCO2cave = np.concatenate((pCO2cave_before,pCO2cave),axis=0)
    Ca_drip = np.concatenate((Ca_drip_before,Ca_drip),axis=0)

    TK =  TC + 273.15
    
    #fixed Parameters
    V = 1e-7       # water volume of the drop [m³], after "Processes in Karst Systems: Physics, Chem- istry, and Geology (Band 4); Dreybrodth; 1988"
    delta = 1e-4   # water film thickness [m]
    
    
    alpha = alpha_poly(TC,delta)
    
    # Calculation of the excess calzit concentration, if the CO2 concentration in the Cave is higher than in the drop there is no growth,
    # calculations after "Regular Satalgmites: The theory behind their shape; Dreybrodt; 2008"
    CaEx = np.zeros(len(TC))
    for i in range(0,len(TC)):
        tmp_cave = Konstanten(TC[i], pCO2cave[i])
 
        if Ca_drip[i] > (tmp_cave / np.sqrt(0.8)):
            CaEx[i] = ( Ca_drip[i] - tmp_cave  / np.sqrt(0.8) ) #excess calcium, apparent after "Chemical kinetics, speleothem growth and climate, Dreybrodt, 1999" [mol/m³]
        else:
            CaEx[i] = 0
    
    Z=delta/alpha
    
    #Mixing prozess like "Modelling stalagmite growth and δ 13C as a function of drip interval and temperature; Mühlinghaus; 2007"
    #We take lim n -> inf because we look at big Time arrays
    Lambda = phi/(1 - (1-phi)*np.exp(-tau*alpha/delta))

    # Equilibrium radius
    R=np.sqrt(V/(np.pi*delta*(1-np.exp(-tau/Z))))/ np.sqrt(Lambda) 

    #Growth
    W_0 = 1.168 * 10**3 * CaEx * delta  / tau * ( 1 - np.exp( - tau / Z)) * Lambda  #[m/y]
    Growth = W_0
    
    #inital grid points
    N=2000 #defines the resolution of the stalagmite and the numbers of the polygon
    l=5*R/N   #defines the length of the Polygon, has just visiual influences
    
    X=np.zeros(shape=(len(tau)+1,N))
    Y=np.zeros(shape=(len(tau)+1,N))

    X[0,] = np.linspace(0,max(l)*N,N)
    Y[0,] = 0
    
    
    # Growing the Stalgmite
    for i in range(0,len(tau)):
        #Calcualtion of the growth of the surface
        g = np.zeros(N)
        g[0] = W_0[i]*timestep

        beta = np.zeros(N)
        L = np.zeros(N)

        # Calculation of the next layer
        X[i+1,0]=0
        Y[i+1,0]=Y[i,0]+g[0]

        for j in range(1,N):
            L[j]=np.sqrt((X[i,j]-X[i,j-1])**2+(Y[i,j]-Y[i,j-1])**2)
            g[j] = g[0]*np.exp(-sum(L)**2 * np.pi / (4 * R[i]**2))
            beta=np.arctan((Y[i,j-1]-Y[i,j])/(X[i,j]-X[i,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i,j]+g[j]*np.sin(beta)
            Y[i+1,j]=Y[i,j]+g[j]*np.cos(beta)

        #Putting the grid points at the same distance again
        for j in range(1,N):
            beta=np.arctan((Y[i+1,j-1]-Y[i+1,j])/(X[i+1,j]-X[i+1,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i+1,j-1]+max(l)*np.cos(beta)
            Y[i+1,j]=Y[i+1,j-1]-max(l)*np.sin(beta)

        Y[Y<0]=0

    # Reversal of the expansion
    X=X[250:]
    Y=Y[250:]-max(Y[250])
    
    
    fig = plt.figure(num=1,figsize=(3,3),dpi=300)
     
    ax = fig.add_subplot(111)
    plt.scatter(X[-1],Y[-1],marker=".",color='coral',s=1, linewidths=0)
    plt.scatter(-X[-1],Y[-1],marker=".",color='coral',s=1, linewidths=0)
    plt.ylim(1.185, 1.325)
    plt.xlim(-0.1,.1)
    ax.set_aspect('equal', 'box')
    plt.axis('off')
    plt.text(0.80, 0.97, 'Gauss',transform=plt.gca().transAxes,fontsize=8,verticalalignment = 'top',horizontalalignment = 'right',color='coral')
    
    return


# In[5]:


def Simulator_comparsion_exp_upper(TC, tau, pCO2cave, Ca_drip,phi=1, timestep=1):
    # This Code shall provide an Algorithm to simulate the shape of a speleothem 
    # using the EXP-Modell as statet in: Modeling stalagmite growth 
    # by first principles of chemistry and physics of calcite precipitation, Kaufmann, 2008"
    # and starts growing in equilibrium and just displays the upper layer, 
    # no output values 
    #
    # Input Variables:
    # TC: Temperatur [°C]
    # tau: Drip Rate [s]
    # pCO2cave: CO2 concentration in the Cave [atm]
    # Ca_drip: Ca concentration of the water drop [mol/m^3]
    # timestep: defines how many years pass between a layer, default is 1 [y]
    
    # Expansion of all time series by 250 times the first value
    TC_before = []
    tau_before = []
    pCO2cave_before = []
    Ca_drip_before = []
    for i in range(0,250):
        TC_before.append(TC[0])
        tau_before.append(tau[0])
        pCO2cave_before.append(pCO2cave[0])
        Ca_drip_before.append(Ca_drip[0])

    TC = np.concatenate((TC_before,TC),axis=0)
    tau = np.concatenate((tau_before,tau),axis=0)
    pCO2cave = np.concatenate((pCO2cave_before,pCO2cave),axis=0)
    Ca_drip = np.concatenate((Ca_drip_before,Ca_drip),axis=0)

    TK =  TC + 273.15
    
    #fixed Parameters
    V = 1e-7       # water volume of the drop [m³], after "Processes in Karst Systems: Physics, Chem- istry, and Geology (Band 4); Dreybrodth; 1988"
    delta = 1e-4   # water film thickness [m]
    
    
    alpha = alpha_poly(TC,delta)
    
    # Calculation of the excess calzit concentration, if the CO2 concentration in the Cave is higher than in the drop there is no growth,
    # calculations after "Regular Satalgmites: The theory behind their shape; Dreybrodt; 2008"
    CaEx = np.zeros(len(TC))
    for i in range(0,len(TC)):
        tmp_cave = Konstanten(TC[i], pCO2cave[i])
 
        if Ca_drip[i] > (tmp_cave / np.sqrt(0.8)):
            CaEx[i] = ( Ca_drip[i] - tmp_cave  / np.sqrt(0.8) ) #excess calcium, apparent after "Chemical kinetics, speleothem growth and climate, Dreybrodt, 1999" [mol/m³]
        else:
            CaEx[i] = 0
    
    Z=delta/alpha
    
    #Mixing prozess like "Modelling stalagmite growth and δ 13C as a function of drip interval and temperature; Mühlinghaus; 2007"
    #We take lim n -> inf because we look at big Time arrays
    Lambda = phi/(1 - (1-phi)*np.exp(-tau*alpha/delta))

    # Equilibrium radius
    R=np.sqrt(V/(np.pi*delta*(1-np.exp(-tau/Z))))/ np.sqrt(Lambda) 

    #Growth
    W_0 = 1.168 * 10**3 * CaEx * delta  / tau * ( 1 - np.exp( - tau / Z)) * Lambda  #[m/y]
    Growth = W_0
    
    #inital grid points
    N=3000 #defines the resolution of the stalagmite and the numbers of the polygon
    l=7*R/N   #defines the length of the Polygon, has just visiual influences
    
    X=np.zeros(shape=(len(tau)+1,N))
    Y=np.zeros(shape=(len(tau)+1,N))

    X[0,] = np.linspace(0,max(l)*N,N)
    Y[0,] = 0
    
    
    # Growing the Stalgmite
    for i in range(0,len(tau)):
        #Calcualtion of the growth of the surface
        g = np.zeros(N)
        g[0] = W_0[i]*timestep

        beta = np.zeros(N)
        L = np.zeros(N)

        # Calculation of the next layer
        X[i+1,0]=0
        Y[i+1,0]=Y[i,0]+g[0]

        for j in range(1,N):
            L[j]=np.sqrt((X[i,j]-X[i,j-1])**2+(Y[i,j]-Y[i,j-1])**2)
            g[j] = g[0] * np.exp( - np.sum(L) / R[i])
            beta=np.arctan((Y[i,j-1]-Y[i,j])/(X[i,j]-X[i,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i,j]+g[j]*np.sin(beta)
            Y[i+1,j]=Y[i,j]+g[j]*np.cos(beta)

        #Putting the grid points at the same distance again
        for j in range(1,N):
            beta=np.arctan((Y[i+1,j-1]-Y[i+1,j])/(X[i+1,j]-X[i+1,j-1]))
            if beta <0 :
                beta = beta +np.pi
            X[i+1,j]=X[i+1,j-1]+max(l)*np.cos(beta)
            Y[i+1,j]=Y[i+1,j-1]-max(l)*np.sin(beta)

        Y[Y<0]=0

    # Reversal of the expansion
    X=X[250:]
    Y=Y[250:]-max(Y[250])
    
    
    fig = plt.figure(num=1,figsize=(3,3),dpi=300)
     
    ax = fig.add_subplot(111)
    plt.scatter(X[-1],Y[-1],marker=".",color='forestgreen',s=1, linewidths=0)
    plt.scatter(-X[-1],Y[-1],marker=".",color='forestgreen',s=1, linewidths=0)
    plt.ylim(1.185, 1.325)
    plt.xlim(-0.1,.1)
    ax.set_aspect('equal', 'box')
    plt.axis('off')
    plt.text(0.80, 0.97, 'Exp',transform=plt.gca().transAxes,fontsize=8,verticalalignment = 'top',horizontalalignment = 'right',color='forestgreen')
    
    return


# In[ ]:




