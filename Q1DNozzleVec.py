# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 09:51:08 2019

Anderson Chapter 7.1
Numerical Solutions of quasi-1d nozzle flows Vectorized

@author: @jayjoeay
"""
import numpy as np
import matplotlib.pyplot as plt

#setting up initial conditions
gamma=1.4
x=np.arange(0,3.1,0.01)
dx=0.1
A=1+2.2*(x-1.5)**2
rho=1-0.3146*x
T=1-0.2314*x
V=(0.1+1.09*x)*T**0.5
M=V/T**0.5

#Calculation of CFL time step
dt=np.min(0.5*dx/(T[1:]**0.5+V[1:]))
size=15000

#Variable preDef
rho_p=np.zeros_like(rho)        #predicted values
drhodt_p=np.zeros_like(rho)

V_p=np.zeros_like(V)
dVdt_p=np.zeros_like(V)

T_p=np.zeros_like(T)
dTdt_p=np.zeros_like(T)

drhodt_c=np.zeros_like(rho)     #Corrected Values
dVdt_c=np.zeros_like(V)
dTdt_c=np.zeros_like(T)

T_inTime=np.zeros(size)         #Used for plotting 
rho_inTime=np.zeros(size)
V_inTime=np.zeros(size)
P_inTime=np.zeros(size)
M_inTime=np.zeros(size)

err=1
iterNo=0
throat=np.where(A==1)
throat=throat[0]


while err>0.01 :
    iterNo=iterNo+1
    rho_n=rho.copy()
    V_n=V.copy()
    T_n=T.copy()
    M_n=M.copy()
    

    drhodt_p[:-1]=-rho[:-1]*(V[1:]-V[:-1])/dx-rho[:-1]*V[:-1]*(np.log(A[1:])-np.log(A[:-1]))/dx - V[:-1]*(rho[1:]-rho[:-1])/dx
    rho_p=rho+drhodt_p*dt
    

    dVdt_p[:-1]=-V[:-1]*(V[1:]-V[:-1])/dx - (1/gamma)*((T[1:]-T[:-1])/dx + (T[:-1]/rho[:-1])*(rho[1:]-rho[:-1])/dx)
    V_p=V+dVdt_p*dt
    
    
    dTdt_p[:-1]=-V[:-1]*(T[1:]-T[:-1])/dx - (gamma-1)*T[:-1]*( (V[1:]-V[:-1])/dx + V[:-1]*(np.log(A[1:])-np.log(A[:-1]))/dx)
    T_p=T+dTdt_p*dt
    


    drhodt_c[1:]=-rho_p[1:]*(V_p[1:]-V_p[:-1])/dx - rho_p[1:]*V_p[1:]*(np.log(A[1:])-np.log(A[:-1]))/dx - V_p[1:]*(rho_p[1:]-rho_p[:-1])/dx
    drhodt_av=0.5*(drhodt_c+drhodt_p)
    rho=rho+drhodt_av*dt
    rho_inTime[iterNo]=rho[throat]
    
    

    dVdt_c[1:]=-V_p[1:]*(V_p[1:]-V_p[:-1])/dx - (1/gamma)*((T_p[1:]-T_p[:-1])/dx + (T_p[1:]/rho_p[1:])*(rho_p[1:]-rho_p[:-1])/dx)
    dVdt_av=0.5*(dVdt_c+dVdt_p)
    V=V+dVdt_av*dt
    V_inTime[iterNo]=V[throat]
    

    dTdt_c[1:]=-V_p[1:]*(T_p[1:]-T_p[:-1])/dx - (gamma-1)*T_p[1:]*( (V_p[1:]-V_p[:-1])/dx + V_p[1:]*(np.log(A[1:])-np.log(A[:-1]))/dx)
    dTdt_av=0.5*(dTdt_c+dTdt_p)
    T=T+dTdt_av*dt
    T_inTime[iterNo]=T[throat]
    
    rho[0]=1
    rho[-1]=2*rho[-2]-rho[-3]
    V[0]=2*V[1]-V[2]
    V[-1]=2*V[-2]-V[-3]
    T[0]=1
    T[-1]=2*T[-2]-T[-3]
    
    P=rho*T
    P_inTime[iterNo]=P[throat]
    
    M=V/T**0.5
    M_inTime[iterNo]=M[throat]
    err=np.sum(np.abs(M-M_n))
    
tt=np.arange(iterNo)













