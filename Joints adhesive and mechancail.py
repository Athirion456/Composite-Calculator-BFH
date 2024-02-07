# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Mechancail fastening
#stress conentration factors 

def stress_conc_edge_hole(e,w,d):
    "e distances from edge, w width strip, d diameter of hole"
    if e/w < 1:
        o = (w/e - 1)
    else:
        o = 0
    kte = (w/d) + (d/w) + (0.5*(1-(d/w))*o)
    return(kte)
    
def stress_conc_unloaded_hole(e,w,d):    
    "e distances from edge, w width, d diameter of holes"
    kte = 2 + ((1-d/w)**3)
    return(kte)

def stress_conc_row_of_holes(e,p,d):    
    "e distances from edge, p pitch, d diameter of holes"
    kte = 1 + (2*(1-d/p)**1.5)
    return(kte)    

def stress_conc_failure(Ftu,w,d,P):
    "Ftu unnotched ultimate tensile strength, w width, d diameter of holes, P failing load"
    kte = (Ftu * (w-d))/P
    return(kte)
    
def stress_conc_edge_loaded_hole(e,p,d):
    "e distances from edge, p pitch, d diameter of hole"
    if e/p < 1:
        o = (p/e - 1)
    else:
        o = 0
    kte = (p/d) + (0.5*(1-(d/p))*o)
    return(kte)   
    
#Adhesive joints 
#Bigwood and Crocombe

#G,E1 mech properties
#Poisson ratio   
#t = thickness
#n = adhesvies thickness
    
alpha = (G*(1-(v1 ** 2)))/(E1*t1*n)

beta = (12*E1*(1-v1 ** 2))/(E1*n*(t1**3))

def T_Loading(T,t1,t2,alpha1,alpha2):
    "Alpha 1 top surface, Alpha 2 bottom surface, T tensile load"
    t = (alpha1*T)/((alpha1 + alpha2)**0.5)
    print(t)
    return(t)
    
def V_Loading(V,t1,t2,beta1,beta2): 
    "Beta 1 top surface Beta 2 bottom surface, V through thickness load"
    sigma = -((beta1*V*(2**0.5))/((beta1 + beta2)**0.75))
    t = (3*V)/(4*t1)
    print(t,sigma)
    return(t,sigma)
    
def M_Loading(M,t1,t2, alpha1, alpha2, beta1, beta2):
    "Beta/Alpha 1 top surface Beta/Alpha 2 bottom surface, M moment load"    
    sigma = -((beta1*M)/((beta1 + beta2)**0.5))
    t = (3*alpha1*M)/(t1*((alpha1 + alpha2)**0.5))
    print(t,sigma)
    return(t,sigma)