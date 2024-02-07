# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 21:30:36 2020

@author: admin
"""

import numpy as np
import math

def Qbar(theta,Q11,Q22,Q12,Q66):
#    Calculates the m and n transformation matrix components.     
    m = math.cos(math.radians(theta))
    n = math.sin(math.radians(theta))
#    Defines the transformation matrix 
    T = np.matrix([[m**2, n**2, 2*m*n],[n**2, m**2, -2*m*n],[-m*n, m*n, (m**2 - n**2)]])
#    Manipulates the transformation matrix
    Ti = np.linalg.inv(T)
    Tt = np.transpose(Ti)

#    Defines the Q matrix
    Q = np.matrix([[Q11 , Q12, 0],[Q12, Q22, 0],[0, 0, Q66]])
#    Computes Qbar matrix and returns it 
    qbar = Ti * Q * Tt
    return(qbar)
    
def Function2(E1,E2,V12,G12,pthick,stse,function):
   "Function 2 computes A, B and D matrix. Function 3 computes Exlam, Eylam, Vxylam and Gxylam. Z refers to the round factor you want applied on final matrix answer. "
#   Set the function outputs as global variable so they can be returned from outside the tkinter window. 

#   Setups the variables used throughout the function.
   d = {}
   o = 0
   r =len(stse)+ 1    
#   Setups the arrays used throughout the function as arrays of 0 that will be added to.
   ktable = np.zeros((r,5))
   A = np.zeros((3,3))
   B = np.zeros((3,3))
   D = np.zeros((3,3))
#   Calculates V21
   V21 = (E2 * V12) / E1 
#   Calculates Q11, Q22, Q12 and Q66
   Q11 = E1 / (1 - (V12 * V21))
   Q22 = E2 / (1 - (V12 * V21))
   Q12 = (V12 * E2) / (1 - (V12 * V21))
   Q66 = G12 
#   Legacy code
#   Q11 = 20
#   Q22 = 2
#   Q12 = 1
#   Q66 = 1  
#   Generates a dictionary (python data type) of the Qbar matrix for each orientation. Using the Qbar functions. 
   for i in stse:
       o = 1 + o
       d["{0}".format(o)]= Qbar(i,Q11,Q22,Q12,Q66)
#   Calculates the positive max of the stacking sequences this is inputted into the first ktable cell.
   ktable[0,0] = sum(pthick)/2
#   Calculates the first row of the ktable
   for i in range(1,r):
       ktable[i,0] = ktable[i - 1,0] - pthick[i-1] 
       for i in range(len(stse)):
#           Calculates the second row of the ktable
           ktable[i,1] = ktable[i+1,0]  
#    Calculates the third, fourth and fifth row of the ktable 
   for i in range(r):
       ktable[i,2] = ktable[i,0] - ktable[i,1]
       ktable[i,3] = 0.5 * (((ktable[i,0]) ** 2)- (ktable[i,1]** 2 ))
       ktable[i,4] = (1/3) * ((ktable[i,0] ** 3 )- (ktable[i,1]** 3 ))
#  Ensures the all the columns but the first one in the last row are set to zero in the ktable.  
   ktable[r-1,1] = 0 
   ktable[r-1,2] = 0 
   ktable[r-1,3] = 0 
   ktable[r-1,4] = 0 
#   Using the function input to define which output will be produced. 
#   These could have been joined but this lead to too many outputs to be useful to the user. 
   if function == 2:
       #Removes any properties from core Qbar
       d['2']=np.zeros((3,3))
       print(d)
#       Calculate the A, B and D array using the relevant functions
#       A and B are converted to GN(mm^3)/m^2
       A = Aarray(A,d,ktable) 
       B = Barray(B,d,ktable) 
       D = Darray(D,d,ktable)
#       This sections join A,B and D together to form the full matrix
       #E = np.concatenate((A,B),axis=1)
       #F = np.concatenate((B,D),axis=1)
       #J = np.concatenate((E, F))
#       The full matrix is formatted and printed
#       A, B, D and the full matrix are returned 
       return(D)
   elif function == 3:
#       Calculates the A matrix and converts it to GN(mm^3)/m^2
       A = Aarray(A,d,ktable) 
#       Computes the height of the stack
       h = sum(pthick) 
#       Calculates the laminate elastic engineering properties
       Exlam = (1/h) * (A[0,0] - ((A[0,1] **2)/A[1,1]))
       Eylam = (1/h) * (A[1,1] - ((A[0,1] **2)/A[1,1]))
       Vxylam = A[0,1]/A[1,1]
       Gxylam = A[2,2]/h
#       Returns the laminate elastic engineering properties
       return(Exlam,Eylam,Vxylam,Gxylam) 
   else: 
#       Returns an error message
       return(print("Please choose a function type. Either 2 or 3."))

       
def Aarray(A,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    a = []
    a1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
#        Indexing   
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values
        a.append(ktable[ty,2] * l[0,0]) 
        a1.append(ktable[ty,2] * l[1,1])
        a2.append(ktable[ty,2] * l[2,2])
        a3.append(ktable[ty,2] * l[1,0])
        a4.append(ktable[ty,2] * l[2,0])
        a5.append(ktable[ty,2] * l[2,1])
#    Sums each list and inputs them into the matrix
    A[0,0] = sum(a)
    A[1,1] = sum(a1)
    A[2,2] = sum(a2)
    A[0,1] = sum(a3)
    A[0,2] = sum(a4)
    A[2,1] = sum(a5)
    A[1,0] = sum(a3)
    A[2,0] = sum(a4)
    A[1,2] = sum(a5)
    return(A)
    
def Barray(B,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    a = []
    a1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
#        Indexing 
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values 
        a.append(ktable[ty,3] * l[0,0]) 
        a1.append(ktable[ty,3] * l[1,1])
        a2.append(ktable[ty,3] * l[2,2])
        a3.append(ktable[ty,3] * l[1,0])
        a4.append(ktable[ty,3] * l[2,0])
        a5.append(ktable[ty,3] * l[2,1])
#    Sums each list and inputs them into the matrix
    B[0,0] = sum(a)
    B[1,1] = sum(a1)
    B[2,2] = sum(a2)
    B[0,1] = sum(a3)
    B[0,2] = sum(a4)
    B[2,1] = sum(a5)
    B[1,0] = sum(a3)
    B[2,0] = sum(a4)
    B[1,2] = sum(a5)
    return(B)
    
def Darray(D,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the matrix
    k = []
    k1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
#        Indexing 
        ty += 1 
#        Appends the value of the relevant multiplication of Qbar and ktable values 
        k.append(ktable[ty,4] * l[0,0]) 
        k1.append(ktable[ty,4] * l[1,1])
        a2.append(ktable[ty,4] * l[2,2])
        a3.append(ktable[ty,4] * l[1,0])
        a4.append(ktable[ty,4] * l[2,0])
        a5.append(ktable[ty,4] * l[2,1])
#    Sums each list and inputs them into the matrix
    D[0,0] = sum(k)
    D[1,1] = sum(k1)
    D[2,2] = sum(a2)
    D[0,1] = sum(a3)
    D[0,2] = sum(a4)
    D[2,1] = sum(a5)
    D[1,0] = sum(a3)
    D[2,0] = sum(a4)
    D[1,2] = sum(a5)
    return(D)
       
#Inputs
E1 = 54 # GPa
E2 = 18 # GPa
V12 = 0.28
G23 = 38 * 10**-3 # GPa
G12 = 6 # GPa
G13 = 38 * 10**-3 # GPa
pthick = [0.125,0.125,0.125,0.125]
stse = [0,90,90,0]
a = 750 # mm
b = 750 # mm
hc = 25 # mm

def sandwich(E1,E2,V12,G23,G13,G12,pthick,stse,a,b,hc):
    #Face sheet properties 
    Exlam,Eylam,Vxylam,Gxylam = Function2(E1,E2,V12,G12,pthick,stse,3)
    #Face sheet thickness
    P = sum(pthick)
    pthick1 = [P,hc,P]
    stse1 = [0,0,0]
    D = Function2(Exlam,Eylam,Vxylam,Gxylam,pthick1,stse1,2)
    #Alternative method 
    #Face sheet minor poisson ratio
    #V21 = (Eylam * Vxylam) / Exlam  
    #Calculates Q11, Q22, Q12 and Q66 of face sheets
    #Q11 = Exlam / (1 - (Vxylam * V21))
    #Q22 = Eylam / (1 - (Vxylam * V21))
    #Q12 = (Vxylam * Eylam) / (1 - (Vxylam * V21))
    #Q66 = Gxylam
    #Transforms Qbar matrix. No effect but included for sake of completeness 
    #theta = 0
    #Q = Qbar(theta,Q11,Q22,Q12,Q66)
    #Generates empty D matrix for sandwich structures 
    Ds = D
    print(Ds)
    #Alternative method 
    #Computes D matrix for sandwich structures
    #H = 2 * (((hc +P)/2)**2) *(P)
    #Ds = H * Q
    #Generates m,n table
    table = np.zeros((10,10))
    for m in range(1,10):
        for n in range(1,10):
            #Calculates Forces
            F11 = (((m**2)*(math.pi**2)*(Ds[0,0]))/(a**2)) + (((n**2)*(math.pi**2)*(Ds[2,2]))/(b**2)) + hc*G13
            F12 = (((m*n)*(math.pi**2)*(Ds[0,1] + Ds[2,2]))/(a*b))
            F13 = ((m*hc*G13)*(math.pi))/a
            F23 = ((n*hc*G23)*(math.pi))/b
            F22 = (((n**2)*(math.pi**2)*(Ds[1,1]))/(b**2)) + (((m**2)*(math.pi**2)*(Ds[2,2]))/(a**2)) + hc*G23
            #Unit conversion 
            F11 = F11 * 10 **-3
            F12 = F12 * 10 **-3
            F22 = F22 * 10 **-3
            #Calculates F33
            F33 = (((F11*(F23**2)) + (F22*(F13**2)) - (2*F12*F13*F23)) / ((F11*F22) - (F12**2)))
            No = (hc * (G13 + (((n/m)**2)*((a/b)**2)*G23))*(10**-3)) - (((a**2)/((math.pi **2)*(m**2))) * F33* (10**-6))
            #Method to calculate No without account for transverse shear. 
            #No = ((math.pi**2)/(a**2))*((Ds[0,0] * (m**2)) + (2*(Ds[0,1] + 2*Ds[2,2])*(n**2)) + (Ds[1,1]*((n**4)/(m**2))))
            #Convertes to MN/m
            No = No * 10**3
            #Appends to table
            table[n,m] = No 
    #Removes zero values from table
    GA = table[table != 0]
    #Finds minimum values   
    G = np.min(GA)
    #Finds n and m values of minimum values
    K = np.where(table == G)
    #Returns values and inputs for wrinkling calc
    return(table,G,K,Exlam, Vxylam)

table,G, K, Ef, vf = sandwich(E1,E2,V12,G23,G13,G12,pthick,stse,a,b,hc)   

#Wrinkling
#Q = 0.825
#Alternate Q value
Q = 0.5
Ec = 110 #MPa
Ef = Ef * 10**3  #MPa
Gc = 38  #MPa
#Wrinkling calc
Wrinkling = (Q) * (((Ec*Ef*Gc)/(1-(vf**2)))** (1/3))
