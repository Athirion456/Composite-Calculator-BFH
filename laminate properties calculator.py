# -*- coding: utf-8 -*-
"""
at1u16
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
    
def Function2(E1,E2,V12,G12,pthick,stse,Z,function):
   "Function 2 computes A, B and D matrix. Function 3 computes Exlam, Eylam, Vxylam and Gxylam. Z refers to the round factor you want applied on final matrix answer. "
#   Set the function outputs as global variable so they can be returned from outside the tkinter window. 
   global A
   global B
   global D
   global J
   global Exlam
   global Eylam
   global Vxylam
   global Gxylam
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
#       Calculate the A, B and D array using the relevant functions
#       A and B are converted to GN(mm^3)/m^2
       A = Aarray(A,d,ktable) * 10 ** 6
       B = Barray(B,d,ktable) * 10 ** 3
       D = Darray(D,d,ktable)
#       This sections join A,B and D together to form the full matrix
       E = np.concatenate((A,B),axis=1)
       F = np.concatenate((B,D),axis=1)
       J = np.concatenate((E, F))
#       The full matrix is formatted and printed
       K = np.round(J,Z)
       print('Full matrix =', K)
#       A, B, D and the full matrix are returned 
       return(A)
       return(B)
       return(D)
       return(J) 
   elif function == 3:
#       Calculates the A matrix and converts it to GN(mm^3)/m^2
       A = Aarray(A,d,ktable) * 10 ** 6
#       Computes the height of the stack
       h = sum(pthick) * (10 ** - 3)
#       Calculates the laminate elastic engineering properties
       Exlam = (1/h) * (A[0,0] - ((A[0,1] **2)/A[1,1]))
       Eylam = (1/h) * (A[1,1] - ((A[0,1] **2)/A[1,1]))
       Vxylam = A[0,1]/A[1,1]
       Gxylam = A[2,2]/h
#       Formats and print the laminate elastic engineering properties
       T1 = "Exlam = {Exlamf} GN/m^2".format(Exlamf = round((Exlam / 10 ** 9),Z))
       T2 = "Exlam = {Eylamf} GN/m^2".format(Eylamf = round((Eylam/ 10 ** 9),Z))
       T3 = "Vxylam = {Vxylamf}".format(Vxylamf = round(Vxylam,Z))
       T4 = "Gxylam = {Gxylamf} GN/m^2".format(Gxylamf = round((Gxylam/ 10 ** 9),Z))
       print(T1)
       print(T2)
       print(T3)
       print(T4)
#       Returns the laminate elastic engineering properties
       return(Exlam,Eylam,Vxylam,Gxylam) 
   else: 
#       Returns an error message
       return(print("Please choose a function type. Either 2 or 3."))
       
def Aarray(A,d,ktable):
#    Sets up a start value for indexing 
    ty = -1
#    Set up a list for each parameter in the relevant matrix
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
#    Sums each list and inputs them into the relevant matrix
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
#    Set up a list for each parameter in the relevant matrix
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
#    Sums each list and inputs them into the relevant matrix
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
#    Set up a list for each parameter in the relevant matrix
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
#    Sums each list and inputs them into the relevant matrix
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

#This section deals with the creation of a GUI.

import tkinter 
window = tkinter.Tk()
window.title("Input box")

#This section set up the input box and box labels for the tkinter window. 

tkinter.Label(window, text = "Longitudinal modulus (E1)").grid(row = 0) 
e1 = tkinter.Entry(window)
e1.grid(row = 0, column = 1)
tkinter.Label(window, text = "GPa").grid(row = 0, column = 2,sticky= tkinter.W)

tkinter.Label(window, text = "Transverse modulus (E2)").grid(row = 1)
e2 = tkinter.Entry(window)
e2.grid(row = 1, column = 1) 
tkinter.Label(window, text = "GPa").grid(row = 1, column = 2,sticky= tkinter.W)

tkinter.Label(window, text = "In-plane shear modulus (G12)").grid(row = 2) 
e3 = tkinter.Entry(window)
e3.grid(row = 2, column = 1)
tkinter.Label(window, text = "GPa").grid(row = 2, column = 2,sticky= tkinter.W)

tkinter.Label(window, text = "Major Poisson’s ratio (V12)").grid(row = 3)
e4 = tkinter.Entry(window)
e4.grid(row = 3, column = 1)

tkinter.Label(window, text = "Thickness sequences ").grid(row = 5)
e6 = tkinter.Entry(window)
e6.grid(row = 5, column = 1)
tkinter.Label(window, text = "Input in bottom to top layer order in mm with a comma between values i.e. 0.125, 0.25, 0.125, 0.25").grid(row = 5, column = 2, sticky= tkinter.W,)

tkinter.Label(window, text = "Stacking sequences").grid(row = 6)
e7 = tkinter.Entry(window)
e7.grid(row = 6, column = 1)
tkinter.Label(window, text = "Input in bottom to top layer order with a comma between values i.e. 90,0,90,0").grid(row = 6, column = 2, sticky= tkinter.W)

#Defines a default variable for the amount of decimal points in answer 
v = tkinter.StringVar(window, value='3')

tkinter.Label(window, text = "Amount of decimal points in answer").grid(row = 7)
e8 = tkinter.Entry(window, textvariable=v)
e8.grid(row = 7, column = 1)


#The function input is defined by which tkinter button you press, this is seen below 

# Function = 2
def test():
#    Extract the input from the input boxes.
    E11 = float(e1.get()) 
    E22 = float(e2.get())   
    G12 = float(e3.get()) 
    V12 = float(e4.get()) 
    Z = int(e8.get())
#    This section extracts the list inputs
#    First the string is split
#    Then it is converted to a integer and inputted into a list
    E6 = e6.get().split(",")
    E7 = e7.get().split(",")
    pthick = []
    stse = []
    for i in E7:
        stse.append(float(i))
    for i in E6:
        pthick.append(float(i))
#    Runs the function
    return(Function2(E11,E22,V12,G12,pthick,stse,Z,2))

# Function = 3
def test1():
#    Extract the input from the input boxes.
    E11 = float(e1.get()) 
    E22 = float(e2.get())   
    G12 = float(e3.get()) 
    V12 = float(e4.get()) 
    Z = int(e8.get())
#    This section extracts the list inputs
#    First the string is split
#    Then it is converted to a integer and inputted into a list
    E6 = e6.get().split(",")
    E7 = e7.get().split(",")
    pthick = []    
    stse = []
    for i in E7:
        stse.append(float(i))
    for i in E6:
        pthick.append(float(i))
#    Runs the function
    return(Function2(E11,E22,V12,G12,pthick,stse,Z,3))

#Set up the two button which run the function from the tkinter window    
B = tkinter.Button(text="Compute [A],[B] and [D]", command=test)
B.grid(row = 10, column = 1)
C = tkinter.Button(text="Compute Exlam,Eylam,Vxylam, Gxylam", command=test1)
C.grid(row = 11, column = 1)
window.mainloop()