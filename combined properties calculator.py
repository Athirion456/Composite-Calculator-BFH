# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 00:18:13 2020

@author: admin
"""
import numpy as np
import math

def Function3 (Ef, Ff, Rhof, Nuf, Gf, Em, Rhom, Num, Gm, method, fraction, stse, pthick, Z):
   "method M for rule of mixture and H for halpin tsai,fraction W for weight fraction and V for volume fraction"
   if fraction == "W":
            Vf = ((Ff / Rhof) / ((Ff / Rhof) + ((1-Ff) / Rhom)))
   elif fraction == "V":
            Vf = Ff 
   else:
        return(print("Incorrect fraction input"))
   Vm = 1 - Vf 
   if method == "M":
        "Rule of mixturs"
        E11 = (Ef * Vf) + (Em * Vm) 
        E22 = (Ef * Em)/ ((Ef * Vm) + (Em * Vf))
        V12 = (Vf * Nuf) +  (Vm * Num) 
        V21 = V12 * (E22 / E11)
        if (Gf * Gm) == 1:
            Gf = Ef / (2 * ( 1 + Nuf))
            Gm = Em / (2 * ( 1 + Num))
            G12 =(Gf * Gm) / ((Gf * Vm) + (Gm * Vf))
        else: 
            G12 =(Gf * Gm) / ((Gf * Vm) + (Gm * Vf))
   elif method == "H":
        "Halpin Tsai"
        E11 = (Ef * Vf) + (Em * Vm) 
        n1 = ((Ef / Em) - 1)/( (Ef / Em) + 0.2)
        E22 = Em * ((1 + (0.2 * n1 * Vf))/(1 - (n1 * Vf)))
        V12 = (Vf * Nuf) +  (Vm * Num) 
        V21 = V12 * (E22 / E11)
        if (Gf * Gm) == 1:
                 Gf = Ef / (2 * ( 1 + Nuf))
                 Gm = Em / (2 * ( 1 + Num))
                 n2 = ((Gf / Gm) - 1)/( (Gf / Gm) + 1)
                 G12 =Gm * ((1 + (1 * n2 * Vf))/(1 - (n1 * Vf)))
        else: 
                 n2 = ((Gf / Gm) - 1)/( (Gf / Gm) + 1)
                 G12 =Gm * ((1 + (1 * n2 * Vf))/(1 - (n1 * Vf)))
   else:
        return(print("Incorrect method input"))    
   d = {}
   o = 0
   r =len(stse)+ 1    
   ktable = np.zeros((r,5))
   A = np.zeros((3,3))
   B = np.zeros((3,3))
   D = np.zeros((3,3))
   Q11 = E11 / (1 - (V12 * V21))
   Q22 = E22 / (1 - (V12 * V21))
   Q12 = (V12 * E22) / (1 - (V12 * V21))
   Q66 = G12 
   #Q11 = 20
   #Q22 = 2
   #Q12 = 1
   #Q66 = 1  
   for i in stse:
       o = 1 + o
       d["{0}".format(o)]= Qbar(i,Q11,Q22,Q12,Q66)
   ktable[0,0] = sum(pthick)/2
   for i in range(1,r):
       ktable[i,0] = ktable[i - 1,0] - pthick[i-1] 
       for i in range(len(stse)):
           ktable[i,1] = ktable[i+1,0]  
   for i in range(r):
       ktable[i,2] = ktable[i,0] - ktable[i,1]
       ktable[i,3] = 0.5 * (((ktable[i,0]) ** 2)- (ktable[i,1]** 2 ))
       ktable[i,4] = (1/3) * ((ktable[i,0] ** 3 )- (ktable[i,1]** 3 ))
   ktable[r-1,1] = 0 
   ktable[r-1,2] = 0 
   ktable[r-1,3] = 0 
   ktable[r-1,4] = 0 
   A = Aarray(A,d,ktable)
   B = Barray(B,d,ktable)
   D = Darray(D,d,ktable)
   E = np.concatenate((A,B),axis=1)
   F = np.concatenate((B,D),axis=1)
   J = np.concatenate((E, F))
   global K
   K = np.round(J,Z)
   #print(ktable)
   print(K)
   return(K)
   
def Qbar(theta,Q11,Q22,Q12,Q66):    
    m = math.cos(math.radians(theta))
    n = math.sin(math.radians(theta))
    T = np.matrix([[m**2, n**2, 2*m*n],[n**2, m**2, -2*m*n],[-m*n, m*n, (m**2 - n**2)]])
    Ti = np.linalg.inv(T)
    Tt = np.transpose(Ti)
    Q = np.matrix([[Q11 , Q12, 0],[Q12, Q22, 0],[0, 0, Q66]])
    qbar = Ti * Q * Tt
    return(qbar)
    
def Aarray(A,d,ktable):
    ty = -1
    a = []
    a1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
        ty += 1 
        a.append(ktable[ty,2] * l[0,0]) 
        a1.append(ktable[ty,2] * l[1,1])
        a2.append(ktable[ty,2] * l[2,2])
        a3.append(ktable[ty,2] * l[1,0])
        a4.append(ktable[ty,2] * l[2,0])
        a5.append(ktable[ty,2] * l[2,1])
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
    ty = -1
    a = []
    a1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
        ty += 1 
        a.append(ktable[ty,3] * l[0,0]) 
        a1.append(ktable[ty,3] * l[1,1])
        a2.append(ktable[ty,3] * l[2,2])
        a3.append(ktable[ty,3] * l[1,0])
        a4.append(ktable[ty,3] * l[2,0])
        a5.append(ktable[ty,3] * l[2,1])
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
    ty = -1
    k = []
    k1 = [] 
    a2 = [] 
    a3 = []
    a4 = []
    a5 = []
    for i,l in d.items():
        ty += 1 
        k.append(ktable[ty,4] * l[0,0]) 
        k1.append(ktable[ty,4] * l[1,1])
        a2.append(ktable[ty,4] * l[2,2])
        a3.append(ktable[ty,4] * l[1,0])
        a4.append(ktable[ty,4] * l[2,0])
        a5.append(ktable[ty,4] * l[2,1])
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
    
import tkinter 
window = tkinter.Tk()
window.title("Input box")

tkinter.Label(window, text = "Ef").grid(row = 0) 
e1 = tkinter.Entry(window)
e1.grid(row = 0, column = 1)

tkinter.Label(window, text = "Nuf").grid(row = 1)
e2 = tkinter.Entry(window)
e2.grid(row = 1, column = 1) 

tkinter.Label(window, text = "Ff").grid(row = 2) 
e3 = tkinter.Entry(window)
e3.grid(row = 2, column = 1)

tkinter.Label(window, text = "Em").grid(row = 3)
e4 = tkinter.Entry(window)
e4.grid(row = 3, column = 1)

tkinter.Label(window, text = "Num").grid(row = 4) 
e5 = tkinter.Entry(window)
e5.grid(row = 4, column = 1)

v = tkinter.StringVar(window, value='-1')

tkinter.Label(window, text = "Gf").grid(row = 5)
e6 = tkinter.Entry(window,textvariable=v)
e6.grid(row = 5, column = 1)
tkinter.Label(window, text = "If unkown input -1").grid(row = 5, column = 2)

tkinter.Label(window, text = "Gm").grid(row = 6)
e7 = tkinter.Entry(window, textvariable=v)
e7.grid(row = 6, column = 1)
tkinter.Label(window, text = "If unkown input -1").grid(row = 6, column = 2)

tkinter.Label(window, text = "Rhof").grid(row = 7)
e8 = tkinter.Entry(window)
e8.grid(row = 7, column = 1)

tkinter.Label(window, text = "Rhom").grid(row = 8) 
e9 = tkinter.Entry(window)
e9.grid(row = 8, column = 1)

var1 = tkinter.IntVar()
var2 = tkinter.IntVar()
tkinter.Checkbutton(window, text = "Weight fraction",variable = var1).grid(row=2,column = 2)
tkinter.Checkbutton(window, text = "Volume fraction", variable = var2).grid(row=2,column = 3)

tkinter.Label(window, text = "Method").grid(row = 9)

var3 = tkinter.IntVar()
var4 = tkinter.IntVar()
tkinter.Checkbutton(window, text = "Halpin-Tsai",variable = var3).grid(row=9,column = 1)
tkinter.Checkbutton(window, text = "Rule of Mixtures", variable = var4).grid(row=9,column = 2)

tkinter.Label(window, text = "Thickness").grid(row = 10)
e10 = tkinter.Entry(window)
e10.grid(row = 10, column = 1)
tkinter.Label(window, text = "Input in bottom to top layer order in mm and with a comma detween values i.e. 0.125,0.25,0.125,0.25").grid(row = 10, column = 2)

tkinter.Label(window, text = "Stacking seqences").grid(row = 11)
e11 = tkinter.Entry(window)
e11.grid(row = 11, column = 1)
tkinter.Label(window, text = "Input in bottom to top layer order and with a comma detween values i.e. 90,0,90,0").grid(row =11, column = 2)

v = tkinter.StringVar(window, value='3')
tkinter.Label(window, text = "Amount of sig figs in answear").grid(row = 12)
e12 = tkinter.Entry(window, textvariable=v)
e12.grid(row = 12, column = 1)

def test():
    Ff = float(e3.get())
    Rhof = float(e8.get()) 
    Rhom = float(e9.get())
    Z = int(e8.get())
    if (var1.get() == 1) & (var2.get() == 0):
        fraction = "W"              
    elif (var1.get() == 0) & (var2.get() == 1):
        fraction = "V"
    else:
        tkinter.messagebox.showinfo("Error", "please pick either weight or volume fraciton")
    if (var3.get() == 1) & (var4.get() == 0):
        method = "H"               
    elif (var3.get() == 0) & (var4.get() == 1):
        method = "M"
    else:
        tkinter.messagebox.showinfo("Error", "please pick a method of calculation")    
    Ef = float(e1.get()) 
    Nuf = float(e2.get())   
    Em = float(e4.get()) 
    Num = float(e5.get()) 
    Gf = float(e6.get()) 
    Gm = float(e7.get())
    E6 = e11.get().split(",")
    Z = int(e12.get())
    E7 = e10.get().split(",")
    pthick = []
    stse = []
    for i in E7:
        stse.append(int(i))
    for i in E6:
        pthick.append(int(i))
    return(Function3 (Ef, Ff, Rhof, Nuf, Gf, Em, Rhom, Num, Gm, method, fraction, stse, pthick, Z))
  
B = tkinter.Button(text="Compute", command=test)
B.grid(row = 13, column = 1)
window.mainloop()