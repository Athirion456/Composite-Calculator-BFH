# -*- coding: utf-8 -*-
"""
at1u16
"""
def Function1(Ef, Ff, Rhof, Nuf, Gf, Em, Rhom, Num, Gm, method, fraction,Z):
    "method M for rule of mixture and H for halpin tsai,fraction W for weight fraction and V for volume fraction, Z sets amount of decimal points in answer "
#   Set the function outputs as global variable so they can be returned from outside the tkinter window.  
    global E1
    global E2
    global V12
    global G12
#    This section converts Ff to volume fractions, the method depends on the fraction input.
    if fraction == "W":
            Vf = ((Ff / Rhof) / ((Ff / Rhof) + ((1-Ff) / Rhom)))
    elif fraction == "V":
            Vf = Ff 
    else:
#        Returns an error message
        return(print("Incorrect fraction input"))
#    Calculates the matrix volume fractions.
    Vm = 1 - Vf 
#    This section calculates the laminate properties. 
#    Firstly, the method of calculations is defined using the method input.
    if method == "M":
        "Rule of mixturs"
        E1 = (Ef * Vf) + (Em * Vm) 
        E2 = (Ef * Em)/ ((Ef * Vm) + (Em * Vf))
        V12 = (Vf * Nuf) +  (Vm * Num) 
#        Legacy code
        #V21 = V12 * (E22 / E11)
        if (Gf * Gm) == 1:
#            If Gf and Gm unknown they are calculated using Ef, Em, Num and Nuf.
#            Then G12 is calculated using the relevant method  
            Gf = Ef / (2 * ( 1 + Nuf))
            Gm = Em / (2 * ( 1 + Num))
            G12 =(Gf * Gm) / ((Gf * Vm) + (Gm * Vf))
        else:
#            If Gf and Gm are known G12 is calculated directly using the relevant method
            G12 =(Gf * Gm) / ((Gf * Vm) + (Gm * Vf))
    elif method == "H":
        "Halpin Tsai"
        E1 = (Ef * Vf) + (Em * Vm) 
#        Halpin-Tsai parameter is calculated 
        n1 = ((Ef / Em) - 1)/( (Ef / Em) + 0.2)
        E2 = Em * ((1 + (0.2 * n1 * Vf))/(1 - (n1 * Vf)))
        V12 = (Vf * Nuf) +  (Vm * Num) 
#        Legacy code
        #V21 = V12 * (E22 / E11)
        if (Gf * Gm) == 1:
#            If Gf and Gm unknown they are calculated using Ef, Em, Num and Nuf.
#            Then G12 is calculated using the relevant method 
                 Gf = Ef / (2 * ( 1 + Nuf))
                 Gm = Em / (2 * ( 1 + Num))
#                 Halpin-Tsai parameter is calculated 
                 n2 = ((Gf / Gm) - 1)/( (Gf / Gm) + 1)
                 G12 =Gm * ((1 + (1 * n2 * Vf))/(1 - (n1 * Vf)))
        else: 
#           If Gf and Gm are known G12 is calculated directly using the relevant method    
#            Halpin-Tsai parameter is calculated 
            n2 = ((Gf / Gm) - 1)/( (Gf / Gm) + 1)
            G12 =Gm * ((1 + (1 * n2 * Vf))/(1 - (n1 * Vf)))
    else:
#        Return an error message
        return(print("Incorrect method input"))
#    This section formats the result and applies the decimal factor prior to printing.
    T1 = "E1 = {E1f} MPa".format(E1f = round(E1,Z))
    T2 = "E2 = {E2f} MPa".format(E2f = round(E2,Z))
    T3 = "V12 = {V12f}".format(V12f = round(V12,Z))
    T4 = "G12 = {G12f} MPa".format(G12f = round(G12,Z))
    print(T1)
    print(T2)
    print(T3)
    print(T4)
#    The results are returned so they can be used in subsequent calculations 
    return(E1,E2,V12,G12) 

#This section deals with the creation of a GUI.

import tkinter 
window = tkinter.Tk()
window.title("Input box")

#This section set up the input box and box labels for the tkinter window. 

tkinter.Label(window, text = "Young’s modulus of fibre material (Ef)").grid(row = 0) 
e1 = tkinter.Entry(window)
e1.grid(row = 0, column = 1)
tkinter.Label(window, text = "GN/m^2").grid(row = 0, column = 2,sticky= tkinter.W)

tkinter.Label(window, text = "Poisson’s ratio of fibre material (Nuf)").grid(row = 1)
e2 = tkinter.Entry(window)
e2.grid(row = 1, column = 1) 

tkinter.Label(window, text = "Fibre fraction (Ff)").grid(row = 2) 
e3 = tkinter.Entry(window)
e3.grid(row = 2, column = 1)

tkinter.Label(window, text = "Young’s modulus of matrix material (Em)").grid(row = 3)
e4 = tkinter.Entry(window)
e4.grid(row = 3, column = 1)
tkinter.Label(window, text = "GN/m^2").grid(row = 3, column = 2,sticky= tkinter.W)

tkinter.Label(window, text = "Poisson’s ratio of matrix material (Num)").grid(row = 4) 
e5 = tkinter.Entry(window)
e5.grid(row = 4, column = 1)

#Defines a default variable for Gm and Gf input box
v = tkinter.StringVar(window, value='-1')

tkinter.Label(window, text = "Shear modulus fibre material (Gf)").grid(row = 5)
e6 = tkinter.Entry(window,textvariable=v)
e6.grid(row = 5, column = 1)
tkinter.Label(window, text = "GN/m^2 (If unkown input -1)").grid(row = 5, column = 2)

tkinter.Label(window, text = "Shear modulus matrix material (Gm)").grid(row = 6)
e7 = tkinter.Entry(window, textvariable=v)
e7.grid(row = 6, column = 1)
tkinter.Label(window, text = "GN/m^2 (If unkown input -1)").grid(row = 6, column = 2)

tkinter.Label(window, text = "Density of fibre material (Rhof)").grid(row = 7)
e8 = tkinter.Entry(window)
e8.grid(row = 7, column = 1)
tkinter.Label(window, text = "kg/m^3").grid(row = 7, column = 2,sticky= tkinter.W)

tkinter.Label(window, text = "Density of matrix material (Rhom)").grid(row = 8) 
e9 = tkinter.Entry(window)
e9.grid(row = 8, column = 1)
tkinter.Label(window, text = "kg/m^3").grid(row = 8, column = 2,sticky= tkinter.W)

#Defines a default variable for the amount of decimal points in answer 
v = tkinter.StringVar(window, value='3')

tkinter.Label(window, text = "Amount of decimal points in answer").grid(row = 9) 
e10 = tkinter.Entry(window, textvariable=v)
e10.grid(row = 9, column = 1)

tkinter.Label(window, text = "Method").grid(row = 10)

#Set up tick box to define fraction and method

var1 = tkinter.IntVar()
var2 = tkinter.IntVar()
tkinter.Checkbutton(window, text = "Weight fraction",variable = var1).grid(row=2,column = 2)
tkinter.Checkbutton(window, text = "Volume fraction", variable = var2).grid(row=2,column = 3)

var3 = tkinter.IntVar()
var4 = tkinter.IntVar()
tkinter.Checkbutton(window, text = "Halpin-Tsai",variable = var3).grid(row=10,column = 1)
tkinter.Checkbutton(window, text = "Rule of Mixtures", variable = var4).grid(row=10,column = 2)

def test():
#    Using the tick box inuts to define fraction and method.
    if (var1.get() == 1) & (var2.get() == 0):
        fraction = "W"               
    elif (var1.get() == 0) & (var2.get() == 1):
        fraction = "V" 
    else:
#        Return an error message
        tkinter.messagebox.showinfo("Error", "please pick either weight or volume fraciton")
    if (var3.get() == 1) & (var4.get() == 0):
        method = "H"               
    elif (var3.get() == 0) & (var4.get() == 1):
        method = "M"
    else:
#        Return an error message
        tkinter.messagebox.showinfo("Error", "please pick a method of calculation")    
#    Extract the input from the input boxes.
    Ff = float(e3.get())
    Rhof = float(e8.get()) 
    Rhom = float(e9.get()) 
    Ef = float(e1.get()) 
    Nuf = float(e2.get())   
    Em = float(e4.get()) 
    Num = float(e5.get()) 
    Gf = float(e6.get()) 
    Gm = float(e7.get()) 
    Z = int(e10.get())
#    Runs the function 1
    return(Function1(Ef, Ff, Rhof, Nuf, Gf, Em, Rhom, Num, Gm, method, fraction, Z))

#Set up the button which runs function 1 from the tkinter window   
B = tkinter.Button(text="Compute E11, E22, V12 and G12", command=test)
B.grid(row = 11, column = 1)
window.mainloop()