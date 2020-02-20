"""
Created on Sun Nov 24 14:43:23 2019
@author: Kieran & Elliot 
"""

import numpy as np
import pylab 
import random 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from random import gauss
import math
import time
 

k=0.4
rmax= 1
n = 30 # number of beads
T = 7 # Total simulation time
dt =1 # time step
setup1 = np.zeros((n,4)) # empty data ### make 4
#print (setup1)

PP = 3
QQ = 2
segment = 2*math.pi/(n)

def x(tt):
    """
    Computes and returns the source term (right-hand side)
    of the Poisson equation.
    
    Parameters
    ----------
    x : numpy.ndarray
        The gridline locations in the x direction
        as a 1D array of floats.
    y : numpy.ndarray
        The gridline locations in the y direction
        as a 1D array of floats.
    Lx : float
        Domain length in the x direction.
    Ly : float
        Domain length in the y direction.
    
    Returns
    -------
    b : numpy.ndarray of floats
        The forcing function as a 2D array.
    """
    x=(2+math.cos(PP*tt))*math.cos(QQ*tt)
    return(x)
    
def y(tt):
    y=(2+math.cos(PP*tt))*math.sin(QQ*tt)
    return(y)
    
def z(tt):
    z=math.sin(PP*tt)
    return(z)
    
for i in range(0, n):
    tt = i*segment 
    setup1[i,0] = x(tt)
    setup1[i,1] = y(tt)# starting position of beads (chain along x axis)
    setup1[i,2] = z(tt)



####   
print(setup1)

epsilon = 0.0000001
# min to equate minimum and lenth of knot or something
min1 = 1

def LJ(a,b):
    #print(a,b,"a,b")
    r=b-a
    #print(r,"r")
    rr = np.dot(r,r)
    rr = math.sqrt(rr)
    #print(rr,"rr")
    LJforce = (r/rr)*48*epsilon*((min1**12/(rr**13))-((min1**6/(2*rr**7))))
    return(LJforce)
    
forcevalues = np.zeros((T + dt, n, 4))
forcevalues[0,0:,0:]=setup1 # Need fix 
    
setup2 = setup1 # for time step
for i in range(dt,T + dt,dt): # dt, T+dt,dt
    setup1 = setup2
    setup3 = setup1
    for j in range(0, n):
      a = gauss(0, 1)/100
      b = gauss(0, 1)/100
      c = gauss(0, 1)/100

      if j==n-1:
          q = 0
      else:
          q = j+1
         
      spring1 = -setup1[j,0:] + setup1[j-1,0:] 
      #print(spring1,"initial")
      k21 = np.dot(spring1,spring1)
      sqrt1 = math.sqrt( k21 )
      #print(sqrt1,"sqrt1")
      if sqrt1 == 1:    
        spring1 =np.zeros((1,4))
        #print(spring1,"if")
      else:
          spring1 = spring1*(sqrt1-1)/sqrt1
          #print(spring1,"else")
          
      spring2 = (spring1*k/(1-(k21)/rmax**2)/m)*dt**2
      #print(spring2,"spring2")
      
      spring21 = -setup1[j,0:] + setup1[q,0:]
      #print(spring21,"initial")
      k22 = np.dot(spring21,spring21)
      sqrt2 = math.sqrt( k22 )
      #print(sqrt2,"sqrt2")
      if sqrt2 == 1:    
        spring21 =np.zeros((1,4))
        #print(spring21,"if")
      else:
          spring21 = spring21*(sqrt2-1)/sqrt2
          #print(spring1,"else")
      
      LJtotal = np.zeros(4)
      for h in range(0, n):
          if h!=j:
            #print(np.zeros((1,3)),"zeros")
            #print(setup1[j,0:],"setup1[j,0:]")
            #print(setup1[h,0:],"setup1[h,0:])")
            #print(LJ(setup1[j,0:],setup1[h,0:]),"LJ")  
            LJtotal = LJtotal + LJ(setup1[j,0:],setup1[h,0:])  # for all n here 
      spring22 = (spring21*k/(1-(k22)/rmax**2)/m)*dt**2
      
      
      ###print(a,"kick")    
      ###print(spring2,"spring2 F")    
      ###print(spring22,"spring22 F")
      setup2[j,0:] = setup1[j,0:] + [a+0.15,b,c,0] +  spring22 + spring2 + LJtotal
      ###print(setup2,"FINAL")
    setup2[0:,3] = i
    #print(setup2,i,"SETTTUP END")
    forcevalues[i,0:,0:] = setup2[0:,0:]


"""
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(setup2[0:,0],setup2[0:,1],setup2[0:,2]) 
"""

for j in range(0, n-1):
    if abs(setup2[j,0] -setup2[j+1,0])>1.5:
       print(j,setup2[j,0],"large x error")
    if abs(setup2[j,1] -setup2[j+1,1])>1.5:
       print(j,setup2[j,1],"large y error")
    if abs(setup2[j,2] -setup2[j+1,2])>1.5:
       print(j,setup2[j,2],"large z error")
    dotter = setup2[j,0:] -setup2[j+1,0:]
    dotter2 = np.dot(dotter,dotter)
    dotter3 = math.sqrt( dotter2 )
    if abs(dotter3)>1.5:
       print(dotter3,setup2[j,0:], "check error")
print("start",forcevalues)


forcevalues2 = np.zeros((T + dt, n+1, 4))
forcevalues2[:, :n, :] = forcevalues
forcevalues2[:, n, :] = forcevalues[:, 0, :]

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim3d([-3, 3.0])
ax.set_xlabel('X')

ax.set_ylim3d([-3, 3.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-1, 1.0])
ax.set_zlabel('Z')
    
for i in range(dt,T + dt,dt):
    ##forcevalues2[i,0:,0]+=i/10 difusion hack term

    ax.plot3D(forcevalues2[i,0:,0],forcevalues2[i,0:,1],forcevalues2[i,0:,2])
     