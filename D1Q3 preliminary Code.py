import matplotlib.pyplot as plt
from matplotlib import animation
import math
import numpy as np


sT = 1024 
dt = 1

weights = np.array([2 / 3, 1 / 6, 1 / 6])
velocities = np.array([0, -1, 1])
cs2 = 1 / 3
rt=1.2

pos = 64


u = np.zeros(pos)
rho = np.array([2. for i in range(0, int(pos/2))] + [1. for i in range(int(pos/2), pos)])

f0 = np.zeros(pos)
f1 = np.zeros(pos)
f2 = np.zeros(pos)


def feq(i,m,n):
    equal=weights[i]*m*(1+3*velocities[i]*n+9/2*velocities[i]**2*n**2-3/2*n*n)
    return equal


def collideAndStream(rho, u):
    newU = np.zeros(pos) #declare auxilary variable for u
    newRho = np.zeros(pos) #auxilary variable for rho
    g0=np.zeros(pos)
    g1=np.zeros(pos)
    g2=np.zeros(pos)
    h0=np.zeros(pos)
    h1=np.zeros(pos)
    h2=np.zeros(pos)

    for x in range(0, pos):
        g0[x]=feq(0, rho[x], u[x])
        g1[x]=feq(1, rho[x], u[x])
        g2[x]=feq(2, rho[x], u[x])
        h0[x]=g0[x]*1/rt+(1-1/rt)*f0[x]
        h1[x]=g1[x]*1/rt+(1-1/rt)*f1[x]
        h2[x]=g2[x]*1/rt+(1-1/rt)*f2[x]
        #calculation of equlibrium and intermediate values
    for x in range(0, pos):    
        if(x == 0):
            f0[x]=h0[x]
            f1[x]=h1[x+1]
            f2[x]=h1[x] # streaming with boundary
        elif(x == pos - 1):
            f0[x]=h0[x]
            f1[x]=h2[x]
            f2[x]=h2[x-1]# streaming with boundary
        else:
            f0[x]=h0[x]
            f1[x]=h1[x+1]
            f2[x]=h2[x-1] # plain streaming 
    for x in range(0,pos):
        newU[x]=(velocities[0]*f0[x]+velocities[1]*f1[x]+velocities[2]*f2[x])/(f0[x]+f1[x]+f2[x]) #fill auxilary u
        newRho[x]=f0[x]+f1[x]+f2[x] #fill auxilary rho
    newURho=np.stack((newRho, newU), axis=0)
    return newURho


uData=[]
rhoData=[]

t = 0
uprho=rho
upu=u
uprho1=rho #all of those are auxilary variables that will be used in the loop
print(u)


while t <= sT:
    uprho1=collideAndStream(uprho, upu)[0,:]
    upu=collideAndStream(uprho, upu)[1,:]
    uprho=uprho1
    rhoData.append(uprho)
    uData.append(upu)
    t += dt
    print(max(abs(uprho)))
    if(max(abs(uprho))>5):
        print("error")
        print(t)


def rhoplot(t):
    rhopoints=rhoData[t]
    line=plt.plot(rhopoints)
    return line

rhoplot(32)
rhoplot(64)
plt.show()
