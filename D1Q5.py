import matplotlib.pyplot as plt
from matplotlib import animation
import math
import numpy as np


sT = 1000
dt = 1
pt = 571
weights = np.array([1 / 2, 1 / 6, 1 / 6,1/12,1/12])
vel = np.array([0, -1, 1,-2,2])
cs2 = 1  #standard D1Q5 setting
rt =1 # to compare with the analytical solution

pos= 2000
rho1=2
rho2=1

u = np.zeros(pos)
rho = np.zeros(pos)
rho[0:int(pos/2)]=rho1
rho[int(pos/2):pos]=rho2

f=np.zeros([5,pos])


def feq(i,m,n):
    equal=weights[i]*m*(1+1*vel[i]*n+1/2*vel[i]**2*n**2-1/2*n*n)
    return equal

def collideAndStream(rho, u):
    newU = np.zeros(pos) #declare auxilary variable for u
    newRho = np.zeros(pos) #auxilary variable for rho
    Mom = np.zeros([5,pos]) #auxilary variable for urho
    vuse = np.zeros([5,pos])
    for z in range(0,5):
        vuse[z,:]=vel[z]
    g=np.zeros([5,pos])
    h=np.zeros([5,pos])
    for z in range(0, 5):
        g[z,:]=feq(z,rho[:],u[:])

        h[z,:]=g[z,:]*(1/rt)+f[z,:]*(1-1/rt)

        #calculation of equlibrium and intermediate values
    f[0,:]=h[0,:]
    f[1,0:pos-1]=h[1,1:pos]
    f[1,pos-1]=h[2,pos-1]
    f[2,1:pos]=h[2,0:pos-1]
    f[2,0]=h[1,0]
    f[3,0:pos-2]=h[3,2:pos]
    f[3,pos-2]=h[4,pos-1]
    f[3,pos-1]=h[4,pos-2]
    f[4,2:pos]=h[4,0:pos-2]
    f[4,1]=h[3,0]
    f[4,0]=h[3,1]
#  streaming and collision
    Mom=f*vuse
    newRho=f.sum(axis=0) #fill auxilary rho
    newu=Mom.sum(axis=0)/newRho
    newURho=np.stack((newRho, newu), axis=0)
    return newURho
print(collideAndStream(rho, u))


uData=[]
rhoData=[]

t = 0
uprho=rho
upu=u
uprho1=rho
upu1=u #all of those are auxilary variables that will be used in the loop



while t <= sT:
    uprho1=collideAndStream(uprho, upu)[0,:]
    upu1=collideAndStream(uprho, upu)[1,:]
    uprho=uprho1
    upu=upu1
    rhoData.append(uprho)
    uData.append(upu)
    t += dt



def rhoplot(t):
    rhopoints=rhoData[t]
    line=plt.plot(rhopoints)
    return line
def uplot(t):
    upoints=uData[t]
    line=plt.plot(upoints)
    return line


rho3=math.sqrt(rho2*rho1)
def solverhoa(a,b):
    init=(a/b)**(1/2)
    for i in range(0,10):
        intercept=init**0.5-init**-0.5+math.log(init)-math.log(a/b)
        slope=(init**0.5+1)**2/init**1.5
        newinit=init-intercept/slope
        init=newinit
    return init
rhoa=solverhoa(rho1,rho2)
print(rhoa)
ua=cs2**0.5*(math.log(rho1/rhoa))
print(ua)
ub=cs2**0.5*((rhoa/rho2)**0.5)
print(ub)


def anauplot(k):
    T1=-cs2**0.5*k+int(pos/2)
    T2=-cs2**0.5*k+ua*k+int(pos/2)
    T3=ub*k+int(pos/2)
    anau=np.zeros(pos)
    for T in range(int(T1),int(T2)):
        anau[T]=ua*(T-int(T1))/(int(T2)-1-int(T1))
    anau[int(T2):int(T3)]=ua
    line=plt.plot(anau)
    return line
def anarhoplot(k):
    T1=-cs2**0.5*k+int(pos/2)
    T2=-cs2**0.5*k+ua*k+int(pos/2)
    T3=ub*k+int(pos/2)
    anarho=np.zeros(pos)
    anarho[0:int(T1)]=rho1
    anarho[int(T2):int(T3)]=rhoa
    anarho[int(T3):pos]=rho2
    for T in range(int(T1),int(T2)):
        anarho[T]=rho1*((rhoa/rho1)**((T-int(T1))/(int(T2)-1-int(T1))))
    line=plt.plot(anarho)
    return line


uplot(pt)
anauplot(pt)
plt.show()


rhoplot(pt)
anarhoplot(pt)
plt.show()
