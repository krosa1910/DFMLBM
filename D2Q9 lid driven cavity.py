import matplotlib.pyplot as plt
from matplotlib import animation
import math
import numpy as np

sT=100000
dt=1
tt=1000
pt=int(sT*0.8/tt)

weights = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
uvel = np.array([0,1,0,-1,0,1,-1,-1,1])
vvel = np.array([0,0,1,0,-1,1,1,-1,-1])
cs2 = 1 / 3 #standard D2Q9 setting
rt=1 # to compare with the analytical solution


625
301
748
#5=1+2,6=3+2,7=3+4,8=1+4


xpos=500
ypos=500
#size of seting
speed=0.1
#speed of lid
u = np.zeros([xpos,ypos])
v = np.zeros([xpos,ypos])
rho = np.ones([xpos,ypos])
#macroscopic velocities 
f = np.zeros([9,xpos,ypos])
# create a three dimensional array with each layer being an f0

def feq(i,d,l,n):
    dot=uvel[i]*l+vvel[i]*n
    equal=weights[i]*d*(1+3*dot+9/2*(dot)**2-3/2*(l*l+n*n))
    return equal
#define equilibrium
fk=np.zeros(9)
for z in range(0,9): 
    fk[z]=feq(z,1,0,0)
    f[z,:,:]=fk[z]
#initialize

def collideAndStream(rho,u,v):
    newU = np.zeros([xpos,ypos]) 
    #declare auxilary variable for u
    newV = np.zeros([xpos,ypos])
    #declare auxilary variable for v
    newRho = np.zeros([xpos,ypos]) 
    #auxilary variable for rho
    h=np.zeros([9,xpos,ypos])
    for z in range(9):
        h[z,:,:]=feq(z,rho[:,:],u[:,:],v[:,:])*1/rt+f[z,:,:]*(1-1/rt)
        #calculation of equlibrium and intermediate values
    f[0,:,:]=h[0,:,:]
    f[1,:,1:ypos]=h[1,:,0:ypos-1]
    f[2,1:xpos,:]=h[2,0:xpos-1,:]
    f[3,:,0:ypos-1]=h[3,:,1:ypos]
    f[4,0:xpos-1,:]=h[4,1:xpos,:]
    f[5,1:xpos,1:ypos]=h[5,0:xpos-1,0:ypos-1]
    f[6,1:xpos,0:ypos-1]=h[6,0:xpos-1,1:ypos]
    f[7,0:xpos-1,0:ypos-1]=h[7,1:xpos,1:ypos]
    f[8,0:xpos-1,1:ypos]=h[8,1:xpos,0:ypos-1]
    #regular streaming, boundaries will be overwritten
    f[1,:,0]=h[3,:,0]
    f[5,:,0]=h[7,:,0]
    f[8,:,0]=h[6,:,0]
    #bounce back at left boundary
    f[3,:,ypos-1]=h[1,:,ypos-1]
    f[7,:,ypos-1]=h[5,:,ypos-1]
    f[6,:,ypos-1]=h[8,:,ypos-1]
    #bounce back at right boundary
    f[2,0,:]=h[4,0,:]
    f[5,0,:]=h[7,0,:]
    f[6,0,:]=h[8,0,:]
    #bounce back at upper moving boundary
    f[4,xpos-1,:]=h[2,xpos-1,:]-6*weights[2]*rho[xpos-1,:]*(uvel[2]*speed)
    f[7,xpos-1,:]=h[5,xpos-1,:]-6*weights[5]*rho[xpos-1,:]*(uvel[5]*speed)
    f[8,xpos-1,:]=h[6,xpos-1,:]-6*weights[6]*rho[xpos-1,:]*(uvel[6]*speed)
    #bounce back at lower boundary


    newRho=f.sum(axis=0)
    momU=f*np.reshape(uvel,[9,1,1])
    momV=f*np.reshape(vvel,[9,1,1])

    newU=momU.sum(axis=0)/newRho
    newV=momV.sum(axis=0)/newRho
    newUVRho=np.stack((newRho, newU,newV), axis=0)
    return newUVRho

uData=[]
vData=[]
rhoData=[]

t = 0
uprho=rho
upu=u
upv=v
q=np.zeros([3,xpos,ypos])
#all of those are auxilary variables that will be used in the loop

while t <= sT:
    q=collideAndStream(uprho,upu,upv)
    uprho=q[0,:,:]
    upu=q[1,:,:]
    upv=q[2,:,:]
    if(t%tt==0):
        rhoData.append(uprho)
        uData.append(upu)
        vData.append(upv)
        print(t)
    t += dt

x = np.arange(0,xpos,1)
y = np.arange(0,ypos,1)
 
X, Y = np.meshgrid(x, y)
u=uData[100]
v=vData[100]
 
# creating plot
fig, ax = plt.subplots(figsize =(10, 10))
ax.quiver(X, Y, u, v)
 
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
ax.axis([-0.5, xpos+0.5, -0.5, ypos+0.5])
ax.set_aspect('equal')
 
# show plot
plt.show()

x = np.arange(0,xpos,1)
y = np.arange(0,ypos,1)
 
X, Y = np.meshgrid(x, y)
u=uData[100]
v=vData[100]
z=u
# creating plot
plt.subplot(1, 2, 1)
cs = plt.contourf(X, Y, z) 
  
cbar = plt.colorbar(cs) 
  
plt.title('100k') 
x = np.arange(0,xpos,1)
y = np.arange(0,ypos,1)
 
X, Y = np.meshgrid(x, y)
u=uData[100]
v=vData[100]
z=v
# creating plot
plt.subplot(1, 2, 2)
cs = plt.contourf(X, Y, z)

cbar = plt.colorbar(cs) 
plt.title('100k') 
plt.show() 

x = np.arange(0,xpos,1)
y = np.arange(0,ypos,1)
 
X, Y = np.meshgrid(x, y)
u=uData[50]
v=vData[50]
z=u
# creating plot
plt.subplot(1, 2, 1)
cs = plt.contourf(X, Y, z) 
  
cbar = plt.colorbar(cs) 
  
plt.title('50k') 
x = np.arange(0,xpos,1)
y = np.arange(0,ypos,1)
 
X, Y = np.meshgrid(x, y)
u=uData[50]
v=vData[50]
z=v
# creating plot
plt.subplot(1, 2, 2)
cs = plt.contourf(X, Y, z) 
  
cbar = plt.colorbar(cs) 
  
plt.title('50k') 
plt.show() 

x = np.arange(0,xpos,1)
y = np.arange(0,ypos,1)
 
X, Y = np.meshgrid(x, y)
u=uData[50]
v=vData[50]
 
fig = plt.figure(figsize = (10, 10))
 
# Plotting stream plot
plt.streamplot(X, Y, u, v, density = 5)
 
# show plot
plt.title("T=50000,x=500,y=500,lid velocity=0.1") 
plt.show()





