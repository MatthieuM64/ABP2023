#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.integrate as Int
import scipy.optimize as opt
import os
import sys
import time
import multiprocessing
from matplotlib.ticker import FormatStrFormatter

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

##################
### PARAMETERS ###
##################

vp=2
vg=0.5
LY=20.
LX=20.

xfig=-0.75*LX
xfigzoom=-0.562*LX

for arg in sys.argv[1:]:
	if "-vp=" in arg:
		vp=float(arg[4:])
	elif "-vg=" in arg:
		vg=float(arg[4:])
	elif "-LX=" in arg:
		LX=float(arg[4:])
	elif "-LY=" in arg:
		LY=float(arg[4:])
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
cmap=plt.get_cmap('jet')

###########################
### IMPORTATION OF DATA ###
###########################

if not os.path.isfile('data_npz/ABP2d_ideal_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.npz'%(vp,vg,LX,LY)):
	RHO=np.loadtxt('data_ABP2d_ideal/ABP2d_ideal_rho_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.txt'%(vp,vg,LX,LY))
	PX=np.loadtxt('data_ABP2d_ideal/ABP2d_ideal_px_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.txt'%(vp,vg,LX,LY))
	PY=np.loadtxt('data_ABP2d_ideal/ABP2d_ideal_py_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.txt'%(vp,vg,LX,LY))
	
	rho0=np.mean(RHO)
	RHO/=rho0
	PX/=rho0
	PY/=rho0
	
	RHO=0.5*(RHO+RHO[:,::-1])
	PX=0.5*(PX-PX[:,::-1])
	PY=0.5*(PY+PY[:,::-1])
	
	os.system('mkdir -p data_npz')
	np.savez('data_npz/ABP2d_ideal_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.npz'%(vp,vg,LX,LY),RHO=RHO,PX=PX,PY=PY)
else:
	RHO=np.load('data_npz/ABP2d_ideal_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.npz'%(vp,vg,LX,LY))['RHO']
	PX=np.load('data_npz/ABP2d_ideal_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.npz'%(vp,vg,LX,LY))['PX']
	PY=np.load('data_npz/ABP2d_ideal_vp=%.8g_vg=%.8g_LX=%.8g_LY=%.8g.npz'%(vp,vg,LX,LY))['PY']

NX=len(RHO[0,:])
NY=len(RHO[:,0])
x=np.linspace(-LX*0.5,LX*0.5,NX)
y=np.linspace(0,LY,NY)

##################################
### CALCULATION OF THE CURRENT ###
##################################

dx=LX/(NX-1)
dy=LY/(NY-1)

JX=-(RHO[:,1:]-RHO[:,:-1])/dx + 0.5*vp*( PX[:,:-1] + PX[:,1:] )
JY=-(RHO[1:,:]-RHO[:-1,:])/dy + 0.5*vp*( PY[:-1,:] + PY[1:,:] )  - 0.5*vg*( RHO[:-1,:] + RHO[1:,:] )

JX=0.5*(JX[:-1,:]+JX[1:,:])
JY=0.5*(JY[:,:-1]+JY[:,1:])
X=0.5*(x[1:]+x[:-1])
Y=0.5*(y[1:]+y[:-1])

A=0.5*(JY[1:,1:]+JY[:-1,1:]-JY[1:,:-1]-JY[:-1,:-1])/dx - 0.5*(JX[1:,1:]+JX[1:,:-1]-JX[:-1,1:]-JX[:-1,:-1])/dy
XX=0.5*(X[1:]+X[:-1])
YY=0.5*(Y[1:]+Y[:-1])

##############
### FIGURE ###
##############

fig=plt.figure(figsize=(8,6))
gs=matplotlib.gridspec.GridSpec(2,2,width_ratios=[1,1],height_ratios=[1,1],left=0.095,right=0.94,bottom=0.06, top=0.97, hspace=0.4,wspace=0.4)

###############
### DENSITY ###
###############

ax=plt.subplot(gs[0,0])	

plt.pcolormesh(x,y,RHO,rasterized=True,cmap=cmap,norm=colors.LogNorm(vmin=0.1, vmax=10))	
cb=plt.colorbar()
cb.solids.set_rasterized(True)

plt.contour(x,y,RHO,[1],colors='k',linewidths=1,linestyles=['-'])

rhoy=RHO[:,int(NX/2)+1]
rhoyy=rhoy[rhoy>=1]
h0=y[len(rhoyy)-1]
plt.plot([-0.5*LX,0.5*LX],[h0,h0],'--r',lw=1)

#plt.axis('equal')
plt.xlim([-0.5*LX,0.5*LX])
plt.ylim([0,LY])
plt.xticks([-0.5*LX,-0.25*LX,0,0.25*LX,0.5*LX])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,LY,'{\\bf (a)}',ha='center',va='center')

####################
### POLARIZATION ###
####################

ax=plt.subplot(gs[0,1])

THETA=np.arctan2(PY,PX)
plt.pcolormesh(x,y,THETA,vmin=-np.pi,vmax=np.pi,rasterized=True,cmap=plt.get_cmap('hsv'))
cb=plt.colorbar(ticks=[-np.pi,-0.5*np.pi,0,0.5*np.pi,np.pi])
cb.ax.set_yticklabels(['$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'])
cb.solids.set_rasterized(True)

st=plt.streamplot(x,y,PX,PY, color='k', linewidth=0.5, density=1., arrowstyle='->', arrowsize=1.)

#plt.axis('equal')
plt.xlim([-0.5*LX,0.5*LX])
plt.ylim([0,LY])
plt.xticks([-0.5*LX,-0.25*LX,0,0.25*LX,0.5*LX])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,LY,'{\\bf (b)}',ha='center',va='center')

######################
### CURL AMPLITUDE ###
######################

ax=plt.subplot(gs[1,0])	

Amax=np.abs(A).max()
plt.pcolormesh(XX,YY,A,vmin=A.min(),vmax=A.max(),rasterized=True,cmap=cmap)	
cb=plt.colorbar(ticks=[-Amax,-0.5*Amax,0,0.5*Amax,Amax])
cb.ax.set_yticklabels(['$%.2g$'%(-Amax),'$%.2g$'%(-0.5*Amax),'$0$','$%.2g$'%(0.5*Amax),'$%.2g$'%(Amax)])
cb.solids.set_rasterized(True)

st=plt.streamplot(X,Y,JX,JY, color='k', linewidth=0.5, density=1.25, arrowstyle='->', arrowsize=1.)

#plt.axis('equal')
plt.xlim([-0.5*LX,0.5*LX])
plt.ylim([0,LY])
plt.xticks([-0.5*LX,-0.25*LX,0,0.25*LX,0.5*LX])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,LY,'{\\bf (c)}',ha='center',va='center')

########################
### VORTEX SCHEMATIC ###
########################

ax=plt.subplot(gs[1,1])

Amax=A[:250,:250].max()
Amin=A[:250,:250].min()
cmap=plt.get_cmap('bwr')
plt.pcolormesh(XX,YY,A,rasterized=True,cmap=cmap,norm=colors.TwoSlopeNorm(vmin=Amin,vcenter=0,vmax=Amax))	
cb=plt.colorbar(ticks=[Amin,0.5*Amin,0,0.5*Amax,Amax])
cb.ax.set_yticklabels(['$%.2g$'%(-Amax),'$%.2g$'%(-0.5*Amax),'$0$','$%.2g$'%(0.5*Amax),'$%.2g$'%(Amax)])
cb.solids.set_rasterized(True)

st=plt.streamplot(X,Y,JX,JY, color='k', linewidth=0.5, density=3., arrowstyle='->', arrowsize=1.)

#plt.axis('equal')
plt.xlim([-0.5*LX,-0.25*LX])
plt.ylim([0,0.5*LY])
plt.xticks([-0.5*LX,-0.45*LX,-0.4*LX,-0.35*LX,-0.3*LX,-0.25*LX])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(-9.25,5.5,'$3$',ha='center',va='center',color='r',bbox=dict(facecolor='w', alpha=0.8, edgecolor='None',boxstyle="circle,pad=0.1"),fontsize=fontsize)
plt.text(-8.25,0.75,'$1$',ha='center',va='center',color='r',bbox=dict(facecolor='w', alpha=0.8, edgecolor='None',boxstyle='circle,pad=0.1'),fontsize=fontsize)
plt.text(-7,3.5,'$2$',ha='center',va='center',color='b',bbox=dict(facecolor='w', alpha=0.8, edgecolor='None',boxstyle='circle,pad=0.1'),fontsize=fontsize)

plt.text(xfigzoom,0.5*LY,'{\\bf (d)}',ha='center',va='center')

plt.savefig('figure_ABP2d_ideal_FEM.pdf')
plt.savefig('figure_ABP2d_ideal_FEM.png',dpi=250)
plt.close()	

print('OK - time=%d sec'%(time.time()-clock))
