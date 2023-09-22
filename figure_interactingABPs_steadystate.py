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
matplotlib.use('TkAgg')

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

##################
### PARAMETERS ###
##################

Pe=30
alpha=0.2
F0=100
LX=100
LY=400
init=1
Nran=100
N=5000

xfig=-0.31*LX
ly=100

for arg in sys.argv[1:]:
	if "-Pe=" in arg:
		Pe=float(arg[4:])
	elif "-alpha=" in arg:
		alpha=float(arg[7:])
	elif "-F0=" in arg:
		F0=float(arg[4:])
	elif "-LX=" in arg:
		LX=float(arg[4:])
	elif "-LY=" in arg:
		LY=float(arg[4:])
	elif "-init=" in arg:
		init=float(arg[6:])
	elif "-Nran=" in arg:
		Nran=int(arg[6:])
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)

Jmax=0.1
Amax=0.04
cmap=plt.get_cmap('jet')

###########################
### IMPORTATION OF DATA ###
###########################
	
if not os.path.isfile('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init)):
	RHO=0
	PX=0
	PY=0
	JX=0
	JY=0
	for ran in range(Nran):
		print(ran)
		RHO+=np.loadtxt('data_ABP2d_int_density/ABP2d_int_rho_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d.txt'%(Pe,alpha,F0,LX,LY,init,ran))
		PX+=np.loadtxt('data_ABP2d_int_polarization/ABP2d_int_px_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d.txt'%(Pe,alpha,F0,LX,LY,init,ran))
		PY+=np.loadtxt('data_ABP2d_int_polarization/ABP2d_int_py_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d.txt'%(Pe,alpha,F0,LX,LY,init,ran))
		JX+=np.loadtxt('data_ABP2d_int_current/ABP2d_int_jx_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d.txt'%(Pe,alpha,F0,LX,LY,init,ran))
		JY+=np.loadtxt('data_ABP2d_int_current/ABP2d_int_jy_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d.txt'%(Pe,alpha,F0,LX,LY,init,ran))

	rho0=N/(LX*LY*np.mean(RHO))

	RHO*=rho0
	PX*=rho0
	PY*=rho0
	JX*=rho0
	JY*=rho0
	
	os.system('mkdir -p data_npz')
	np.savez('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init),RHO=RHO,PX=PX,PY=PY,JX=JX,JY=JY)
else:
	RHO=np.load('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init))['RHO']
	PX=np.load('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init))['PX']
	PY=np.load('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init))['PY']
	JX=np.load('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init))['JX']
	JY=np.load('data_npz/ABP2d_intAv_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.npz'%(Pe,alpha,F0,LX,LY,init))['JY']
	
NX=len(RHO[0,:])
NY=len(RHO[:,0])
x=np.linspace(0,LX,NX+1)
y=np.linspace(0,LY,NY+1)

##############
### FIGURE ###
##############

fig=plt.figure(figsize=(8,5.5))
gs=matplotlib.gridspec.GridSpec(2,2,width_ratios=[1,1],height_ratios=[1,1],left=0.11,right=0.94,bottom=0.06, top=0.97, hspace=0.28,wspace=0.45)

###############
### DENSITY ###
###############

ax=plt.subplot(gs[0,0])

plt.pcolormesh(x,y,RHO,rasterized=True,cmap=cmap,norm=colors.LogNorm(vmin=0.01*RHO.max(), vmax=RHO.max()))
cb=plt.colorbar()
cb.solids.set_rasterized(True)

#plt.axis('equal')
plt.xlim([0,LX])
plt.ylim([0,ly])
plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
plt.yticks([0,0.25*ly,0.5*ly,0.75*ly,ly])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,ly,'{\\bf (a)}',ha='center',va='center')

#########################
### CURRENT AMPLITUDE ###
#########################

ax=plt.subplot(gs[0,1])

plt.pcolormesh(x,y,np.sqrt(JX**2+JY**2),vmin=0,vmax=Jmax,rasterized=True,cmap=cmap)
cb=plt.colorbar(ticks=[0,0.2*Jmax,0.4*Jmax,0.6*Jmax,0.8*Jmax,Jmax])
cb.solids.set_rasterized(True)
cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

#plt.axis('equal')
plt.xlim([0,LX])
plt.ylim([0,ly])
plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
plt.yticks([0,0.25*ly,0.5*ly,0.75*ly,ly])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,ly,'{\\bf (b)}',ha='center',va='center')

######################
### CURL AMPLITUDE ###
######################

ax=plt.subplot(gs[1,0])

dx=LX/(NX-1)
dy=LY/(NY-1)

A=0.5*(JY[1:,1:]+JY[:-1,1:]-JY[1:,:-1]-JY[:-1,:-1])/dx - 0.5*(JX[1:,1:]+JX[1:,:-1]-JX[:-1,1:]-JX[:-1,:-1])/dy
XX=0.5*(x[1:]+x[:-1])
YY=0.5*(y[1:]+y[:-1])

plt.pcolormesh(XX,YY,A,vmin=-Amax,vmax=Amax,rasterized=True,cmap=cmap)
cb=plt.colorbar(ticks=[-Amax,-0.5*Amax,0,0.5*Amax,Amax])
cb.solids.set_rasterized(True)
cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

st=plt.streamplot(XX,YY,JX,JY, color='k', linewidth=0.5, density=2, arrowstyle='->', arrowsize=1.)

#plt.axis('equal')
plt.xlim([0,LX])
plt.ylim([0,ly])
plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
plt.yticks([0,0.25*ly,0.5*ly,0.75*ly,ly])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,ly,'{\\bf (c)}',ha='center',va='center')

####################
### POLARIZATION ###
####################

ax=plt.subplot(gs[1,1])

THETA=np.arctan2(PY,PX)
THETA[RHO==0]=np.nan
THETA=np.ma.masked_where(np.isnan(THETA),THETA)

plt.pcolormesh(x,y,THETA,vmin=-np.pi,vmax=np.pi,rasterized=True,cmap=plt.get_cmap('hsv'))
cb=plt.colorbar(ticks=[-np.pi,-0.5*np.pi,0,0.5*np.pi,np.pi])
cb.ax.set_yticklabels(['$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'])
cb.solids.set_rasterized(True)

st=plt.streamplot(XX,YY,PX,PY, color='k', linewidth=0.5, density=1.5, arrowstyle='->', arrowsize=1.)

#plt.axis('equal')
plt.xlim([0,LX])
plt.ylim([0,ly])
plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
plt.yticks([0,0.25*ly,0.5*ly,0.75*ly,ly])
ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))

plt.text(xfig,ly,'{\\bf (d)}',ha='center',va='center')

os.system('mkdir -p figures')
plt.savefig('figures/figure_ABP2d_int_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.png'%(Pe,alpha,F0,LX,LY,init),dpi=250)
plt.savefig('figures/figure_ABP2d_int_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d.pdf'%(Pe,alpha,F0,LX,LY,init),transparent=True)
plt.close()

print('OK - time=%d sec'%(time.time()-clock))
