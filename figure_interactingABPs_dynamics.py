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
ran=0

DT=2.
teq=0
tmax=4000

ly=100

multi=True
NCPU=6
movie=False

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
	elif "-ran=" in arg:
		ran=float(arg[5:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie" in arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU>1:
	multi=True
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)
	
#############################
### CREATION OF SNAPSHOTS ###
#############################
		
def Snapshot(i):
	if not os.path.isfile('snapshots/snapshot_ABP2d_int_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d_t=%.8g.png'%(Pe,alpha,F0,LX,LY,init,ran,t)):
		t=teq+i*DT
		data=np.loadtxt('data_ABP2d_int_position/ABP2d_int_position_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d_t=%.8g.txt'%(Pe,alpha,F0,LX,LY,init,ran,t))

		XX=data[:,1]
		YY=data[:,2]
		RR=data[:,3]
		tt=data[:,4]

		fig=plt.figure(figsize=(6,6))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.09,right=0.97, bottom=0.06,top=0.94,hspace=0.3,wspace=0.3)

		ax=plt.subplot(gs[0,0])		
		for k,(x,y) in enumerate(zip(XX,YY)):
			if y<ly+0.6:
				circle=plt.Circle((x,y), RR[k], color='None', ec='r',lw=0.5)
				ax.add_patch(circle)
				plt.arrow(x,y,1.5*RR[k]*np.cos(tt[k]),1.5*RR[k]*np.sin(tt[k]), head_width=0.3, head_length=0.3,fc='k', ec='k',lw=0.5,overhang=0.3,zorder=10)
			
		for axis in ['top','bottom','left','right']:
			ax.spines[axis].set_linewidth(0.5)
		ax.tick_params(width=0.5)

		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,ly])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*ly,0.5*ly,0.75*ly,ly])
		ax.xaxis.set_major_formatter(FormatStrFormatter('%.8g'))
		ax.yaxis.set_major_formatter(FormatStrFormatter('%.8g'))
		
		plt.text(0,1.03*ly,'$t=%d$'%(t),ha='left',va='center')
		plt.text(LX,1.03*ly,'${\\rm Pe}_s=%.8g$, $\\alpha=%.8g$, $F_0=%.8g$'%(Pe,alpha,F0),ha='right',va='center')
		
		plt.savefig('snapshots/snapshot_ABP2d_int_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d_%d.png'%(Pe,alpha,F0,LX,LY,init,ran,i),dpi=180)
		plt.close()
		
		print("t=%d -tcpu=%d sec"%(t,time.time()-clock))

ARG=[]
t=teq
i=0
while t<=tmax:
	ARG.append(i)
	t+=DT
	i+=1
	
os.system('mkdir -p snapshots')
print('%d snapshots'%len(ARG))
	
if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG)
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)
		
#############################
### CREATION OF THE MOVIE ###
#############################

if movie:
	os.system('mkdir -p movies')
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/snapshot_ABP2d_int_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -crf 30 -s 1080x1080 movies/movie_ABP2d_int_Pe=%.8g_alpha=%.8g_F0=%.8g_LX=%.8g_LY=%.8g_init=%d_ran=%d.mp4'%(Pe,alpha,F0,LX,LY,init,ran,Pe,alpha,F0,LX,LY,init,ran))

print('OK - time=%d sec'%(time.time()-clock))
