#!/usr/bin/env python

##
## Generating a polar plot of the electric potential 
##  calculated by the SMRAI 2.1 emulator
##
##  by S. Nakano (Apr. 2025)
##

import numpy as np
import matplotlib.pylab as plt
from multiprocessing import Pool

imgtype = 'eps'
#imgtype = 'png'

nprocs = 4 ## Number of the processes
#nprocs = 24 ## Number of the processes

nlon = 72
nlat = 20

theta=np.arange(0.0, 361.0, 5.0)*np.pi/180.0

r=np.arange(50.0, 90.0, 2.0)
r=90-r

tt, rr = np.meshgrid(theta,r) 

cb_min, cb_max,cb_div  = -60,60,24
interval_of_cf = np.linspace(cb_min, cb_max, cb_div+1)
interval_of_cf2 = np.linspace(0, 30, 16)


def hourlymap(hour):
  icount = 1

  fig = plt.figure(figsize=(4, 18))
#  fig = plt.figure(figsize=(10, 18))
  plt.subplots_adjust(wspace=0.1)
  plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95)

  for imin in range(0,60,10):
    pcmap = np.zeros((nlat, nlon+1))
    psigmap = np.zeros((nlat, nlon+1))

    adir = 'emltr/'

    pc = np.loadtxt(adir + "potest{0:0>2}{1:0>2}.dat".format(hour,imin))
    pcmap[:,0:nlon] = pc.reshape(nlat,nlon)
    pcmap[:,nlon] = pcmap[:,0]

    ax = plt.subplot(6,1,icount,polar=True)

    ax.set_title( "{0:0>2}:{1:0>2}".format(hour,imin), x=0.8,y=1 )
    ctf = plt.contourf(tt, rr, pcmap, interval_of_cf, cmap='bwr' )
    ax.set_theta_zero_location("S")
    ax.axes.set_rgrids((np.arange(4)+1)*10.0, ['80','70','60','50'])
    ax.axes.set_thetagrids(np.arange(4)*90.0, ['0','6','12','18'])

    cbar=plt.colorbar(shrink=0.6, pad=0.15, norm=ctf.norm, cmap=ctf.cmap)
    cbar.set_label('(kV)', fontsize=8, labelpad=-10, y=1.15, rotation=0)

    icount += 1

  plt.savefig("pot{0:0>2}.".format(hour,imin) + imgtype)
  plt.close()

  return


hrarr = np.arange(0,23)
nimg = len(hrarr)

if __name__ == '__main__':
  p = Pool(nprocs)
  p.map(hourlymap, range(0,nimg))
