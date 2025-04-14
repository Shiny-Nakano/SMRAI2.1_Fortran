#!/usr/bin/env python

## 
## Quasi-uniformly determining locations of 
## radial basis functions on a sphere 
## according to the method by Saff & Kuijlaars (1997).
## 
## Referece:
##  Saff & Kuijlaars, Math. Intelligencer, 11, 5-11, 1997. 
##
## Usage: 
##  > python3 genpoles_spiral.py
##
## Output: polepoints_sp.dat
##

import numpy as np
import matplotlib.pylab as plt

latbound = 50.0
N = 2000

NN = int(N * 2.0 / (1.0 - np.cos(np.pi*(90.0-latbound)/180.0)))

phi = 0.0
h = -1.0
theta = np.arccos(h)
xlambda = theta - np.pi / 2.0

fout = open("polepoints_sp.dat", 'w')

fout.write("{0:18.14f}{1:18.14f}\n".format(xlambda, phi))


for k in range(2,NN):
  h = -1.0 + 2.0*(k-1.0)/(NN-1)
  theta = np.arccos(h)
  xlambda = theta - np.pi / 2.0
  phi  = phi + 3.6 / ( np.sqrt(1.0*NN) * np.sqrt(1.0 - h*h) )

  while phi > 2.0*np.pi:
    phi -= 2.0*np.pi

  if 180.0*xlambda/np.pi < 50.0:
    break

  fout.write("{0:18.14f}{1:18.14f}\n".format(xlambda, phi))

fout.close()
