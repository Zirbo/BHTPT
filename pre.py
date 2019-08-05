#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# calcola la pressione SW dai dati di El Mendoub

from math import pi


rho, T, lmbda = 0., 0., 1.2
gc, gp, gm = 0,0,0

k=4

if( k==1):
    gc, gp, gm = 3.45, 0.75, 3.12
    T = 0.7; rho = 0.5

if( k==2):
    gc, gp, gm = 2.8,.89,2.4
    T = 1.0; rho = 0.5

if( k==3):
    gc, gp, gm = 2.67, 0.97, 2.63
    T = 1.0; rho = 0.1

if( k==4):
    gc, gp, gm = 3.17, 0.87, 2.37
    T = 1.0; rho = 0.7

print( 'SW = ',1+(2*pi*rho/3) * (gc + (gp-gm)*lmbda**3))