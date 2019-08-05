#! /usr/bin/env python3
# -*- coding: utf-8 -*-
from funcs import samplepotential, comp_betaF, compute_gVW, compute_state
from math import pi
from scipy.optimize import broyden1
#requires both scipy and numpy!!!

# final version, fully computes the phase diagram during the loop

# 30n
# rho T    0.8172442322341,  -25.20951227009,   228.9672928480,   0.32,   0.28,   0.9468762331
# 45n
# rho T 0.2457275669234  -3.116945803881   21.22982043666   0.22   0.38   0.1424945265
# ISMENE
# rho T    55.82222168158  -315.2713205326   1281.080975680   0.3   0.325   15.71009391543
# altri
# rho T 6.494914733343  -73.05252812912   487.1847101358   0.40   0.35   15.12477874235
#
# YuKal (max)
# delta .1
# rho T   2.86,   -74.6,   661.0,   0.3,   0.25,  0.6683
# delta .2
# rho T   .674,   -16.9,   144.6,   0.3,   0.30,  0.6683
# delta .3
# rho T   .283,   -6.85,   571.0,   0.3,   0.35,  0.6683
#
# YuKal (tot)
# delta .1
# rho T   7.2158,   -193.095,   1721.9,   0.3,   0.25,  1.7327
# delta .2
# rho T   0.84315,   -21.5265,   185.89,   0.3,   0.30,  0.85293
# delta .3
# rho T   0.23453,   -5.7583,   48.469,   0.3,   0.35,  0.56169


# model
class param:
    #IPC
    e_bb, e_bs, e_ss, ecc, sRadius, e_min = 2.86,   -74.6,   661.0,   0.3,   0.25,  0.6683
    bRadius = sRadius + ecc; Range = 2.*bRadius
    #KF
    #chi, lmbda = .5, 1.5
    #c0 = 1.-2.*chi; Range = lmbda
    #Every model
# parameters
dr, Nr, nNodes, drho = .01, 1024, 20, 0.01
# states
kT = .32
rho_start, rho_end = .01,.6

rhoset = [rho_start]; r = rho_start
while(r < rho_end):
    r += drho
    rhoset.append(r)
#print(rhoset)

x, phi, phisq = [], [], []
rhobraketing = 0.0
samplepotential(dr,nNodes, x,phi,phisq, param)
cacca = open('pot','w')
for i in range(0,len(x)):
    cacca.write(repr(x[i]).ljust(20)+'\t'+repr(phi[i]).ljust(20)+'\t'+repr(phisq[i]).ljust(20)+'\n')


print('Potenziale calcolato')

for kT in [.5, .32, .23, .18]:
    risult = open('TPTsM1_T'+str(int(kT*100)),'w')
    risult.write(str("#rho").ljust(20)+'\t'+str("p").ljust(20)+'\t'+str("mu").ljust(20)+'\t'+str("betaF").ljust(24)+str("Uxc_pert").ljust(24)+'\tUxc_term\n')
    for rho in rhoset:
        class A: P = 0.; mu = 0.; betaF = 0.; U_p = 0.; U_t = 0.
        compute_state(rho,kT, A, Nr,dr,nNodes,phi,phisq,x)
        risult.write(repr(rho).ljust(20)+'\t'+repr(A.P).ljust(20)+'\t'+repr(A.mu).ljust(20)+'\t'+repr(A.betaF).ljust(24)+repr(A.U_p).ljust(24)+'\t'+repr(A.U_t).ljust(24)+'\n')
    #    if(A.P <= .0): rhobraketing = rho + drho
    print('Isoterma calcolata')#; estremi = ', rho_start, rhobraketing)

#def delta(rho):
#    class state: P = 0.; mu = 0.; betaF = 0.
#    compute_state(rho[0],kT, state, Nr,dr,nNodes,phi,phisq,x)
#    Dp = state.P;    Dmu = state.mu
#    state.P, state.mu, state.betaF = 0., 0., 0.
#    compute_state(rho[1],kT, state, Nr,dr,nNodes,phi,phisq,x)
#    Dp -= state.P;    Dmu -= state.mu
#    return  [Dp, Dmu]

#if(rhobraketing != 0):
#    rhoguess = [rho_start, rhobraketing]
#    sol = broyden1(delta, rhoguess,verbose=True)
#    ddf = open('coexC','a')
#    ddf.write(str(sol[0]).ljust(20)+'\t'+str(kT).ljust(20)+'\t'+str(sol[1]).ljust(20)+'\n')
#else:
#    print('No solution')



#ALWAYS REMEMBER that rho d/drho = eta d/deta
#print(rho)
#gigia = compute_gVW(rho, Nr, dr)
#cacca = open('g','w')
#for i in range(0,len(gigia)):
    #cacca.write(repr(i*dr).ljust(20)+'\t'+repr(gigia[i]).ljust(20)+'\n')
