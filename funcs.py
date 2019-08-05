# -*- coding: utf-8 -*-
from math import cos, pi, sqrt, exp, log
from scipy import fft, ifft








def gauleg(LeftExtremum, RightExtremum, nNodes, Nodes, Weights):
    """Compute Legende polynomial nodes and weights for gaussian quadrature

    loosely copied from Numerical Recipes third edition"""
    eps = 1e-14
    M = int((nNodes+1)/2)
    Pari = 1
    if( nNodes %2 ): Pari = 0
    xm = .5*(RightExtremum+LeftExtremum)
    xl = .5*(RightExtremum-LeftExtremum)
    for i in range(0,M):
        z = cos(pi*(i+0.75)/(nNodes+0.5))
        extrema = 1.0
        while extrema > eps :
            p1 = 1.0
            p2 = 0.0
            for j in range(0,nNodes):
                p3 = p2
                p2 = p1
                p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1)
            pp = nNodes*(z*p1-p2)/(z*z-1.0)
            z1 = z
            z = z1-p1/pp
            extrema = abs(z-z1)
        Nodes.insert(i, xm-xl*z)
        if(i != M-1 or Pari):
            Nodes.insert(len(Nodes)-i, xm+xl*z)
        Weights.insert(i, 2.0*xl/((1.0-z*z)*pp*pp) )
        if(i != M-1 or Pari):
            Weights.insert(len(Weights)-i, Weights[i] )



def omega(Ra,Rb,rab):
    """Overlap volume of two spheres of radii Ra and Rb.

    from Bianchi, Kahl & Likos, Soft Matter 7, 8313 (2011), formula 18"""
    omega = 0
    if ( rab >= Ra+Rb ):
        pass
    elif ( rab <= abs(Ra-Rb) ):
        omega = 8.0*(min(Ra,Rb)**3)
    else:
        cfrza = (Ra*Ra-Rb*Rb)/(2*rab)
        omega = 2*( (2*Ra+cfrza+rab/2)*((Ra-cfrza-rab/2)**2) \
                + (2*Rb-cfrza+rab/2)*((Rb+cfrza-rab/2)**2) )
    return omega


def calcpotIPC(r,c1,c2,phi, param):
    """Compute CGDH potential at the input coordinates

    r is the CM distance, so it should be rBiBj.
    c1, c2 are the cosines of the two polar angles and phi the azimuth angle difference.
    Uses omega, same reference."""
    calcpot = 0
    if ( r < param.Range ):
        s1 = sqrt(1-c1**2)
        s2 = sqrt(1-c2**2)
        rBiS1j    = sqrt( (r+param.ecc*c2)**2 + (param.ecc*s2)**2 )
        rBiS2j    = sqrt( (r-param.ecc*c2)**2 + (param.ecc*s2)**2 )
        rBjS1i    = sqrt( (r-param.ecc*c1)**2 + (param.ecc*s1)**2 )
        rBjS2i    = sqrt( (r+param.ecc*c1)**2 + (param.ecc*s1)**2 )
        tonda1    = s1*s2*cos(phi)+c1*c2
        due_e     = 2*param.ecc
        due_e2    = param.ecc*due_e
        rS1iS1j = sqrt( due_e2*(1-tonda1)+r*(r+due_e*(c2-c1)) )
        rS1iS2j = sqrt( due_e2*(1+tonda1)+r*(r-due_e*(c2+c1)) )
        rS2iS1j = sqrt( due_e2*(1+tonda1)+r*(r+due_e*(c2+c1)) )
        rS2iS2j = sqrt( due_e2*(1-tonda1)+r*(r+due_e*(c1-c2)) )
        # computing the potential
        w_bb = omega(param.bRadius,param.bRadius,r)
        w_bs = omega(param.bRadius,param.sRadius,rBiS1j) \
             + omega(param.bRadius,param.sRadius,rBiS2j) \
             + omega(param.bRadius,param.sRadius,rBjS1i) \
             + omega(param.bRadius,param.sRadius,rBjS2i)
        w_ss = omega(param.sRadius,param.sRadius,rS1iS1j) \
             + omega(param.sRadius,param.sRadius,rS1iS2j) \
             + omega(param.sRadius,param.sRadius,rS2iS1j) \
             + omega(param.sRadius,param.sRadius,rS2iS2j)
        calcpot = (param.e_bb*w_bb + param.e_bs*w_bs + param.e_ss*w_ss)/param.e_min

    return calcpot


def compute_gVW(rho, Nr, dr):
    """Compute the Verlet Weis HS g(r)

    Compute the HS g(r) from the Verlet-Weis correction to the Wertheim fit,
    see the comments below. References:
    Verlet & Weis,Phis. Rev. A 5 issue 2, 939
    """
    # compute the effective quantities for the VW corrections
    eta = pi*rho/6
    etaEFF = eta*(1.-eta/16)
    drEFF = dr/(etaEFF/eta)**(1./3.)
    dk = 2.*pi/(2.*Nr*drEFF)
    # compute the Wertheim c(r)
    lambda1 = (1+2*etaEFF)**2/(1-etaEFF)**4
    lambda2 = -(1+.5*etaEFF)**2/(1-etaEFF)**4
    cWr = []
    NdI = int(1.0/dr)
    Nd  = int(1.0/drEFF)

    for i in range(0,Nr):
        r = i*drEFF
        # To compute the 3d fourier transform of it, it has to be multiplied it by r!
        if ( i <= Nd ):
            cWr.append( (-lambda1 -6*etaEFF*lambda2*r - .5*etaEFF*lambda1*r**3)*r )
        else:
            cWr.append(0.0)
    cWr.append(0)
    for i in range(1,Nr):
        cWr.append(-cWr[Nr-i])
    cWk = fft(cWr)
    for i in range(1,len(cWk)):
        # 4pi/k and a 1/2i sine transform correction
        cWk[i] =  -cWk[i].imag*2*pi*drEFF/(i*dk)
    # cWk now contains the fourier transform.
    # hWk can be obtained solving the Ornstein Zernike which is now algebric:
    hWk = [0.0]
    for i in range(1,len(cWk)):
        hWk.append( (cWk[i]/(1-rho*cWk[i]))*i*dk )
    hWr = fft(hWk)
    gW = [ 0.0 ]
    A = .75*(etaEFF**2)*(1-etaEFF*(0.7117-0.114*etaEFF))/(1-etaEFF)**4
    g0W = 1.-hWr[NdI].imag*dk/( (2*pi**2)*NdI*drEFF )
    mu = 24*A/(etaEFF*g0W)

    for i in range(1,Nr):
        if(i <= Nd):
            gW.append(0)
        else:
            # 1/2pi^2r and a 1/2i sine transform correction
            r = i*dr
            coreC = 1.*(A/r)*exp(-mu*(r-1.0))*cos(-mu*(r-1.0))
            gW.append( 1.-hWr[i].imag*dk/( (2*pi**2)*(i*drEFF) ) + coreC )
    ## Now gW contains the Wertheim g!
    return gW


def calcpotKF(r,c1,c2,phi, param):
    """Kern Frenkel potential for testing"""
    if( r > param.lmbda ):
        return 0
    elif( c1 >= param.c0 and -c2 >= param.c0):
        return -1.0
    else:
        return 0


def samplepotential(dr,nNodes, x,pot,pot2, param):
    """Angularly averages the potential

    Stores <V(r)> in pot and <V^2(r)> in pot2;
    it is done at Core(=1./dr),Core+dr,Core+2dr,.... up to Range included
    """
    # Compute extrema
    Core = int(1.0/dr)
    End = int(param.Range/dr) - Core + 1
    # Sample node and weights
    cWeights, cNodes, pWeights, pNodes = [], [], [], []
    gauleg(-1., 1., nNodes, cNodes, cWeights)
    gauleg(0, 2*pi, nNodes, pNodes, pWeights)
    # Integrate
    for r in range(0,End):
        x.append(1.0+r*dr)
        Vci, Vci2 = 0., 0.
        for i in range(0,nNodes):
            Vcj, Vcj2 = 0., 0,
            for j in range(0,nNodes):
                Vphi, Vphi2 = 0., 0.
                for k in range(0,nNodes):
                    Vloc = calcpotIPC(x[r],cNodes[i],cNodes[j],pNodes[k], param)
                    Vphi  += Vloc*pWeights[k]
                    Vphi2 += Vloc*Vloc*pWeights[k]
                #return None
                Vcj  += Vphi*cWeights[j]
                Vcj2 += Vphi2*cWeights[j]
                #
            Vci  += Vcj*cWeights[i]
            Vci2 += Vcj2*cWeights[i]
            #
        pot.append(Vci/(8.*pi))
        pot2.append(Vci2/(8.*pi))
        #
    #fine




def comp_betaF(rho,kT, Nr,dr,nNodes,phi,phisq,x):
    """Computes betaFxc"""
    eta = pi*rho/6.
    beta = 1./kT
    kTdrhodp = ( (1-eta)**4 )/( 1+eta*(4+eta*(4+eta*(-4+eta))) )
    Core = int(1.0/dr)

    g = compute_gVW(rho, Nr, dr)
    #cacca = open('gVW','w')
    #for i in range(0,len(g)):
        #cacca.write( str(i*dr).rjust(5)+str(g[i].real).rjust(20)+'\n' )

    last = len(x)-1
    betaFxc = .5*x[0]*x[0]*g[Core]*( 2.*phi[0] - beta*kTdrhodp*phisq[0] )
    U = .5*x[0]*x[0]*g[Core]* phi[0]
    for i in range(1,last):
        betaFxc += x[i]*x[i]*g[i+Core]*( 2.*phi[i] - beta*kTdrhodp*phisq[i] )
        U       += x[i]*x[i]*g[i+Core]*phi[i]
    betaFxc += .5*x[last]*x[last]*g[last+Core]*( 2.*phi[last] - beta*kTdrhodp*phisq[last] )
    U       += .5*x[last]*x[last]*g[last+Core]*phi[last]
    betaFxc *= pi*rho*beta*dr #6.*eta*beta*dr
    U       *= 2.*pi*rho*dr
    
    try:
        betaFxc += log(rho) - 1. + eta*(4-3*eta)/(1-eta)**2
    except ValueError:
        print('Negative rho!')

    return betaFxc, U


def compute_state(rho,kT, outs, Nr,dr,nNodes,phi,phisq,x):
    Fleft , Uleft  = comp_betaF(rho*0.99,kT, Nr,dr,nNodes,phi,phisq,x)
    Fright, Uright = comp_betaF(rho*1.01,kT, Nr,dr,nNodes,phi,phisq,x)

    outs.P = rho*kT*(Fright - Fleft)/0.02
    outs.mu = kT*(1.01*Fright - .99*Fleft)/0.02
    outs.betaF = (Fright+Fleft)*.5
    outs.U_p = (Uright+Uleft)*.5

    Fleft , Uleft  = comp_betaF(rho,kT*0.99, Nr,dr,nNodes,phi,phisq,x)
    Fright, Uright = comp_betaF(rho,kT*1.01, Nr,dr,nNodes,phi,phisq,x)
    outs.U_p += (Uright+Uleft)*.5
    outs.U_t = -kT*kT*(Fright-Fleft)/.02

#   Check: since F = -pV + muN and also dF = mu dN - pdV,
#   d(rhoF/N)/drho = F/N + p/rho
    print(rho,outs.mu)
    if( abs(1 - (kT*outs.betaF+outs.P/rho)/outs.mu) > 1e-9):
        print("ERROR: Thermodynamic inconsistency!!!")
        print(T,rho)
        exit()



# Old version
#def compute_gVW(rho, Nr, drI):
    #"""Compute the Verlet Weis HS g(r)

    #Compute the HS g(r) from the Verlet-Weis correction to the Wertheim fit,
    #see the comments below. References:
    #Verlet & Weis,Phis. Rev. A 5 issue 2, 939
    #"""
    #etaI = pi*rho/6
    #eta = etaI*(1.-etaI/16)
    #dr = drI/(eta/etaI)**(1./3.)
    #dk = 2.*pi/(2.*Nr*dr)
    ## compute the Wertheim c(r)
    #lambda1 = (1+2*eta)**2/(1-eta)**4
    #lambda2 = -(1+.5*eta)**2/(1-eta)**4
    #cWr = []
    #NdI = int(1.0/drI)
    #Nd  = int(1.0/dr)
    #print Nd
    #print "etaI  = ", etaI, "eta  = ", eta, "dr = ", dr, "dk = ", dk
    #S = (1-eta)**4/( 1+eta*(4+eta*(4+eta*(-4+eta))) )#(1-8*eta**3+5*eta**4)  # from CS; from PY :(1-eta)**4/(1-4*eta**3)
    #print "S(0) = ", S
    #print "h(0) = ", (S-1.)/rho
    #print "c(0) = ", (1.-1./S)/rho

    #for i in range(0,Nr):
        #r = i*dr
        ## To compute the 3d fourier transform of it, it has to be multiplied it by r!
        #if ( i <= Nd ):
            #cWr.append( (-lambda1 -6*eta*lambda2*r - .5*eta*lambda1*r**3) )
            #cWr[i] *= r
        #else:
            #cWr.append(0.0)
    #cWr.append(0)
    #for i in range(1,Nr):
        #cWr.append(-cWr[Nr-i])
## fino innoi andasa
    #cWk = fft(cWr)
    #for i in range(1,len(cWk)):
        ## 4pi/k and a 1/2i sine transform correction
        #cWk[i] =  -cWk[i].imag*2*pi*dr/(i*dk)
    ## cWk now contains the fourier transform.
    ## hWk can be obtained solving the Ornstein Zernike which is now algebric:
    #hWk = [0.0]
    #for i in range(1,len(cWk)):
        #hWk.append( (cWk[i]/(1-rho*cWk[i])) )
        ##it is NOT already multiplied for k!
        #hWk[i] *= i*dk
    #hWr = fft(hWk)
    #gW = [ 0.0 ]
    #A = .75*(eta**2)*(1-eta*(0.7117-0.114*eta))/(1-eta)**4
    #g0W = 1.-hWr[NdI].imag*dk/( (2*pi**2) )
    #print "g0W = ",g0W
    #mu = 24*A/(eta*g0W)
    #print "mu = ", mu, "A = ",A

    #for i in range(1,Nr):
        #if(i <= Nd):
            #gW.append(0)
        #else:
            ## 1/2pi^2r and a 1/2i sine transform correction
            #r = i*drI
            #coreC = 1.*(A/r)*exp(-mu*(r-1.0))*cos(-mu*(r-1.0))
            #gW.append( 1.-hWr[i].imag*dk/( (2*pi**2)*(i*dr) ) + coreC )
    ### Now gW contains the Wertheim g!
    #return gW





#def compute_gW(rho, Nr, dr):
    #"""Compute the Wertheim HS g(r)
    #No longer in use

    #Compute the HS g(r) from the Verlet-Weis correction to the Wertheim fit,
    #see the comments below. References:
    #Verlet & Weis,Phis. Rev. A 5 issue 2, 939
    #"""
    ## compute the Wertheim c(r)
    #eta = pi*rho/6
    #dk = 2.*pi/(2.*Nr*dr) # ?
    #lambda1 = (1+2*eta)**2/(1-eta)**4
    #lambda2 = -(1+.5*eta)**2/(1-eta)**4
    #cWr = []
    #Nd = int(1.0/dr)
    ##print "eta  = ", eta, "dr = ", dr, "dk = ", dk
    ##S = (1-eta)**4/( 1+eta*(4+eta*(4+eta*(-4+eta))) )# from PY :(1-eta)**4/(1-4*eta**3)
    ##print "S(0) = ", S
    ##print "h(0) = ", (S-1.)/rho
    ##print "c(0) = ", (1.-1./S)/rho

    #for i in range(0,Nr):
        #r = i*dr
        ## To compute the 3d fourier transform of it, it has to be multiplied it by r!
        #if ( i <= Nd ):
            #cWr.append( (-lambda1 -6*eta*lambda2*r - .5*eta*lambda1*r**3) )
            #cWr[i] *= r
        #else:
            #cWr.append(0.0)
    #cWr.append(0)
    #for i in range(1,Nr):
        #cWr.append(-cWr[Nr-i])
## fino innoi andasa
    #cWk = fft(cWr)
    #for i in range(1,len(cWk)):
        ## 4pi/k and a 1/2i sine transform correction
        #cWk[i] =  -cWk[i].imag*2*pi*dr/(i*dk)
    ## cWk now contains the fourier transform.
    ## hWk can be obtained solving the Ornstein Zernike which is now algebric:
    #hWk = [0.0]
    #for i in range(1,len(cWk)):
        #hWk.append( (cWk[i]/(1-rho*cWk[i])) )
        ##it is NOT already multiplied for k!
        #hWk[i] *= i*dk
    #hWr = fft(hWk)
    #gW = [ 0.0 ]
    #for i in range(1,Nr):
        #if(i < Nd):
            #gW.append(0)
        #else:
            ## 1/2pi^2r and a 1/2i sine transform correction
            #gW.append( 1-(.5/pi**2)*hWr[i].imag*dk/( i*dr) )
    ### Now gW contains the Wertheim g!
    #return gW
