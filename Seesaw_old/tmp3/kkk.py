#!/usr/bin/env python
from pylab import *
import commands
import schmidt_spread
import sys

def Ifn(x):
    '''Loop function'''
    IFN=(x/(1-x))*(1-(x*log(x)/(1-x)))
    return IFN

def Lambda(Mi,Meta):
    '''Lambda definition for each right neutrino'''
    deltaMsol=2.43E-21
    h32=0.3**2
    Lambda3=deltaMsol/(3*h32)
    M3=6000
    vev=246
    lambda5=(Lambda3*M3*(4*pi)**2)/(2.*Ifn((M3/Meta)**2)*vev**2)
    return ((2.*lambda5*vev**2)/(Mi*(4*pi)**2))*Ifn((Mi/Meta)**2)

def h1func(psi,M1,Meta):
    deltaMatm=7.59E-23
    Lambda1=Lambda(M1,Meta)
    #print (deltaMatm**2/(4*Lambda1**2)+4.*psi**2)
    if (deltaMatm**2/(4*Lambda1**2)+4.*psi**2)>=0:
        h1_plus=((deltaMatm/(2.*Lambda1)+sqrt((deltaMatm**2/(4*Lambda1**2)+4.*psi**2))))/2
        h1_minus=((deltaMatm/(2.*Lambda1)-sqrt((deltaMatm**2/(4*Lambda1**2)+4.*psi**2))))/2
        if h1_plus>=0:
            return sqrt(h1_plus)
        else:
            print 'Error 1: Negative solution'
            pass
    else:
        print 'Error 2: Imaginary solution'
        pass

if  __name__=='__main__':
    gap=array([2.5,5.0,7.5,9.5]) #N1/eta
    mn1min=2050
    mn1max=3050
    mn1pts=6000
    MN1=linspace(mn1min,mn1max,mn1pts)
    relic=array([zeros_like(MN1),zeros_like(MN1),zeros_like(MN1),zeros_like(MN1)])
    psi=linspace(1E-2,10.,1000)
    ForSaving=[]
    #h1=h2=0.1
    #for j in xrange(len(gap)):
    #    for i in xrange(len(MN1)):
    #        relic[j][i]=schmidt.OMEGA(gap[j],MN1[i])
    #        print relic[j][i]
    #        ForSaving.append([MN1[i],0.01,relic[j][i],gap[j]])
    #np.save('schmidt_n_restrictions1',np.array(ForSaving))
    #h1=0.8,h2=0.125
    h1=0; h2=0
    kk=0
    for j in xrange(len(gap)):
        for i in xrange(len(MN1)):
            for k in xrange(len(psi)):
                h1=h1func(psi[k],MN1[i],gap[j]*MN1[i])
                h2=psi[k]/h1
                if kk%1000==0: print h1,h2, h1*h2
                relic[j][i]=schmidt_spread.OMEGA_minimal_complex(gap[j],MN1[i],h1,h2)
                if kk%1000==0: print relic[j][i]
                ForSaving.append([MN1[i],psi[k],relic[j][i],gap[j]])
                np.save('schmidt_Neu_Res_Spread',np.array(ForSaving))
                kk=kk+1
                if kk==10:
                    sys.exit()

    X=np.array(ForSaving)
    Om=0.1109
    DOm=0.0056
    nDOm=3
    XX=X[np.logical_and(X[:,2]>Om-nDOm*DOm,X[:,2]<Om+nDOm*DOm)]
    fig=figure()
    ax1=fig.add_subplot(1,1,1)
    #plotting=['rH','bv','k^','go']
    ax1.plot(XX[:,0],XX[:,1],'b^')
    #for h in xrange(len(gap)):
        #ax1.plot(XX[XX[:,3]==gap[h],plotting[h]])
    #ax1.hlines(0.2,11,1000)
    #ax1.hlines(0.08,11,1000)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$M_{1}$',size='x-large')
    ax1.set_ylabel(r'$\xi$',size='x-large')
    fig.savefig('Spread_Sector_A')
    #show()
    print 'Finish'
