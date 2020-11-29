#!/usr/bin/python


"""
Bethe Bloch equation and calculation of energy loss along trajectory

"""

from math import *

#Define some natural constants

mp = 1.672621637e-27        #kg
me = 9.1093821499999992e-31 #kg
qe = 1.602176487e-19        #C
na = 6.02214179e23          #mol^-1
eps0 = 8.854187817e-12      #F/m
c = 299792458               #m/s


#and now some problem specific ones
def get_joule_from_gev(gev):
    return gev * 1.6e-10
def get_gev_from_joule(joule):
    return joule / 1.6e-10

# projectile properties
initial_energy = get_joule_from_gev(10) #1 GeV in J
m0   = 1.88e-28                   #mass of projectile (muon) in kg
Z1   = 1                    #electric charge of projectile
# material properties
EB   = 23.6*qe                #75 eV in J, i.e., average ionization energy of LAr
Z = 18
A = 39.948 # g/mol
rho = 1396*1000 # g/m^3
Mu = 1.0 # constante de masse molaire g/mol
ne   = na * Z * rho / (A * Mu)  # electron density in 1/m^3
print "Electron density, ", ne, " (water is 3.34e29)"
length_traversed_per_gap = 1.618377e-3 # mm 
n_gap_per_cell = 4

Ekin = initial_energy
print "Initial energy: %f GeV"%get_gev_from_joule(Ekin)

def beta(v):
    return v/c

def gamma(v):               #Lorentz gamma
    return 1./sqrt(1-v*v/c/c)

def v_of_Ekin_m0(Ekin, m0): #invert kinetic energy, E_kin, for speed, v.
    b2 = 1.-1./(1.+Ekin/m0/c/c)**2
    return sqrt(b2)*c

def dEdx(Z1,Ekin,m0,EB,ne): #Bethe-Bloch equation
    v = v_of_Ekin_m0(Ekin, m0)
    b2 = beta(v)**2
    C = Z1**2*qe**4/4/pi/eps0**2/me
    ln_term = log(2.*me*v**2/EB)
    return C/v**2*ne*(ln_term  - log(1.-b2) - b2)


#Energy loss in first layer
print str(dEdx(1,Ekin,m0,EB,ne)/(qe*1.e9)) + ' in MeV/mm'

xMax = 1.618377e-3 # mm 
#initialize position, energy loss, and dx
x=0         #position in mm
dE = 0.     #energy loss
dx = 1.e-4  #0.1mm
#bbf = open('bb_in_water.dat','w')
while x < xMax :
    string = str(x) + ', ' + str(Ekin/(qe*1e6)) + ', ' + str(dE/(qe*1.e9)/dx) + '\n'
    #print x, Ekin/(qe*1e6), dE/(qe*1.e9)/dx
    #print string
    #bbf.write(string)
    dE = dEdx(Z1,Ekin,m0,EB,ne)*dx     #units J/m*dx
    x = x+dx
    Ekin = Ekin - dE

#bbf.close()

print "Final energy: %f GeV"%get_gev_from_joule(Ekin)
print "Energy lost per gap: %f MeV"%(1000*get_gev_from_joule(initial_energy - Ekin))
print "Energy lost per cell: %f MeV"%(1000*get_gev_from_joule(initial_energy - Ekin)*4)
