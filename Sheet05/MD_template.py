import numpy as np
from numpy import empty, zeros

"""
=================================================================================================
                     PROGRAM STRUCTURE

==============             ======                =============
| Thermostat | ----------- | MD | -------------- | Potential |
==============             ======                =============
                             |
                          ========
                          | main |
                          ========

- The MD class contains the primary logic of the algorithm
- Thermostat and Potential class are separated from it to allow different
thermostats and potentials can be implemented by inheritance and used flexibly.
can be used
- main() calls MD with the parameters needed for the task parts
=================================================================================================
"""


# ================================ Potential-class ================================================

# Virtual class from which concrete potentials can be inherited
# (Here only Lennard-Jones necessary, but so you can quickly implement other potentials).
class Potential:
    def V(r2):
        None

    def F(r):
        None

class PotentialLJ(Potential):
    def V(r2):
        # TODO
        # gets squared distance, returns potential
        None
    
    def F(r):
        # TODO
        # gets position returns force
        None

# ------------------------------ End of Potential-class -------------------------------------------

# ================================ Thermostat class ===============================================

# Virtual class from which concrete thermostats can be inherited
class Thermostat:
    def __init__(self,T_target = None):
        self.T_target = T_target

    def rescale(self,v, T):
        # gets list of 2d Vectors and temperature
        None

# No thermostat
class NoThermostat(Thermostat):
    def rescale(self,v, T):
        ...
        # TODO
        # gets list of 2d Vectors and temperature


# Isokinetic thermostat for task d)
class IsokinThermostat(Thermostat):
    def rescale(self,v, T):
        ...
        # TODO
        # gets list of 2d Vectors and temperature
        # returns new velocities 


# ------------------------------ End of Thermostat class ------------------------------------------
#


# In python there are no strict access modifiers (like private and protected in C)
# that restrict access to class members.
# But there is a commonly used naming convention to indicate with a leading _
# that a member is for internal use only
class MD:
    def __init__(self, L, N, T, potential, thermostat, numBins):
        self.N = N
        self.L = L
        self.numBins = numBins
        self.potential = potential()
        self.thermostat = thermostat(T)

        self.t = 0

        self.r = None # TODO  # shape(r) = (N,2)
        self.v = None # TODO  # shape(v) = (N,2)

    # Integration without data acquisition for pure equilibration
    def equilibrate(self, dt, n):
        ...
        # TODO
        # no return

    # Integration with data acquisition
    def measure(self, dt, n):
        self._init_data(n)
        self._update_data(0)

        # TODO

        return self.data


    def _init_data(self, n):
        self.data = {
            "t": empty(n+1 , dtype="float32"),
            "T": empty(n+1 , dtype="float32"),
            "Ekin": empty(n+1 , dtype="float32"),
            "Epot": empty(n+1 , dtype="float32"),
            "vCOMx": empty(n+1 , dtype="float32"),
            "vCOMy": empty(n+1 , dtype="float32"),
            "g": zeros(self.numBins, dtype="float32")#,
#            "r": empty((n+1, *self.r.shape), dtype="float32") # not needed unless you feel the urge to animate the simulation which is not required and does not provide additional points
            }

    def _update_data(self, i):
        ...
        # TODO
        # saves current state in self.data

    def _calc_T(self, Ekin=None):
        ...
        # TODO
        # returns temperature of current state

    def _calc_Ekin(self):
        ...
        # TODO
        # returns kinetic energy

    def _calc_Epot(self, r2):
        ...
        # TODO
        # gets squared distance, returns potential energy

    def _calc_vCOM(self):
        ...
        # TODO
        # returns center of mass velocity 

    # if you are not too familiar with numpy and advanced indexing it is easier to define the next two functions
    # so that they get two parameters i and j and only calculate the distance between two particles
    def _calc_DistanceVectors(self):
        ...
        # TODO
        # returns the distance vectors between all particles
        # shape(return) = (N, N-1, 2)

    def _calc_SquareDistances(self):
        ...
        # TODO
        # returns the squared distance between all particles
        # shape(return) = (N, N-1)


    def _calc_Acc(self):
        ...
        # TODO
        # calcualtes accelaration and either stores it in self.a or returns it

    def _verlet_step(self, dt):
       ...
       # TODO
       # performs a single interaction with the verlet algorithm



partPerRow = None # TODO
N = None # TODO
L = None # TODO
numBins = None # TODO Number of bins for pair correlation function


# b) Equilibration test
do_b = False
if do_b:
    T = None # TODO
    dt = None # TODO
    steps = None # TODO
    md = MD(L, N, T, LJ, noThermo, numBins)
    data = md.measure(dt, steps)

    # TODO Visualisation


# c) Pair correlation function
do_c = False
if do_c:
    for T in [0.01, 1, 100]:
        dt = None # TODO
        equiSteps = None  # TODO
        steps = None # TODO

        md = MD(L, N, partPerRow, T, LJ, noThermo, numBins)
        md.equilibrate(dt, equiSteps)
        data = md.meausure(dt, steps)
        

    # TODO Visualisation

# d) Thermostat
# TODO
