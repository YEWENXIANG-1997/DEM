# import the division module for compatibility between Python 2 and Python 3
from __future__ import division
# import the appropriate ESyS-Particle modules:
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *

import os
import sys
import time
import numpy as np

class tensileClass(Runnable):
    def __init__(self,
                 LsmMpi=None,
                 speed=1e-5,
                 endTs=1000,
                 iniD=1e-3,
                 r_p=1e-3
                 ):
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.speed = speed
        self.endTs = endTs
        self.iniD = iniD
        self.r_p = r_p
        self.dt = self.sim.getTimeStepSize()
        self.dx = self.dt*self.speed
        self.p0 = self.sim.getParticlePosn(1)
        self.p1 = self.sim.getParticlePosn(2)
        self.Lini = self.p0.toList()[2] - self.p1.toList()[2]
        self.strain = 0.0

    def run(self):
        self.ts = self.sim.getTimeStep()
        print(self.p0.toList()[2])
        print(self.p1.toList()[2])
        print(self.iniD)
        print(self.p0.toList()[2] - self.p1.toList()[2])
        Lnow = self.p0.toList()[2] - self.p1.toList()[2]
        self.strain = (Lnow - self.Lini) / (self.Lini) * 100
        if self.strain > 1.0:
            self.sim.exit()

        print("________________________")
        print("ts            :{:>8}".format(self.ts))
        print("end ts        :{:>8}".format(self.endTs))
        print("dx       [m/s]:{:>8}".format(self.dx))
        print("strain rate[%]:{:>8}".format(self.strain))
        
        self.p0 = self.p0+Vec3(0,0,+self.dx)
        self.p1 = self.p1+Vec3(0,0,-self.dx)

        self.sim.moveParticleTo(1, self.p0)
        self.sim.moveParticleTo(2, self.p1)
