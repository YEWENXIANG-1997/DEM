# https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
from __future__ import division
from esys.lsm import *
from esys.lsm.util import *
import sys
import time
import random
import math
import numpy as np


class packingRunnable(Runnable):
    def __init__(self,
                LsmMpi=None,
                porosity=0.4,
                R=1e-3,
                H=4e-3,
                threshold_overlap=None,
                incT_snap=100,
                iniFac=0.1,
                prad=1.0
                 ):
        """
        Subroutine to initialise the Runnable and store parameter values.
        """
        Runnable.__init__(self)
        self.sim = LsmMpi
        self.sim.setRadiusExpansion(1)
        self.dt = self.sim.getTimeStepSize()
        self.tporo = porosity
        self.poro = 1.0
        self.ts = 0
        self.R = R
        self.H = H
        self.V_all = np.pi * self.R**2 * self.H
        self.endTs = 10000000
        self.addTs = 10000
        self.bswitch = True
        self.incTs_snap = incT_snap
        self.bStop = False
        self.nParticles = 100
        self.iniFac = iniFac
        self.cumlFac = 1.0
        self.prad = prad*iniFac
        self.beta = 0.3
        self.gamma = 1.001
        self.sim.setRadiusExpansionParams(beta=self.beta,gamma=self.gamma)
        self.sim.setParticleRadiusInitFactor(self.iniFac)

        self.threshold_overlap_l = threshold_overlap[:,1].tolist()
        self.threshold_index_l = threshold_overlap[:,0].tolist()
        self.nThreshold = len(self.threshold_index_l)
        self.threshold_overlap_l.append(-1.0)
        self.threshold_index_l.append(self.threshold_index_l[-1]+1)
        self.nThreshold = self.nThreshold + 1
        self.bFirstThreshold = []
        for i in range(0, self.nThreshold):
            self.bFirstThreshold.append(True)


    def run(self):
        self.ts = self.sim.getTimeStep()
        self.nParticles = self.sim.getNumParticles()
        if self.bStop == True:
            print("the distance of all particles is below the threshold.")
            self.sim.exit()
            sys.exit("stop the packing.")  

        if self.ts%100 == 0:
            print("________________________")
            print("ts              :{:>8}".format(self.ts))
            print("num of particles:{:>8}".format(self.nParticles))
            print("target porosity :{:>8}".format(self.tporo))
            print("current porosity:{:>8}".format(self.poro))
            if self.poro > self.tporo:
                print("current status  :{:>8}".format("scaling"))
            else:
                print("current status  :{:>8}".format("easing"))

        if (1./self.cumlFac >= self.iniFac):
            # if self.poro > self.tporo:
            self.scaling()
        else:
            self.isWithinDistance()
            self.sim.setRadiusExpansion(0)
            if self.bswitch is True:
                self.setEndTs()
                self.bswitch = False

        if self.ts == self.endTs:
            self.sim.exit()
            sys.exit("end Ts.")  

    def isWithinDistance(self):
        minDist = 0.0
        maxDist = 0.0
        dist = self.sim.getMinMaxDistance(name="friction").toList()
        minDist = dist[0] 
        maxDist = dist[1]

        for i in range(0, self.nThreshold-1):
            # print("_______")
            # print(self.threshold_overlap_l[i+1], abs(minDist),  self.threshold_overlap_l[i])
            # print(i,self.bFirstThreshold[i])
            # print("_______")
            if (self.threshold_overlap_l[i+1] <= abs(minDist)) and (abs(minDist) < self.threshold_overlap_l[i]):
                if self.bFirstThreshold[i] is True:
                    print("minimum overlap is abs(minDist)", abs(minDist))
                    # add a CheckPointer to store simulation data:
                    self.sim.createCheckPointer(
                        RestartCheckPointPrms(
                            fileNamePrefix="snapshot_threshold_idx"+str(int(self.threshold_index_l[i])),
                            beginTimeStep=0,
                            endTimeStep=self.ts,
                            timeStepIncr=1
                        )
                    )
                    self.bFirstThreshold[i] = False
                    if i == self.nThreshold-2:
                        self.bStop = True

    def scaling(self):
        # V = self.sim.getTotalVolume()
        V = self.getTotalVolume()
        self.poro = (self.V_all - V)/self.V_all

    def setEndTs(self):
        self.entTs = self.ts + self.addTs

    def getTotalVolume(self):
        fac = 1.0 + self.beta / self.ts**self.gamma
        self.cumlFac = self.cumlFac * fac
        self.prad = self.prad * fac
        V = self.nParticles * 4.0 / 3.0 * np.pi * self.prad**3
        return V

