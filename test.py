import os
import shutil
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *


# instantiate a simulation object and
# initialise the neighbour search algorithm:
sim = LsmMpi (numWorkerProcesses = 1, mpiDimList = [1,1,1])
sim.initVerletModel (
   particleType = "NRotSphere",
   gridSpacing = 0.5,
   verletDist = 0.04
)

# specify the number of timesteps and the timestep increment:
sim.setNumTimeSteps (1000)
sim.setTimeStepSize (1.0e-4)

# the length of the domain 
sim.setSpatialDomain(BoundingBox(Vec3(-1,-1,-1),Vec3(1,1,1)))

# add the first sand lay to the domain:
geoHCP = HexagBlock (
   dimCount = [10,10,3],
   radius = 0.2
)
geoHCP.translate ( translation = Vec3(0,0,-0.7) )

sim.createParticles(geoHCP)

# add the hydrate lay to the domain:
geoHydrate = RandomBoxPacker(
   minRadius = 0.02,
   maxRadius = 0.05,
   cubicPackRadius = 0.12,
   maxInsertFails = 2000,
   bBox = BoundingBox(
      Vec3(-1,-1,-0.7),
      Vec3(1,1,0.7)
   ),
   circDimList = [False, False, False],
   tolerance = 0.001
)
geoHydrate.generate()
geoRandomBlock_particles = geoHydrate.getSimpleSphereCollection()
sim.createParticles(geoRandomBlock_particles)

# add the second sand lay to the domain:
geoHCP2 = HexagBlock (
   dimCount = [10,10,3],
   radius = 0.2
)
geoHCP2.translate ( translation = Vec3(0,0,0.7))
sim.createParticles(geoHCP2)

# set the tag of first sand as #1
plist = []
for p1 in geoHCP:
   p1.setTag(1)
   plist.append(p1.getId())
maxId = max(plist)

# set the tag of hydrate as #2
for p2 in geoRandomBlock_particles:
   p2.setTag(2)
   oldId_2 = p2.getId() + 1
   p2.setId(oldId_2 + maxId)
   plist.append(p2.getId())
maxId = max(plist)

# set the tag of second sand as #1
for p3 in geoHCP2:
   p3.setTag(1)
   oldId_3 = p3.getId() + 1
   p2.setId(oldId_3 + maxId)
   plist.append(p3.getId())

print "Initial number particles:",sim.getNumParticles()
for n in range(1000):
   sim.runTimeStep()
   if (n%100==0): print "number particles:",sim.getNumParticles()
print "Final number particles:",sim.getNumParticles()
sim.exit()