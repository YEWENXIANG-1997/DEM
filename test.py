import os
import shutil
import numpy as np
from esys.lsm import *
from esys.lsm.util import *
from esys.lsm.geometry import *
from tensileRunnable import tensileClass
from WallLoader import WallLoaderRunnable


def makedirectory(dirname):
    if os.path.exists(dirname):
        print(dirname, "exist. Cleared contents.")
        shutil.rmtree(dirname)
        os.mkdir(dirname)
    else:
        print(dirname, "doesn't exist. made")
        os.mkdir(dirname)


bBond = False

out_dat_dir = "out_data"
disp_dir = "displacement"
force_dir = "nforce"
position_dir = "position"
sigmad_dir = "sigmaD"

dirname = "./"+out_dat_dir
makedirectory(dirname)
# dirname = "./"+out_dat_dir + "/" + disp_dir
# makedirectory(dirname)
dirname = "./"+out_dat_dir + "/" + force_dir
makedirectory(dirname)
# dirname = "./"+out_dat_dir + "/" + position_dir
# makedirectory(dirname)
# dirname = "./"+out_dat_dir + "/" + sigmad_dir
# makedirectory(dirname)


con_d = np.genfromtxt(
    "conditions.csv", delimiter=" ", dtype="U", autostrip=True)
for c in con_d:
    c = c.split(",")
    print(c)
    if c[0] == "dT[s]":
        dT = float(c[1])
    if c[0] == "incT[s]":
        incT_data = float(c[1])
        incTs_data = int(incT_data/dT)
    if c[0] == "incT_snap[s]":
        incT_snap = float(c[1])
        incTs_snap = int(incT_snap/dT)
    if c[0] == "speed[m/s]":
        speed = float(c[1])
    if c[0] == "L[m]":
        L = float(c[1])
    if c[0] == "R[m]":
        R = float(c[1])
    if c[0] == "initialDistance[m]":
        iniD = float(c[1])
    if c[0] == "endstrain":
        t = iniD*float(c[1]) / (2.*speed)
        endTs = int(t / dT)
        # endTs = 100000
    if c[0] == "r_p[m]":
        r_p = float(c[1])
    if c[0] == "r_s[m]":
        r_s =float(c[1])
    if c[0] == "density[kg/m^3]":
        density = float(c[1])
    if c[0] == "maxTensile[Pa]":
        maxTensile = float(c[1])
    if c[0] == "Sh":
        Sh = float(c[1])
    if c[0] == "k0":
        k0 = float(c[1])
    if c[0] == "k1":
        k1 = float(c[1])
    if c[0] == "bondModulus":
        bondModulus = float(c[1])
    if c[0] == "bondPoissonRatio":
        bondPoissonRatio = float(c[1])
    if c[0] == "bondCohesion":
        bondCohesion = float(c[1])
    if c[0] == "bondTanAngle":
        bondTanAngle = np.tan(np.radians(float(c[1])))
    if c[0] == "beta1":
        beta1 = float(c[1])
    if c[0] == "beta2":
        beta2 = float(c[1])
    if c[0] == "normalK_sand[N/m]":
        normalK_sand = float(c[1])
    if c[0] == "kn/ks":
        alpha = float(c[1])
        shearK_sand = normalK_sand / alpha
    if c[0] == "dynamicMu":
        dynamicMu = float(c[1])
    if c[0] == "staticMu":
        staticMu = float(c[1])
    if c[0] == "viscosity_damp_nr":
        viscosity_damp_nr = float(c[1])


# instantiate a simulation object and
# initialise the neighbour search algorithm:
sim = LsmMpi (numWorkerProcesses = 1, mpiDimList = [1,1,1])
sim.initVerletModel (
   particleType = "RotSphere",
   gridSpacing = 2*r_p+0.2*r_p,
   verletDist = 0.2*r_p
)

# specify the number of timesteps and the timestep increment:
sim.setNumTimeSteps (2000)
sim.setTimeStepSize (1.0e-4)

# the length of the domain 
ep = 1e-5
sim.setSpatialDomain(BoundingBox(Vec3(-4e-3-ep,-4e-3-ep,-4.0e-3-ep),Vec3(4e-3+ep,4e-3+ep,4.0e-3+ep)))

# add the first sand lay to the domain:
# set the tag as #1
count = int(L/r_s)
# geoHCP = HexagBlock (
#    dimCount = [count,count,3],
#    radius = r_s
# )
# geoHCP.translate ( translation = Vec3(0,0,-4.0e-3) )
# sim.createParticles(geoHCP)

geoHydrate = RandomBoxPacker(
   minRadius = r_p,
   maxRadius = r_p+1e-5,
   cubicPackRadius = 2*(r_p+3e-5)+r_p,
   maxInsertFails = 2000,
   bBox = BoundingBox(
      Vec3(-4e-3,-4e-3,-4e-3),
      Vec3(4e-3,4e-3,4e-3)
   ),
   circDimList = [False, False, False],
   tolerance = 0.001
)

geoHydrate.generate()
geoRandomBlock_particles = geoHydrate.getSimpleSphereCollection()
sim.createParticles(geoRandomBlock_particles)


# geoHCP2 = HexagBlock (
#    dimCount = [count,count,3],
#    radius = r_s
# )

# geoHCP2.translate ( translation = Vec3(0,0,4e-3))
# sim.createParticles(geoHCP2)

plist = []
# for p1 in geoHCP:
#    p1.setTag(1)
#    plist.append(p1.getId())
# maxId = max(plist)

for i,p2 in enumerate(geoRandomBlock_particles):
   p2.setTag(i)
   plist.append(p2.getId())
maxId = max(plist)

# for p2 in geoRandomBlock_particles:
#    p2.setTag(2)
#    oldId_2 = p2.getId() + 1
#    p2.setId(oldId_2 + maxId)
#    plist.append(p2.getId())
# maxId = max(plist)

# for p3 in geoHCP2:
#    p3.setTag(3)
#    oldId_3 = p3.getId() + 1
#    p2.setId(oldId_3 + maxId)
#    plist.append(p3.getId())

# sim.setParticleDensity(tag=1, mask=-1, Density=density)
sim.setParticleDensity(tag=2, mask=-1, Density=density)
# sim.setParticleDensity(tag=3, mask=-1, Density=density)


sim.createInteractionGroup(
    RotFrictionPrms(
        name="friction",
        normalK=normalK_sand,
        dynamicMu=dynamicMu,
        staticMu=staticMu,
        shearK=shearK_sand,
        scaling=True,
        rigid=False,
        meanR_scaling=True
    )
)

# #add bond to first layer sand
# sim.createConnections(
#     ConnectionFinder(
#         maxDist=r_p,
#         bondTag=0,
# 	    pList=geoHCP
#     )
# )
# sim.createInteractionGroup(
#     BrittleBeamPrms(
#         name="ppbond_sand_1",
#         youngsModulus=1.1e10,
#         poissonsRatio=0.31,
#         cohesion=6.73e6,
#         tanAngle=bondTanAngle,
#         tag=0,
#         meanR_scaling=True,
#         truncated=2*r_s*2.86e8,
# 	    beta1=1.0,
#         beta2=1.0
#     )
# )

sim.createExclusion(
    interactionName1="ppbond_sand_1",
    interactionName2="friction"
)

# #add bond to second layer sand
# sim.createConnections(
#     ConnectionFinder(
#         maxDist=r_p,
#         bondTag=1,
# 	    pList=geoHCP2
#     )
# )
# sim.createInteractionGroup(
#     BrittleBeamPrms(
#         name="ppbond_sand_2",
#         youngsModulus=1.1e10,
#         poissonsRatio=0.31,
#         cohesion=6.73e6,
#         tanAngle=bondTanAngle,
#         tag=1,
#         meanR_scaling=True,
#         truncated=2*r_s*2.86e8,
# 	    beta1=1.0,
#         beta2=1.0
#     )
# )

# sim.createExclusion(
#     interactionName1="ppbond_sand_2",
#     interactionName2="friction"
# )

#add bond to hydrate
sim.createConnections(
    ConnectionFinder(
        maxDist=r_p*50,
        bondTag=2,
	    pList=geoRandomBlock_particles
    )
)
sim.createInteractionGroup(
    BrittleBeamPrms(
        name="ppbond_hydrate",
        youngsModulus=bondModulus*k1*Sh,
        poissonsRatio=bondPoissonRatio,
        cohesion=bondCohesion*k0*Sh,
        tanAngle=bondTanAngle,
        tag=2,
        meanR_scaling=True,
        truncated=maxTensile,
	    beta1=1.0,
        beta2=1.0
    )
)
sim.createExclusion(
    interactionName1="ppbond_hydrate",
    interactionName2="friction"
)
#add bond between first layer and hydrate
# sim.createInteractionGroupTagged(
#     BrittleBeamPrms(
#         name="ppbond_sand_hydrate1",
#         youngsModulus=1.1e10,
#         poissonsRatio=0.31,
#         cohesion=6.73e6,
#         tanAngle=bondTanAngle,
#         tag=3,
#         meanR_scaling=True,
#         truncated=2*r_s*2.86e8,
# 	    beta1=1.0,
#         beta2=1.0
#     ),
#     tag1=1,
#     mask1=-1,
#     tag2=2,
#     mask2=-1

# )


# sim.createExclusion(
#     interactionName1="ppbond_sand_hydrate1",
#     interactionName2="friction"
# )

# #add bond between second layer and hydrate
# sim.createInteractionGroupTagged(
#     BrittleBeamPrms(
#         name="ppbond_sand_hydrate1",
#         youngsModulus=1.1e10,
#         poissonsRatio=0.31,
#         cohesion=6.73e6,
#         tanAngle=bondTanAngle,
#         tag=3,
#         meanR_scaling=True,
#         truncated=2*r_s*2.86e8,
# 	    beta1=1.0,
#         beta2=1.0
#     ),
#     tag1=2,
#     mask1=-1,
#     tag2=3,
#     mask2=-1
# )

# print "Initial number particles:",sim.getNumParticles()
# for n in range(2000):
#    sim.runTimeStep()
#    if (n%100==0): print "number particles:",sim.getNumParticles()
# print "Final number particles:",sim.getNumParticles()
# sim.exit()

# sim.createInteractionGroup(
#     LinDampingPrms(
#         name="damping",
#         viscosity=viscosity_damp_nr,
#         maxIterations=100
#     )
# )

# sim.createFluidForce(
#     FluidForcePrms(
#         name="lbm"
#     )
# )


# add a back wall to the model:
sim.createWall(
    name="top",
    posn=Vec3(0,0,1.35e-3+6*r_s),
    normal=Vec3(0,0,-1)
)

# add a back wall to the model:
sim.createWall(
    name="bottom",
    posn=Vec3(0,0,-1.35e-3-6*r_s),
    normal=Vec3(0,0,+1)
)

sim.createInteractionGroup(
    NRotBondedWallPrms(
        name="ptop",
        wallName="top",
        normalK=1e+8,
        particleTag=2
    )
)
sim.createInteractionGroup(
    NRotBondedWallPrms(
        name="pbottom",
        wallName="bottom",
        normalK=1e+8,
        particleTag=2
    )
)

#add a wall loader to move the top wall:
wall_loader1 = WallLoaderRunnable(
   LsmMpi = sim,
   wallName = "top_wall",
   vPlate = Vec3 (0.0,0.0,-4e-3),
   startTime = 0,
   rampTime = 50000
)
sim.addPreTimeStepRunnable (wall_loader1)

#add a wall loader to move the bottom wall:
wall_loader2 = WallLoaderRunnable(
   LsmMpi = sim,
   wallName = "bottom_wall",
   vPlate = Vec3 (0.0,0.0,4e-3),
   startTime = 0,
   rampTime = 50000
)
sim.addPreTimeStepRunnable (wall_loader2)



# # # add a wall loader to move the back wall:
# tensile = tensileClass(
#     LsmMpi=sim,
#     speed=speed,
#     endTs=endTs,
#     iniD=iniD,
#     r_p=r_p
# )
# sim.addPreTimeStepRunnable(tensile)


# add a CheckPointer to store simulation data:

sim.createFieldSaver (
    InteractionScalarFieldSaverPrms(
        interactionName="ppbond_hydrate",
        fieldName="count",
        fileName=out_dat_dir + "/nbonds.dat",
        fileFormat="SUM",
        beginTimeStep=0,
        endTimeStep=endTs,
        timeStepIncr=incTs_data
   )
)

sim.createFieldSaver (
    InteractionVectorFieldSaverPrms(
        interactionName="ppbond_hydrate",
        fieldName="normal_force",
        fileName=out_dat_dir + "/" + force_dir + "/nforce",
        fileFormat="RAW2",
        beginTimeStep=0,
        endTimeStep=endTs,
        timeStepIncr=incTs_data
   )
)



# # sim.createFieldSaver(
# #     ParticleVectorFieldSaverPrms(
# #         fieldName="displacement",
# #         fileName=out_dat_dir + "/" + disp_dir + "/displacement",
# #         fileFormat="RAW2",
# #         beginTimeStep=0,
# #         endTimeStep=endTs,
# #         timeStepIncr=incTs_data
# #     )
# # )

# # sim.createFieldSaver(
# #     ParticleVectorFieldSaverPrms(
# #         fieldName="force",
# #         fileName=out_dat_dir + "/" + force_dir + "/force",
# #         fileFormat="RAW2",
# #         beginTimeStep=0,
# #         endTimeStep=endTs,
# #         timeStepIncr=incTs_data
# #     )
# # )

# # sim.createFieldSaver(
# #     ParticleVectorFieldSaverPrms(
# #         fieldName="position",
# #         fileName=out_dat_dir + "/" + position_dir + "/position",
# #         fileFormat="RAW2",
# #         beginTimeStep=0,
# #         endTimeStep=endTs,
# #         timeStepIncr=incTs_data
# #     )
# # )

# # sim.createFieldSaver(
# #     ParticleScalarFieldSaverPrms(
# #         # fieldName="sigma_xx_2d",
# #         fieldName="sigma_d",
# #         fileName=out_dat_dir + "/" + sigmad_dir + "/sigma_d",
# #         fileFormat="RAW_WITH_POS_ID",
# #         beginTimeStep=0,
# #         endTimeStep=endTs,
# #         timeStepIncr=incTs_data
# #     )
# # )


# # # create a FieldSaver to wall forces:
# # force_saver = WallVectorFieldSaverPrms(
# #     wallName=["top", "bottom"],
# #     fieldName="Force",
# #     fileName=out_dat_dir+"/out_Force.dat",
# #     fileFormat="RAW_SERIES",
# #     beginTimeStep=0,
# #     endTimeStep=endTs,
# #     timeStepIncr=incTs_data
# # )
# # sim.createFieldSaver(force_saver)
# add a CheckPointer to store simulation data:
sim.createCheckPointer(
    CheckPointPrms(
        fileNamePrefix="snapshot",
        beginTimeStep=0,
        endTimeStep=10,
        # endTimeStep=endTs,
        timeStepIncr=10
        #timeStepIncr=incTs_snap
    )
)




# execute the simulation:
sim.run()
