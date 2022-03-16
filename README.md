# Contents

- [Packing](#Packing)
- [Tensile](#Tensile)

# Packing
## Table of Contents

- [About](#about)
- [Usage](#usage)

## About <a name = "about"></a>

This program is for packing hydrate-particles in cylindrical shape.

## Usage <a name = "usage"></a>
The basical usage is as follows:
1. [Make a cylindrical mesh ](#1)
2. [Set the condition in condition.csv](#2)
3. [Run packing.py](#3)
4. [Output a geo file](#4)

### Make a cylindrical mesh <a name = "1"></a>
Since this program packs particles into a cylindrical mesh, a cylindrical mesh must be created in advance.


In creating a mesh, [gmsh](http://gmsh.info) is used.
Please install gmsh to the PC.

After installe gmsh, open ./mesh/cylinder.geo by a texteditor.
Radius and height of a cylinder are defined in the top of .geo file like below. 
```cylinder.geo
//radius [m]
r = 0.0025;

//height [m]
h = 0.0027;
```

Then, open cylinder.geo by gmsh.
After opening geo file, select mesh->2D, then cylinder is splitted as a mesh. And you can increase the number of nodes by clicking the "Refine by splitting".

<img width="747" alt="mesh" src="https://user-images.githubusercontent.com/50572759/94533443-15351c80-027a-11eb-9055-ebd4e635fc58.png">

Next, export the mesh. Select File->Export, and fill the name as "cylinder" and chose "INRIA Medit" in Save As. Click Save button, then "cylinder" will be output.
<img width="953" alt="save" src="https://user-images.githubusercontent.com/50572759/94533573-444b8e00-027a-11eb-91c5-fa916e9183d9.png">

"cylinder" has to be converted so as esysparticle can read.
medit2timesh.py can do this.
Please just run the below command:
```
python medit2timesh.py
```
Then, cylinder.lsm, which is a mesh file for esysparticle, will be output.




### Set the condition in condition.csv <a name = "2"></a>


|variable|explanation|
|---|---|
|endTs|the time step stops the simulation|
|incTs_snap|the interval time step for output of snapshot|
|dT[s]|delta t[s]|
|R[m]|the radius of the cylinder[m]|
|H[m]|the height of the cylinder[m]|
|r_h[m]|the radius of hydrates[m]|
|normalK_wall[N/m]|Normal k of wall|
|normalK_sand[N/m]|Normal k of particles|
|shearK_sand[N/m]|Shear k of particles|
|dynamicMu|Dynamic frictional coefficient|
|staticMu|Static frictional coefficient|
|porosity|porosity|
|viscosity_damp|Viscousity|
|initialfactor|the factor multipled to the particle size at the start of the calculation|
|threshold|file name of threshold file|


#### Set the threshold for the overlap in threshold.csv <a name = "2"></a>
|index|threshold[m]|
|---|---|
|...|...|

"index" is used for the filename of snapshot files for each threshold.

### Run packing.py <a name = "3"></a>
To run the packing program, please run the below command.
```
mpirun -np #proc esysparticle packing.py
```

In this program, the radius expansion technique is used to prevent particle overlap. First, the number of particles satisfying the target porosity is calculated and randomly placed in a cylinder. The calculation is started with each particle size multiplied by a certain magnification(initialfactor) to reduce it, and the size of the particles is returned to their original size over time. After the particles have settled, the calculation is finished.
The standard output shows "scaling" when the particle is expanding, and "easing" when it is doing nothing and waiting for the particle to settle down.

The function, getMinMaxDistance(), returns the min. and max. distance of all particles. The distance means (the distance between two particles from their center) - (sum of the radius of two particles).
So, if min. distance is below the threshold, the "easing" will be stopped.
After the min. distance is below the threshold, snapshot_threshold_idx#_t=#_#.txt will be generated. Convert this into geo file and use it in your tensile program.

### Output a geo file <a name = "4"></a>
To use the particles packed by this program in the tensile program, you need to output a geo file by
```
dump2geo -i snapshot_packing -o "outputname" -rot -t "target number" 1 1
```
The usage of dump2geo is the same to dump2vtk.

After created geo file, use it in tensile program.


# Tensile

## About <a name = "about"></a>

Tensile test for methane hydrates using esysparticle.

## Usage <a name = "usage"></a>

This project consists of the followings:
 - tensile.py : a main program for the tensile test
 - tensileRunnable.py : a runnable class for moving particles
 - conditions.csv : a input file
 - plot.py: a program for plots.

To run this program, execute the following command.
```
mpirun -np 2 esysparticle tensile.py
```
After a computation, results will be output in "out_data" directory.
Its contents are :
 - nbonds.dat, the number of bonds.
 - wall_Force.dat, the force acting on the walls.
 - wall_Position.dat, the position of the wall the walls.

The bond will be broken when the tensile stress met the limitation, 
but the compuation will continue until the strain reaches the value set in conditions.csv(endstrain).
So, check if it should be stopped by plot the stress-strain curve.

### conditions.csv
The contens in conditions.csv are explained below.

|variable|explanation|
|---|---|
|dT[s]|delta t[s]|
|incT[s]|the interval for output of data[s]|
|incT_snap[s]|the interval for output of snapshot|
|speed[m/s]|the speed of the particle[m/s]|
|H[m]|the height of the cylinder|
|R[m]|the radius of the cylinder|
|r_p[m]|the radius of particles|
|endstrain|max of the strain.|
|density[kg/m^3]|density of particles|
|maxTensile[Pa]|max tensile stress that the bond will break|
|Sh|the saturation of methane hydrates|
|bondModulus|Young's modulus of methane hydrates|
|bondPoissonRatio|Poisson ratio of methane hydrates|
|bondCohesion|Cohesion of methane hydrates|
|bondTanAngle|Angle of MohrCoulomb criterion. Set as 4.0|
|beta1|a param for BrittleBeamIG. Set as 1.0|
|beta2|a param for BrittleBeamIG. Set as 1.0|
|k0|a cohesion coefficient for MHbond. Set as 1.0|
|k1|a Young's modules coefficient for MHbond. Set as 1.0|
|normalK_wall[N/m]|Normal kn of wall. Set as 10.0e9|
|normalK_sand[N/m]|Normal kn of particles. Set as 3.5e9|
|kn/ks|a ratio of kn to ks of particles. Set as 10.0|
|dynamicMu|Dynamic frictional coefficient. Set as 0.001|
|staticMu|Static frictional coefficient. Set as 0.4|
|viscosity_damp_nr|Viscousity. Set as 0.00|


### plot
plot.py is for plotting some figures.
Just exexute the following command.
```
python plot.py
```
Then, some figures are output in "/fig" directory. 
This program uses matplotlib. 
