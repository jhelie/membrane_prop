# membrane_prop
Python script to analyse properties of non planar membranes.
To print the help below: ```python membrane_prop.py --help```

```
[ DESCRIPTION ]
 
This script calculates density profiles along the local normal of bilayer.

Particles for which to calculate the density are defined via a text file specified
by --particles. Charges density can similarly be calculated using another text file
via --charges.

particles density calculation
-----------------------------
Densities are calculated based on transmembrane cylinders oriented along the local
bilayer normal. Particles present in the cylinder are binned into slices based on
their distance along the local normal to the local center of the lipid bilayer.
The cylinders radius and slices thickness are controlled via --slices_radius and
--slices_thick respectively.

charge density calculation
--------------------------
Charge densities are calculated similarly to particles densities with the exception
that they are also binned by groups whereby the charges of the groups components are
summed in order to calculate the net charge contribution of each group.

bilayer sampling
---------------
A 3D grid is applied to downsample the set of points used to describe the bilayer (see
leaflets identification options). This is achieved using Python bindings for the Point
Cloud Library (http://pointclouds.org/). The dimension of the voxels can be controlled
in each dimension. Voxel are then approximated by the center of geometry of the points
they contain. Particles and charge densities are calculated for each populated voxel.


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - numpy
 - sscipy
 - networkX (if option --algorithm is set to 'min' or 'cog')
 - PCL C library and corresponding python bindings


[ NOTES ]

1. You can specify the particles for which to plot the density by supplying a file
   via the --particles option. Each line of this file should follow the following
   format (without quotation marks):
    -> 'group,label,colour,MDAnalysis selection string'
 
   where 'group' is used to normalise the densities of several particles. By default
   following densities definition are used:
    -> peptide,peptide,#262626,protein
    -> CHOL,CHOL,#bdbdbd,resname CHOL and name ROH
    -> POPC,POPC,#41ab5d,resname POPC and name PO4
    -> POPE,POPE,#6a51a3,resname POPE and name PO4
    -> POPS,POPS,#cc4c02,resname POPS and name PO4
    -> water,water,#1d91c0,resname W
    -> ions,Na+,#7bccc4,name NA+
    -> ions,Cl-,#fa9fb5,name CL-

2. You can specify which particles to take into account for the calculation of the total
   charge density by supplying a file via the --charges option. Each line of this file
   should follow the format (without quotation marks):
    -> 'group_name,colour,charge_name,charge,MDAnalysis selection string for charge'

   The absolute charge for each group will be plotted on the charge density profile. The
   group colour must be specified for each charge. By default the charged are defined as
   follows:
    -> ions,#52A3CC,Na+,1,resname NA+
    -> ions,#52A3CC,CL-,-1,resname CL-
    -> lipids,#b2182b,-1,phosphate,name PO4
    -> lipids,#b2182b,1,amine_choline,name NH3 or name NC3
  
   If you do not want to calculate charge density, just use: '--charges no'

3. Colours in the above definition files can be specified using single letter code (rgbcmykw),
   hex code  or the name of a colour map (see the matplotlib website for a list of available
   colour maps). In case a colour map is used, its name must be specified as the only colour.

4. Identification of the bilayer leaflets can be controlled via the 3 options below. Note that
   by default the gro file -f is only used as a topology file and the 1st frame of the -x file
   is used to identify leaflets. If you wish to use the gro file instead, for instance in the
   case that the 1st frame of the xtc is not flat, you need to specify the --use_gro flag but
   be warned that this might take a few minutes longer on large systems.

   (a) beads
    By default, the particles taken into account to define leaflet are:e
    -> name PO4 or name PO3 or name B1A
   
    Note that only lipids which contain one of the beads mentioned in the selection string
    will be taken into account. If you wish to specify your own selection string (e.g. to
    choose different beads or add a bead not in the default list in order to take into
    account a particular lipid specie) you can do so by supplying a file via the --beads
    option. This file should contain a single line that can be passed as the argument
    to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
     -> name PO4 or name PO3 or name B1A or name AM1
        
   (b) leaflet finding method
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
     -> '--leaflet 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflet large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the 1st frame of the xtc
    file supplied in order to get a meaningful outcome. 

   (c) flipflopping lipids
    In case lipids flipflop during the trajectory, a file listing them can be supplied
    with the --flipflops option. Each line of this file should follow the format:
     -> 'resname,resid,starting_leaflet,z_bead'
    where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
    z_bead is used to track the position of the lipid.
    If flipflopping lipids are not specified they may add significant noise to results as
    they prevent the correct identification of TM clusters.
    Note that, even when specified, flipflopping lipids will be taken into account when
    calculating densities and charges.   

5. There are 3 possible options to determine the local normal to the bilayer. These are
   controlled with the flags --normal and --normal_d:
   (a) 'z': the bilayer is assumed flat in the x,y plane and the z axis is taken to be the
    normal. Best for systems without significant curvature and local deformations. In this
    case the --normal_d flag is ignored.

   (b) 'cog': in this case, within each voxel, --beads particles are identified in each leaflet.
    The local normal is then considered to be the vector going from the cog of the lower ones
    to the cog of the upper ones. 

   (c) 'svd': in this case --beads particles are identified in each voxel and the normal of
    the best fitting  plane to these particles is obtained by performing a singular value
    decomposition of their coordinate matrix.

 
[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc or another .gro]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		10	: process every t-frames
 
Density profile options
-----------------------------------------------------
--particles		: definition of particles, see note 1
--charges		: definition of charged particles, note 2
--slices_thick	0.5 	: z thickness of the slices (Angstrom)
--slices_radius	30 	: radius of the slices (Angstrom)
--slices_range	40 	: distance to consider along local normal on either side of the bilayer center

Surface description
-----------------------------------------------------
--voxel_x	50	: voxel dimension along x (Angstrom)
--voxel_y	50	: voxel dimension along y (Angstrom)
--voxel_z	50	: voxel dimension along z (Angstrom)
--voxel_nb	10	: min nb of points within voxels to consider (requires modified PCL bindings)
--normal	svd	: local normal to bilayer ('z', 'cog' or 'svd'), see note 5
--normal_d	50	: distance of points to take into account for local normal, see note 5

Leaflets identification  
-----------------------------------------------------
--beads			: leaflet identification technique, see note 4(a)
--flipflops		: input file with flipflopping lipids, see note 4(c)
--leaflets	optimise: leaflet identification technique, see note 4(b)
--use_gro		: use gro file instead of xtc, see note 4(b)

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
```
