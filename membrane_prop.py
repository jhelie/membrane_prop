#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.4d"
parser = argparse.ArgumentParser(prog = 'membrane_prop', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@gmail.com)
git: https://github.com/jhelie/membrane_prop
**********************************************

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

Convexity plots
-----------------------------------------------------
--nxmax			: upper boundary of y axis for dnx/dx convexity plot
--nxmin			: lower boundary of y axis for dnx/dx convexity plot
--nymax			: upper boundary of y axis for dny/dy convexity plot
--nymin			: lower boundary of y axis for dny/dy convexity plot
--nzmax			: upper boundary of y axis for dnz/dz convexity plot
--nzmin			: lower boundary of y axis for dnz/dz convexity plot
--nbx			: nb ticks on x axis
--nby			: nb ticks on y axis

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)

#density profile options
parser.add_argument('--particles', nargs=1, dest='particlesfilename', default=['mine'], help=argparse.SUPPRESS)
parser.add_argument('--charges', nargs=1, dest='chargesfilename', default=['mine'], help=argparse.SUPPRESS)
parser.add_argument('--slices_thick', nargs=1, dest='slices_thick', default=[0.5], type=float, help=argparse.SUPPRESS)
parser.add_argument('--slices_radius', nargs=1, dest='slices_radius', default=[30], type=float, help=argparse.SUPPRESS)
parser.add_argument('--slices_range', nargs=1, dest='slices_range', default=[40], type=float, help=argparse.SUPPRESS)

#surface description
parser.add_argument('--voxel_x', nargs=1, dest='voxel_x', default=[50], type=float, help=argparse.SUPPRESS)
parser.add_argument('--voxel_y', nargs=1, dest='voxel_y', default=[50], type=float, help=argparse.SUPPRESS)
parser.add_argument('--voxel_z', nargs=1, dest='voxel_z', default=[50], type=float, help=argparse.SUPPRESS)
parser.add_argument('--voxel_nb', nargs=1, dest='voxel_nb', default=[10], type=int, help=argparse.SUPPRESS)
parser.add_argument('--normal', dest='normal', choices=['z','cog','svd'], default='svd', help=argparse.SUPPRESS)
parser.add_argument('--normal_d', nargs=1, dest='normal_d', default=[50], type=float, help=argparse.SUPPRESS)

#lipids identification options
parser.add_argument('--beads', nargs=1, dest='beadsfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)
parser.add_argument('--use_gro', dest='use_gro', action='store_true', help=argparse.SUPPRESS)

#convexity plots
parser.add_argument('--nxmax', nargs=1, dest='nxmax', default=[-1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nxmin', nargs=1, dest='nxmin', default=[-1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nymax', nargs=1, dest='nymax', default=[-1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nymin', nargs=1, dest='nymin', default=[-1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nzmax', nargs=1, dest='nzmax', default=[-1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nzmin', nargs=1, dest='nzmin', default=[-1], type=float, help=argparse.SUPPRESS)
parser.add_argument('--nbx', nargs=1, dest='nbx', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('--nby', nargs=1, dest='nby', default=[-1], type=int, help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]
#density profile options
args.particlesfilename = args.particlesfilename[0]
args.chargesfilename = args.chargesfilename[0]
args.slices_range = args.slices_range[0]
args.slices_thick = args.slices_thick[0]
args.slices_radius = args.slices_radius[0]
#surface description
args.voxel_x = args.voxel_x[0]
args.voxel_y = args.voxel_y[0]
args.voxel_z = args.voxel_z[0]
args.voxel_nb = args.voxel_nb[0]
args.normal_d = args.normal_d[0]
#lipids identification options
args.beadsfilename = args.beadsfilename[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]
args.selection_file_ff = args.selection_file_ff[0]
#convexity plots
args.nxmax = args.nxmax[0]
args.nxmin = args.nxmin[0]
args.nymax = args.nymax[0]
args.nymin = args.nymin[0]
args.nzmax = args.nzmax[0]
args.nzmin = args.nzmin[0]
args.nbx = args.nbx[0]
args.nby = args.nby[0]

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#point cloud library
try:
	import pcl
except:
	print "Error: you must install the pcl C library and the corresponding Python bindings."
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.selection_file_ff != "no" and not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.beadsfilename != "no" and not os.path.isfile(args.beadsfilename):
	print "Error: file " + str(args.beadsfilename) + " not found."
	sys.exit(1)
if args.particlesfilename != "mine" and not os.path.isfile(args.particlesfilename):
	print "Error: file " + str(args.particlesfilename) + " not found."
	sys.exit(1)
if args.chargesfilename != "no" and args.chargesfilename != "mine" and not os.path.isfile(args.chargesfilename):
	print "Error: file " + str(args.chargesfilename) + " not found."
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)
if args.normal != 'z' and args.normal_d <= 0:
	print "Error: --normal_d (" + str(args.normal__d) + " AA) should be greater 0. Or choose 'z' for --normal, see note 5."
	sys.exit(1)

if args.xtcfilename == "no":
	if args.use_gro:
		print "Error: you can only use the --use_gro file if you specify an xtc file."
		sys.exit(1)
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=========================================================================================
# process options
#=========================================================================================

global lipids_ff_nb
global bins_nb
global bins_nb_max
global bins_labels
global slice_volume
lipids_ff_nb = 0
bins_nb = int(np.floor(args.slices_range/float(args.slices_thick))) 			#actually it's twice that as (-bins_nb,bins_nb) has to be filled
bins_nb_max = bins_nb
bins_labels = [str((n+0.5)*args.slices_thick) for n in range(-bins_nb,bins_nb)]
slice_volume = 2 * math.pi * args.slices_radius * args.slices_thick

#leaflet identification
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
		sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "membrane_prop_" + args.grofilename[:-4].split('/')[-1]
	else:
		args.output_folder = "membrane_prop_" + args.xtcfilename[:-4].split('/')[-1]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	#create folders
	os.mkdir(args.output_folder)	
	os.mkdir(args.output_folder + "/angle")
	os.mkdir(args.output_folder + "/density")
	if args.chargesfilename != "no":
		os.mkdir(args.output_folder + "/charge")
	#os.mkdir(args.output_folder + "/pcl") #TO DO: save point clouds data structures...

	#create log
	filename_log = os.getcwd() + '/' + str(args.output_folder) + '/membrane_prop.log'
	output_log = open(filename_log, 'w')		
	output_log.write("[membrane_prop v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python membrane_prop.py"
	for c in sys.argv[1:]:
		tmp_log += " " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	
	#copy input files
	if args.selection_file_ff != "no":
		shutil.copy2(args.selection_file_ff,args.output_folder + "/")	
	if args.beadsfilename != "no":
		shutil.copy2(args.beadsfilename,args.output_folder + "/")
	if args.particlesfilename != "mine":
		shutil.copy2(args.particlesfilename,args.output_folder + "/")
	if args.chargesfilename != "no" and args.chargesfilename != "mine":
		shutil.copy2(args.chargesfilename,args.output_folder + "/")

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():													#DONE

	global leaflet_sele_string
	global leaflet_sele_string_original

	#set default beads
	leaflet_sele_string = "name PO4 or name PO3 or name B1A"
	
	#use users input
	if args.beadsfilename != "no":
		with open(args.beadsfilename) as f:
			lines = f.readlines()
		if len(lines) > 1:
			print "Error: the file " + str(args.beadsfilename) + " should conly ontain 1 line (" + str(len(lines)) + " found), see note 2(a)."
			sys.exit(1)
		else:
			if lines[0][-1] == "\n":
				lines[0] = lines[0][:-1]
			leaflet_sele_string = lines[0]

	leaflet_sele_string_original = leaflet_sele_string
	
	return
def set_particles():													#DONE

	global particles_def
	global particles_def_pres
	global particles_groups
	
	particles_def = {k: {} for k in ["group","colour","sele","sele_string"]}
	colours_map = 'custom'
	colours_map_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']
	
	#use default particles definition
	#--------------------------------
	if args.particlesfilename == "mine":
		particles_def["labels"] = ["peptide","POPC","POPE","POPS","CHOL","water","Na+","Cl-"]
				
		#peptide
		particles_def["group"]["peptide"] = "peptide"
		particles_def["colour"]["peptide"] = "#262626"					#very dark grey
		particles_def["sele_string"]["peptide"] = "protein"
	
		#lipids
		particles_def["group"]["CHOL"] = "CHOL"
		particles_def["group"]["POPC"] = "POPC"
		particles_def["group"]["POPE"] = "POPE"
		particles_def["group"]["POPS"] = "POPS"
		particles_def["colour"]["CHOL"] = "#bdbdbd"						#light grey
		particles_def["colour"]["POPC"] = "#41ab5d"						#dark green
		particles_def["colour"]["POPE"] = "#6a51a3"						#dark purple
		particles_def["colour"]["POPS"] = "#cc4c02"						#dark orange
		particles_def["sele_string"]["POPC"] = "resname POPC and name PO4"		
		particles_def["sele_string"]["POPE"] = "resname POPE and name PO4"
		particles_def["sele_string"]["POPS"] = "resname POPS and name PO4"
		particles_def["sele_string"]["CHOL"] = "resname CHOL and name ROH"

		#solvent
		particles_def["group"]["water"] = "water"
		particles_def["colour"]["water"] = "#1d91c0"					#dark cyan
		particles_def["sele_string"]["water"] = "resname W"

		#ions
		particles_def["group"]["Na+"] = "ions"
		particles_def["group"]["Cl-"] = "ions"
		particles_def["colour"]["Cl-"] = "#fa9fb5"						#light red
		particles_def["colour"]["Na+"] = "#7bccc4"						#light blue
		particles_def["sele_string"]["Na+"] = "name NA+"
		particles_def["sele_string"]["Cl-"] = "name CL-"
	
	#use user's particles definition
	#-------------------------------	
	else:
		particles_def["labels"] = []
		with open(args.particlesfilename) as f:
			lines = f.readlines()
		for l_index in range(0, len(lines)):
			line = lines[l_index]
			if line[-1] == "\n":
				line = line[:-1]
			l_content = line.split(',')
			if len(l_content) != 4:
				print "Error: wrong format on line " + str(l_index + 1) + " of " + str(args.particlesfilename) + ", see DESCRIPTION in membrane_prop --help."
				sys.exit(1)
			tmp_group = l_content[0]
			tmp_label = l_content[1]
			tmp_color = l_content[2]
			tmp_selec = l_content[3]
			if tmp_label in particles_def["labels"]:
				print "Error: the particle " + str(tmp_label) + " is defined mor than once. Check your --particles option."
				sys.exit(1)
			particles_def["labels"].append(tmp_label)
			particles_def["group"][tmp_label] = tmp_group
			particles_def["colour"][tmp_label] = tmp_color
			particles_def["sele_string"][tmp_label] = tmp_selec

		#generate colours if necessary
		if len(np.unique(particles_def["colour"].values())) == 1:
			if np.unique(particles_def["colour"].values())[0] in colours_map_possible:
				colours_map = np.unique(particles_def["colour"].values())[0]
			else:
				print "Error: either the same colour was specified for all particles or the colour map '" + str(np.unique(particles_def["colour"].values())[0]) + "' is not valid."
				sys.exit(1)
		if colours_map != "custom":
			tmp_cmap = cm.get_cmap(colours_map)
			part_nb = len(particles_def["labels"])
			part_colours_value = tmp_cmap(np.linspace(0, 1, part_nb))
			for part_index in range(0, part_nb):
				particles_def["colour"][particles_def["labels"][part_index]] = part_colours_value[part_index]

	#initialise presence of particles
	#--------------------------------
	particles_def_pres = {part: False for part in particles_def["labels"]}

	#build particle groups for normalisation
	#---------------------------------------
	particles_groups = {k: [part for part, g in particles_def["group"].items() if g == k] for k in np.unique(particles_def["group"].values())}	
		
	return
def set_charges():														#DONE
	
	global charges_groups
	global charges_colours
	global charges_groups_pres
	global charges_groups_pres_q
	charges_groups = {}
	charges_colours = {}
	charges_colours["total"] = '#262626'								#very dark grey
	colours_map = 'custom'
	colours_map_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']
	
	#use default
	#-----------
	if args.chargesfilename == "mine":
		#ions
		charges_colours["ions"] = "#52A3CC"								#cyan colour
		charges_groups["ions"] = {}
		charges_groups["ions"]["names"] = ["Na+","Cl-"]
		charges_groups["ions"]["values"] = {}
		charges_groups["ions"]["values"]["Na+"] = 1
		charges_groups["ions"]["values"]["Cl-"] = -1
		charges_groups["ions"]["sele"] = {}
		charges_groups["ions"]["sele_string"] = {}
		charges_groups["ions"]["sele_string"]["Na+"] = "name NA+"
		charges_groups["ions"]["sele_string"]["Cl-"] = "name CL-"
		
		#lipids
		charges_colours["lipids"] = "#b2182b"							#dark red
		charges_groups["lipids"] = {}
		charges_groups["lipids"]["names"] = ["PO4","NH3-NC3"]			#bear in mind that for PO4 xtc the NH3/NC3 are not there to counterbalance the charge...
		charges_groups["lipids"]["values"] = {}
		charges_groups["lipids"]["values"]["PO4"] = -1
		charges_groups["lipids"]["values"]["NH3-NC3"] = 1
		charges_groups["lipids"]["sele"] = {}
		charges_groups["lipids"]["sele_string"] = {}
		charges_groups["lipids"]["sele_string"]["PO4"] = "name PO4"
		charges_groups["lipids"]["sele_string"]["NH3-NC3"] = "name NH3 or name NC3"
	
	#use user supplied
	#-----------------
	else:	
		with open(args.chargesfilename) as f:
			lines = f.readlines()
		for l_index in range(0, len(lines)):
			line = lines[l_index]
			if line[-1] == "\n":
				line = line[:-1]
			l_content = line.split(',')
			if len(l_content) != 5:
				print "Error: wrong format on line " + str(l_index + 1) + " of " + str(args.chargesfilename) + ", see DESCRIPTION in membrane_prop --help."
				sys.exit(1)
			tmp_group = l_content[0]
			tmp_colour = l_content[1]
			tmp_name = l_content[2]
			tmp_value = float(l_content[3])
			tmp_sele = l_content[4]
			if tmp_group not in charges_groups.keys():
				charges_groups[tmp_group] = {k: {} for k in ["values","sele","sele_string"]}
				charges_groups[tmp_group]["names"] = []
				charges_colours[tmp_groups] = tmp_colour
			charges_groups[tmp_group]["names"].append(tmp_name)
			charges_groups[tmp_group]["values"][tmp_name] = tmp_value
			charges_groups[tmp_group]["sele_string"][tmp_name] = tmp_sele

		#generate colours from colour map if necessary
		if len(np.unique(charges_colours.values())) == 1:
			if np.unique(charges_colours.values())[0] in colours_map_possible:
				colours_map = np.unique(charges_colours.values())[0]
			else:
				print "Error: either the same colour was specified for all charge groups or the colour map '" + str(np.unique(residues_colours.values())[0]) + "' is not valid."
				sys.exit(1)
		if colours_map != "custom":
			tmp_cmap = cm.get_cmap(colours_map)
			charges_groups_keys = charges_groups.keys()
			charges_groups_nb = len(charges_groups_keys)
			charges_colours_value = tmp_cmap(np.linspace(0, 1, charges_groups_nb))
			for q_index in range(0, charges_groups_nb):
				charges_colours[charges_groups_keys[q_index]] = charges_colours_value[q_index]

	#initialise presence of charges groups
	#-------------------------------------
	charges_groups_pres = {charge_g: False for charge_g in charges_groups.keys()}
	charges_groups_pres_q = {charge_g: {q: False for q in charges_groups[charge_g]["names"]} for charge_g in charges_groups.keys()}

	return
def load_MDA_universe():												#DONE
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global nb_frames_to_process
	global f_start
	global f_end
	global residues_list
	global water_pres
	global water_sele
	f_start = 0
	f_end = 1
	
	#load universe
	#-------------
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		if args.use_gro:
			global U_gro
			U_gro = Universe(args.grofilename)
		U = Universe(args.grofilename, args.xtcfilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes
		U.trajectory.rewind()
		
		#case: xtc is a gro
		if args.xtcfilename[-3:] == "gro":
			if not args.use_gro:
				print "Error: the trajectory is a gro file but the --use_gro option has not been specified: either use --use_gro or run the script using -f only (see note 4)"
				sys.exit(1)
			frames_to_process = [0]
			nb_frames_to_process = 1
			args.xtcfilename = "no"
		#case: xtc is a trajectory
		else:
			U_timestep = U.trajectory.dt
			#sanity check
			if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
				print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
				sys.exit(1)
			if U.trajectory.numframes < args.frames_dt:
				print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."
	
			#create list of index of frames to process
			if args.t_end != -1:
				f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
				if f_end < 0:
					print "Error: the ending time specified is before the beginning of the xtc."
					sys.exit(1)
			else:
				f_end = nb_frames_xtc - 1
			if args.t_start != -1:
				f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
				if f_start > f_end:
					print "Error: the starting time specified is after the end of the xtc."
					sys.exit(1)
			if (f_end - f_start)%args.frames_dt == 0:
				tmp_offset = 0
			else:
				tmp_offset = 1
			frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
			nb_frames_to_process = len(frames_to_process)

	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no lipid particles selected, check the --leaflets and/or --beads options."
		sys.exit(1)

	#create selection: particles density
	#-----------------------------------
	water_pres = False
	particles_pres_any = False
	
	#check for presence of each particle
	for part in particles_def["labels"]:
		particles_def["sele"][part] = U.selectAtoms(particles_def["sele_string"][part])
		if particles_def["sele"][part].numberOfAtoms() == 0:
			print " ->warning: particle selection string '" + str(particles_def["sele_string"][part]) + "' returned 0 atoms."
		else:
			particles_pres_any = True
			particles_def_pres[part] = True
			if part == "water":
				water_pres = True
	if not particles_pres_any:
		print "Error: none of the specified particles was found in the system, check your --particles option."
		sys.exit(1)
							
	#create charged particles selections
	#-----------------------------------
	charge_pres_any = False
	if args.chargesfilename != "no":
		for charge_g in charges_groups.keys():
			for q in charges_groups[charge_g]["names"]:
				charges_groups[charge_g]["sele"][q] = U.selectAtoms(charges_groups[charge_g]["sele_string"][q])
				if charges_groups[charge_g]["sele"][q].numberOfAtoms() == 0:
					print " ->warning: charge selection string '" + str(charges_groups[charge_g]["sele_string"][q]) + "' returned 0 atoms."
				else:
					charge_pres_any = True
					charges_groups_pres[charge_g] = True
					charges_groups_pres_q[charge_g][q] = True
		if not charge_pres_any:
			print "Error: no charged particles found, use '--charges no' or supply correct charges definition."
			sys.exit(1)

	#check whether the arg max is less than half the box size
	#--------------------------------------------------------
	if args.slices_range > U.dimensions[2]/float(2):
		print "Warning: the --dist_max option (" + str(args.slices_range) + ") is larger than half the box size in the z direction plane (" + str(U.dimensions[2]/float(2)) + "), density profile will be truncated at " + str(U.dimensions[2]/float(2)) + " Angstrom."
		args.slices_range = U.dimensions[2]/float(2)
		bins_nb_max = int(floor(args.slices_range/float(args.slices_thick)))
		
	#rewind 
	U.trajectory.rewind()
	
	return
def identify_ff():														#DONE
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_bead
	global lipids_sele_ff_bonds
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index = []
	lipids_ff_l2u_index = []
	lipids_sele_ff = {}
	lipids_sele_ff_bead = {}
	lipids_sele_ff_bonds = {}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_sele_string + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		if line[-1] == "\n":
			line = line[:-1]
		try:
			line_content = line.split(',')
			if len(line_content) != 4:
				print "Error: wrong format for line " + str(l_index+1) + " in " + str(args.selection_file_ff) + ", see note 4 in bilayer_perturbations --help."
				print " ->", line
				sys.exit(1)
			#read current lipid details
			lip_resname = line_content[0]
			lip_resnum = int(line_content[1])
			lip_leaflet = line_content[2]
			lip_bead = line_content[3]
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet,lip_bead]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "->unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum))
			lipids_sele_ff_bead[l_index] = lipids_sele_ff[l_index].selectAtoms("name " + str(lip_bead))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + line + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_leaflets():												#DONE
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	global leaflet_nb_beads_upper
	global leaflet_nb_beads_lower
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			if args.use_gro:
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, leaflet_sele_string, cutoff_value[0])
			else:
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			if args.use_gro:
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U_gro, leaflet_sele_string, args.cutoff_leaflet)
			else:
				L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			if args.use_gro:
				tmp_up = L.group(0)
				tmp_lw = L.group(1)
			else:
				leaflet_sele["upper"] = L.group(0)
				leaflet_sele["lower"] = L.group(1)
		else:
			if args.use_gro:
				tmp_up = L.group(1)
				tmp_lw = L.group(0)			
			else:
				leaflet_sele["upper"] = L.group(1)
				leaflet_sele["lower"] = L.group(0)

		if args.use_gro:		
			tmp_up_indices = tmp_up.indices()
			tmp_lw_indices = tmp_lw.indices()
			tmp_up_nb_atoms = len(tmp_up_indices)
			tmp_lw_nb_atoms = len(tmp_lw_indices)
			leaflet_sele["upper"] = U.selectAtoms("bynum " + str(tmp_up_indices[0] + 1))
			leaflet_sele["lower"] = U.selectAtoms("bynum " + str(tmp_lw_indices[0] + 1))
			for index in range(1,tmp_up_nb_atoms):
				leaflet_sele["upper"] += U.selectAtoms("bynum " + str(tmp_up_indices[index] +1 ))
				progress = '\r -identifying upper leaflet from gro file... ' + str(round(index/float(tmp_up_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			print ''
			for index in range(1,tmp_lw_nb_atoms):
				leaflet_sele["lower"] += U.selectAtoms("bynum " + str(tmp_lw_indices[index] + 1))
				progress = '\r -identifying lower leaflet from gro file... ' + str(round(index/float(tmp_lw_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			print ''

		leaflet_sele["both"] = leaflet_sele["lower"] + leaflet_sele["upper"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cog and z coordinates in the GRO file supplied:
	else:
		if args.use_gro:
			tmp_all = U_gro.selectAtoms(leaflet_sele_string)
			tmp_lipids_avg_z = tmp_all.centerOfGeometry()[2]
			tmp_up = tmp_all.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
			tmp_lw = tmp_all.selectAtoms("prop z < " + str(tmp_lipids_avg_z))
			tmp_up_indices = tmp_up.indices()
			tmp_lw_indices = tmp_lw.indices()
			tmp_up_nb_atoms = len(tmp_up_indices)
			tmp_lw_nb_atoms = len(tmp_lw_indices)
			leaflet_sele["upper"] = U.selectAtoms("bynum " + str(tmp_up_indices[0] + 1))
			leaflet_sele["lower"] = U.selectAtoms("bynum " + str(tmp_lw_indices[0] + 1))
			for index in range(1,tmp_up_nb_atoms):
				leaflet_sele["upper"] += U.selectAtoms("bynum " + str(tmp_up_indices[index] + 1))
				progress = '\r -identifying upper leaflet from gro file... ' + str(round(index/float(tmp_up_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			print ''
			for index in range(1,tmp_lw_nb_atoms):
				leaflet_sele["lower"] += U.selectAtoms("bynum " + str(tmp_lw_indices[index] + 1))
				progress = '\r -identifying lower leaflet from gro file... ' + str(round(index/float(tmp_lw_nb_atoms)*100,1)) + '%   '
				sys.stdout.flush()
				sys.stdout.write(progress)
			leaflet_sele["both"] = leaflet_sele["upper"] + leaflet_sele["lower"]
			print ''
		else:
			leaflet_sele["both"] = U.selectAtoms(leaflet_sele_string)
			tmp_lipids_avg_z = leaflet_sele["both"].centerOfGeometry()[2]
			leaflet_sele["upper"] = leaflet_sele["both"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
			leaflet_sele["lower"] = leaflet_sele["both"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))			
		print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'

	#store number of beads used
	leaflet_nb_beads_upper = int(leaflet_sele["upper"].numberOfAtoms())
	leaflet_nb_beads_lower = int(leaflet_sele["lower"].numberOfAtoms())
		
	return

#=========================================================================================
# data structures
#=========================================================================================

def struct_time():														#DONE

	global frames_nb
	global frames_time
	frames_nb = np.zeros(nb_frames_to_process)
	frames_time = np.zeros(nb_frames_to_process)

	return
def struct_data():

	#surface
	#-------
	global z_upper, z_lower
	global nb_voxel_processed
	global grid_statistics_upper_coverage, grid_statistics_lower_coverage
	global grid_statistics_upper_nb_beads, grid_statistics_lower_nb_beads
	z_upper = 0
	z_lower = 0
	nb_voxel_processed = 0		
	grid_statistics_upper_coverage = np.zeros((nb_frames_to_process, leaflet_nb_beads_upper))
	grid_statistics_lower_coverage = np.zeros((nb_frames_to_process, leaflet_nb_beads_lower))
	grid_statistics_upper_nb_beads = np.zeros((nb_frames_to_process, 2))
	grid_statistics_lower_nb_beads = np.zeros((nb_frames_to_process, 2))
	
	#angle
	#-----
	global theta_deg_avg, theta_deg_std
	global norm_dnxdx_avg, norm_dnxdx_std
	global norm_dnydy_avg, norm_dnydy_std
	global norm_dnzdz_avg, norm_dnzdz_std
	theta_deg_avg = np.zeros(nb_frames_to_process)
	theta_deg_std = np.zeros(nb_frames_to_process)
	norm_dnxdx_avg = np.zeros(nb_frames_to_process)
	norm_dnxdx_std = np.zeros(nb_frames_to_process)
	norm_dnydy_avg = np.zeros(nb_frames_to_process)
	norm_dnydy_std = np.zeros(nb_frames_to_process)
	norm_dnzdz_avg = np.zeros(nb_frames_to_process)
	norm_dnzdz_std = np.zeros(nb_frames_to_process)

	#particles
	#---------
	global density_particles_nb
	global density_particles_pc
	density_particles_nb = {part: np.zeros(2*bins_nb) for part in particles_def["labels"]}
	density_particles_pc = {part: np.zeros(2*bins_nb) for part in particles_def["labels"]}

	#charges
	#-------
	if args.chargesfilename != "no":
		global density_charges
		density_charges = {q: np.zeros(2*bins_nb) for q in charges_groups.keys() + ["total"]}

	return
	
#=========================================================================================
# core functions
#=========================================================================================

def coords_remove_whole(coords, box_dim):
	#this function ensures the coordinates are within 0 and box_dim
	#convention: coords between 0 and box_dim in all directions

	coords[:,0] -= np.floor(coords[:,0]/float(box_dim[0])) * box_dim[0]
	coords[:,1] -= np.floor(coords[:,1]/float(box_dim[1])) * box_dim[1]
	
	return coords
def coords_center_in_box(coords, center, box_dim):
	
	#this function centers coords around center and apply pbc along x and y
	#convention: coords between -box_dim/2 and +box_dim/2 along x and y
	
	coords_loc = np.copy(coords)
	
	#centering
	#---------
	coords_loc -= center
	
	#pbc applied along x and y directions
	#------------------------------------
	coords_loc[:,0] -= (np.floor(2*coords_loc[:,0]/float(box_dim[0])) + (1-np.sign(coords_loc[:,0]))/2) * box_dim[0]
	coords_loc[:,1] -= (np.floor(2*coords_loc[:,1]/float(box_dim[1])) + (1-np.sign(coords_loc[:,1]))/2) * box_dim[1]
	
	return coords_loc
def detect_clusters_connectivity(dist, box_dim):						#DONE
	
	#use networkx algorithm
	connected=(dist<args.cutoff_connect)
	network=nx.Graph(connected)
	groups=nx.connected_components(network)
	
	return groups
def detect_clusters_density(dist, box_dim):								#DONE
	
	#run DBSCAN algorithm
	dbscan_output = DBSCAN(eps = args.dbscan_dist, metric = 'precomputed', min_samples = args.dbscan_nb).fit(dist)

	#build 'groups' structure i.e. a list whose element are all the clusters identified
	groups = []
	for c_lab in np.unique(dbscan_output.labels_):
		tmp_pos = np.argwhere(dbscan_output.labels_ == c_lab)
		if c_lab == -1:
			groups += map(lambda p:p[0] , tmp_pos)
		else:
			groups.append(map(lambda p:p[0] , tmp_pos))

	return groups
def calculate_cog(tmp_coords, box_dim):									#DONE
	
	#this method allows to take pbc into account when calculcating the center of geometry 
	#see: http://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
	
	cog_coord = np.zeros(3)
	tmp_nb_atoms = np.shape(tmp_coords)[0]
	
	for n in range(0,3):
		tet = tmp_coords[:,n] * 2 * math.pi / float(box_dim[n])
		xsi = np.cos(tet)
		zet = np.sin(tet)
		tet_avg = math.atan2(-np.average(zet),-np.average(xsi)) + math.pi
		cog_coord[n] = tet_avg * box_dim[n] / float(2*math.pi)
	
	return cog_coord
def calculate_properties(box_dim, f_nb):								#DONE
	
	global z_upper, z_lower
	global nb_voxel_processed
	global grid_statistics_upper_coverage, grid_statistics_lower_coverage
	global grid_statistics_upper_nb_beads, grid_statistics_lower_nb_beads
	theta_tot_tmp = []
	loc_z_axis = np.array([0,0,1])
	loc_z_axis = loc_z_axis.reshape((3,1))
	tmp_grid_statistics_upper_nb_beads = []
	tmp_grid_statistics_lower_nb_beads = []
	tmp_pos = np.zeros(((box_dim[0] // args.voxel_x) + 1 , (box_dim[1] // args.voxel_y) + 1, (box_dim[2] // args.voxel_z) + 1, 3))
	tmp_norm = np.zeros(((box_dim[0] // args.voxel_x) + 1 , (box_dim[1] // args.voxel_y) + 1, (box_dim[2] // args.voxel_z) + 1, 3))
	tmp_norm_d = np.zeros(((box_dim[0] // args.voxel_x) + 1 , (box_dim[1] // args.voxel_y) + 1, (box_dim[2] // args.voxel_z) + 1, 3), float)
	
	#create coordinate array of "middle" leaflet
	#===========================================
	#retrieve leaflets coords
	tmp_lip_coords = {l: coords_remove_whole(leaflet_sele[l].coordinates(), box_dim) for l in ["lower","upper"]}
	tmp_upper = np.copy(tmp_lip_coords["upper"])
	tmp_lower = np.copy(tmp_lip_coords["lower"])

	#calculate middle of bilayer and relative coordinate of upper and lower leaflets assuming the z is the normal to the bilayer
	cog_up = np.average(tmp_upper, axis = 0)
	cog_lw = np.average(tmp_lower, axis = 0)
	cog_delta = (cog_up - cog_lw) / float(2)
	tmp_upper -= cog_delta
	tmp_lower += cog_delta
	tmp_both = np.concatenate((np.float32(tmp_upper), np.float32(tmp_lower)))
	
	#downsample using voxel filtering
	#================================
	p = pcl.PointCloud()
	p.from_array(tmp_both)
	voxel_grid = p.make_voxel_grid_filter()
	voxel_grid.set_leaf_size(args.voxel_x, args.voxel_y, args.voxel_z)
	voxel_grid.set_minimum_points_number_per_voxel(args.voxel_nb)
	voxel_coord = voxel_grid.filter().to_array()

	#cache coordinates of particles for which to calculate TM density
	#================================================================
	#particles
	tmp_coord_p = {}
	for part in particles_def["labels"]:
		if particles_def_pres[part]:
			tmp_coord_p[part] = coords_remove_whole(particles_def["sele"][part].coordinates(), box_dim)
	#charges
	if args.chargesfilename != "no":
		tmp_coord_q = {}
		for charge_g in charges_groups.keys():
			tmp_coord_q[charge_g] = {}
			if charges_groups_pres[charge_g]:
				for q in charges_groups[charge_g]["names"]:
					if charges_groups_pres_q[charge_g][q]:
						tmp_coord_q[charge_g][q] = coords_remove_whole(charges_groups[charge_g]["sele"][q].coordinates(), box_dim)

	#process each occupied voxel
	#===========================
	v_counter = 0
	v_nb = np.shape(voxel_coord)[0]
	for v_index in range(0, v_nb):
		#display update
		nb_voxel_processed += 1
		v_counter += 1
		progress = '\r -processing frame ' + str(f_nb + 1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' from ' + str(f_start) + ' to ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ') and voxel ' + str(v_counter) + '/' + str(v_nb) + '              '
		sys.stdout.flush()
		sys.stdout.write(progress)
	
		#retrieve coord of point in voxel
		tmp_voxel_center = voxel_coord[v_index,:]
		
		#calculate local normal to bilayer
		#---------------------------------
		if args.normal != 'z':
			#switch to voxel center referential
			tmp_lip_coords_up_centered = coords_center_in_box(tmp_lip_coords["upper"], tmp_voxel_center, box_dim)
			tmp_lip_coords_lw_centered = coords_center_in_box(tmp_lip_coords["lower"], tmp_voxel_center, box_dim)
											
			#identify neighbouring particles in each leaflet
			selected_upper = tmp_lip_coords_up_centered[:,0]**2 + tmp_lip_coords_up_centered[:,1]**2 + tmp_lip_coords_up_centered[:,2]**2 < args.normal_d**2
			selected_lower = tmp_lip_coords_lw_centered[:,0]**2 + tmp_lip_coords_lw_centered[:,1]**2 + tmp_lip_coords_lw_centered[:,2]**2 < args.normal_d**2
			grid_statistics_upper_coverage[f_nb, selected_upper] += 1
			grid_statistics_lower_coverage[f_nb, selected_lower] += 1
			
			#keep their coords only
			tmp_lip_coords_up_centered_within = tmp_lip_coords_up_centered[selected_upper]
			tmp_lip_coords_lw_centered_within = tmp_lip_coords_lw_centered[selected_lower]
			tmp_grid_statistics_upper_nb_beads.append(int(np.shape(tmp_lip_coords_up_centered_within)[0]))
			tmp_grid_statistics_lower_nb_beads.append(int(np.shape(tmp_lip_coords_lw_centered_within)[0]))
			if tmp_grid_statistics_upper_nb_beads[-1] == 0:
				print "\nWarning: no neighbouring particles found in the upper leaflet for current voxel. Check the normal and voxel options.\n"
				continue
			else:
				cog_up = np.average(tmp_lip_coords_up_centered_within, axis = 0)
			if tmp_grid_statistics_lower_nb_beads[-1] == 0:
				print "\nWarning: no neighbouring particles found in the lower leaflet for current voxel. Check the normal and voxel options.\n"
				continue
			else:
				cog_lw = np.average(tmp_lip_coords_lw_centered_within, axis = 0)
						
			#identify normal vector: case cog
			if args.normal == 'cog':
				norm_vec = cog_up - cog_lw
				norm_vec /= float(np.linalg.norm(norm_vec))
				norm_vec = norm_vec.reshape((3,1))
			#identify normal vector: case svd
			else:
				tmp_lip_coords_centered_within = np.concatenate((coords_center_in_box(tmp_lip_coords_up_centered_within, cog_up, box_dim), coords_center_in_box(tmp_lip_coords_lw_centered_within, cog_lw, box_dim)))
				svd_U, svd_D, svd_V = np.linalg.svd(tmp_lip_coords_centered_within)
				norm_vec = svd_V[2].reshape((3,1))

			#orientate the normal vector so that from inner leaflet towards outer one
			tmp_cog_delta = cog_up - cog_lw
			tmp_cog_delta = tmp_cog_delta.reshape((3,1))
			if np.dot(tmp_cog_delta[:,0],norm_vec[:,0]) < 0:
				norm_vec *= -1

			#identify rotation matrix
			#NB: see http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
			norm_ax = np.cross(loc_z_axis,norm_vec,axis=0)				# z axis ^ local normal = axis of rotation
			norm_cos = np.dot(loc_z_axis[:,0],norm_vec[:,0])			#if the cos is <0 the angle will be >90deg, this is technically possible (e.g. buckling) but should only happen in extreme deformation cases
			norm_sin = np.linalg.norm(norm_ax)							#one way to get around this is to take the abs of the cos (i.e. of np.dot) but that means we assume that the angle has to be <90 deg
			norm_ax_skew_sym = norm_vec*loc_z_axis.T - loc_z_axis*norm_vec.T
			norm_rot = np.identity(3) - norm_ax_skew_sym + (1-norm_cos)/float(norm_sin**2)*np.dot(norm_ax_skew_sym,norm_ax_skew_sym)
		
			#get rotation angle
			theta_tot_tmp.append(np.arctan2(norm_sin, norm_cos) * 180 /float(np.pi))
			
			#store it in the voxel grid
			tmp_index = (tmp_voxel_center[0] // args.voxel_x, tmp_voxel_center[1] // args.voxel_y, tmp_voxel_center[2] // args.voxel_z)
			if np.any(tmp_norm[(tmp_index)]):
				print "Warning: voxel sampled twice in the same frame"
			else:
				tmp_pos[(tmp_index)] = tmp_voxel_center
				tmp_norm[(tmp_index)] = norm_vec[:,0]
								
			#ROTATION
			#rotate neighbouring bilayer in local cluster referential
			tmp_lip_coords_up_centered_within_rotated = np.dot(norm_rot, tmp_lip_coords_up_centered_within.T).T
			tmp_lip_coords_lw_centered_within_rotated = np.dot(norm_rot, tmp_lip_coords_lw_centered_within.T).T
			
			#identify z coord of local middle of bilayer after rotation
			cog_up_rotated_z = np.median(tmp_lip_coords_up_centered_within_rotated[:,2])
			cog_lw_rotated_z = np.median(tmp_lip_coords_lw_centered_within_rotated[:,2])
			norm_z_middle = cog_lw_rotated_z + (cog_up_rotated_z - cog_lw_rotated_z)/float(2)
											
			#TRANSLATION
			tmp_lip_coords_up_centered_within_rotated[:,2] -= norm_z_middle
			tmp_lip_coords_lw_centered_within_rotated[:,2] -= norm_z_middle
			#store relative coordinate of local upper and lower leaflets (once they've been rotated in the x,y plane)
			z_upper += cog_up_rotated_z - norm_z_middle
			z_lower += cog_lw_rotated_z - norm_z_middle
				
		else:
			z_upper += cog_delta[2]
			z_lower -= cog_delta[2]
			norm_z_middle = cog_lw[2] + cog_delta[2]
							
		#density profile: particles
		#--------------------------
		for part in particles_def["labels"]:
			if particles_def_pres[part]:
				#retrieve coords
				tmp_coord = np.copy(tmp_coord_p[part])
								
				#performs centering/rotating of the referential
				if args.normal != 'z':
					#switch to voxel center referential
					tmp_coord = coords_center_in_box(tmp_coord, tmp_voxel_center, box_dim)
													
					#rotate coordinates so that the local normal of the bilayer is // to the z axis
					tmp_coord = np.dot(norm_rot, tmp_coord.T).T
				
					#center around middle of rotated bilayer in z
					tmp_coord[:,2] -= norm_z_middle
				else:					
					#center around cluster in the x and y direction
					tmp_coord[:,0] -= tmp_voxel_center[0]
					tmp_coord[:,1] -= tmp_voxel_center[1]
									
					#center z coordinates on the bilayer center z coordinate
					tmp_coord[:,2] -= norm_z_middle
			
				#keep those within the specified radius
				tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
								
				#add number of particles within each slice					
				tmp_bins_nb = np.zeros(2*bins_nb)
				bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
				bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)			#the + int(bins_nb) allows to only have positive bin indices
				if len(bin_abs) > 0:				
					tmp_bins_nb = np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0]
				density_particles_nb[part] += tmp_bins_nb
						
		#density profile: charges
		#------------------------
		if args.chargesfilename != "no":
			for charge_g in charges_groups.keys():
				if charges_groups_pres[charge_g]:
					tmp_bins_nb = np.zeros(2*bins_nb)
					for q in charges_groups[charge_g]["names"]:
						if charges_groups_pres_q[charge_g][q]:
							#retrieve coords
							tmp_coord = np.copy(tmp_coord_q[charge_g][q])
							
							#performs centering/rotating of the referential
							if args.normal != 'z':
								#switch to voxel center referential
								tmp_coord = coords_center_in_box(tmp_coord, tmp_voxel_center, box_dim)
								
								#rotate coordinates so that the local normal of the bilayer is // to the z axis
								tmp_coord = np.dot(norm_rot, tmp_coord.T).T
								
								#center around middle of rotated bilayer in z
								tmp_coord[:,2] -= norm_z_middle
							else:					
								#center around cluster in the x and y direction
								tmp_coord[:,0] -= tmp_voxel_center[0]
								tmp_coord[:,1] -= tmp_voxel_center[1]
												
								#center z coordinates on the bilayer center z coordinate
								tmp_coord[:,2] -= norm_z_middle
							
							#keep those within the specified radius
							tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
				
							#add number of particles within each slice
							bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
							bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)
							if len(bin_abs) > 0:
								tmp_bins_nb += np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0] * charges_groups[charge_g]["values"][q]
					density_charges[charge_g] += tmp_bins_nb 
					density_charges["total"] += tmp_bins_nb

		#calculate some grid coverage statistics
		#---------------------------------------
		if args.normal != 'z':
			grid_statistics_upper_nb_beads[f_nb, 0] = np.average(tmp_grid_statistics_upper_nb_beads)
			grid_statistics_upper_nb_beads[f_nb, 1] = np.std(tmp_grid_statistics_upper_nb_beads)
			grid_statistics_lower_nb_beads[f_nb, 0] = np.average(tmp_grid_statistics_lower_nb_beads)
			grid_statistics_lower_nb_beads[f_nb, 1] = np.std(tmp_grid_statistics_lower_nb_beads)

	#calculate convexity
	#===================	
	#create grids structures
	#-----------------------
	grid_shape = np.shape(tmp_norm)
	#shifted grids: local normal
	tmp_shifted_norm_x = np.zeros(grid_shape)
	tmp_shifted_norm_y = np.zeros(grid_shape)
	tmp_shifted_norm_z = np.zeros(grid_shape)
	tmp_shifted_norm_xy = np.zeros(grid_shape)
	tmp_shifted_norm_xz = np.zeros(grid_shape)
	tmp_shifted_norm_yz = np.zeros(grid_shape)
	tmp_shifted_norm_xyz = np.zeros(grid_shape)	
	#shifted grids: positions
	tmp_shifted_pos_x = np.zeros(grid_shape)
	tmp_shifted_pos_y = np.zeros(grid_shape)
	tmp_shifted_pos_z = np.zeros(grid_shape)
	tmp_shifted_pos_xy = np.zeros(grid_shape)
	tmp_shifted_pos_xz = np.zeros(grid_shape)
	tmp_shifted_pos_yz = np.zeros(grid_shape)
	tmp_shifted_pos_xyz = np.zeros(grid_shape)		
	#deltas between shifted grids: local normal
	tmp_delta_norm_x = np.zeros(grid_shape)
	tmp_delta_norm_y = np.zeros(grid_shape)
	tmp_delta_norm_z = np.zeros(grid_shape)
	tmp_delta_norm_xy = np.zeros(grid_shape)
	tmp_delta_norm_xz = np.zeros(grid_shape)
	tmp_delta_norm_yz = np.zeros(grid_shape)
	tmp_delta_norm_xyz = np.zeros(grid_shape)	
	#deltas between shifted grids: positions
	tmp_delta_pos_x = np.zeros(grid_shape)
	tmp_delta_pos_y = np.zeros(grid_shape)
	tmp_delta_pos_z = np.zeros(grid_shape)
	tmp_delta_pos_xy = np.zeros(grid_shape)
	tmp_delta_pos_xz = np.zeros(grid_shape)
	tmp_delta_pos_yz = np.zeros(grid_shape)
	tmp_delta_pos_xyz = np.zeros(grid_shape)	

	#fill shifted grids: local normal
	#--------------------------------
	#primary shifts
	tmp_shifted_norm_x[1:,:] = tmp_norm[:-1,:]							#shift along x
	tmp_shifted_norm_y[:,1:,:] = tmp_norm[:,:-1,:]						#shift along y
	tmp_shifted_norm_z[:,:,1:,:] = tmp_norm[:,:,:-1,:]					#shift along z
	tmp_shifted_norm_x[0,:] = tmp_norm[-1,:] 							#pbc: x
	tmp_shifted_norm_y[:,0,:] = tmp_norm[:,-1,:] 						#pbc: y
	tmp_shifted_norm_y[:,:,0,:] = tmp_norm[:,:,-1,:] 					#pbc: z
	#secondary shifts
	tmp_shifted_norm_xy[:,1:,:] = tmp_shifted_norm_x[:,:-1,:]			#shift along y
	tmp_shifted_norm_xz[:,:,1:,:] = tmp_shifted_norm_x[:,:,:-1,:]		#shift along z
	tmp_shifted_norm_yz[:,:,1:,:] = tmp_shifted_norm_y[:,:,:-1,:]		#shift along z
	tmp_shifted_norm_xy[:,0,:] = tmp_shifted_norm_x[:,-1,:]				#pbc: y
	tmp_shifted_norm_xz[:,:,0,:] = tmp_shifted_norm_x[:,:,-1,:]			#pbc: z
	tmp_shifted_norm_yz[:,:,0,:] = tmp_shifted_norm_y[:,:,-1,:]			#pbc: z
	#tertiary shift
	tmp_shifted_norm_xyz[:,:,1:,:] = tmp_shifted_norm_xy[:,:,:-1,:]		#shift along z
	tmp_shifted_norm_xyz[:,:,0,:] = tmp_shifted_norm_xy[:,:,-1,:]		#pbc: z
	
	#fill shifted grids: positions
	#-----------------------------
	#primary shifts
	tmp_shifted_pos_x[1:,:] = tmp_pos[:-1,:]									#shift along x
	tmp_shifted_pos_y[:,1:,:] = tmp_pos[:,:-1,:]								#shift along y
	tmp_shifted_pos_z[:,:,1:,:] = tmp_pos[:,:,:-1,:]							#shift along z
	tmp_shifted_pos_x[0,:] = tmp_pos[-1,:] - box_dim[0]							#pbc: x
	tmp_shifted_pos_y[:,0,:] = tmp_pos[:,-1,:] - box_dim[1]						#pbc: y
	tmp_shifted_pos_y[:,:,0,:] = tmp_pos[:,:,-1,:] - box_dim[2]					#pbc: z
	#secondary shifts
	tmp_shifted_pos_xy[:,1:,:] = tmp_shifted_pos_x[:,:-1,:]						#shift along y
	tmp_shifted_pos_xz[:,:,1:,:] = tmp_shifted_pos_x[:,:,:-1,:]					#shift along z
	tmp_shifted_pos_yz[:,:,1:,:] = tmp_shifted_pos_y[:,:,:-1,:]					#shift along z
	tmp_shifted_pos_xy[:,0,:] = tmp_shifted_pos_x[:,-1,:] - box_dim[1]			#pbc: y
	tmp_shifted_pos_xz[:,:,0,:] = tmp_shifted_pos_x[:,:,-1,:] - box_dim[2]		#pbc: z
	tmp_shifted_pos_yz[:,:,0,:] = tmp_shifted_pos_y[:,:,-1,:] - box_dim[2]		#pbc: z
	#tertiary shift
	tmp_shifted_pos_xyz[:,:,1:,:] = tmp_shifted_pos_xy[:,:,:-1,:]				#shift along z
	tmp_shifted_pos_xyz[:,:,0,:] = tmp_shifted_pos_xy[:,:,-1,:]	- box_dim[2]	#pbc: z

	#fill deltas grids: local normal
	#-------------------------------
	tmp_delta_norm_x = tmp_norm - tmp_shifted_norm_x
	tmp_delta_norm_y = tmp_norm - tmp_shifted_norm_y
	tmp_delta_norm_z = tmp_norm - tmp_shifted_norm_z
	tmp_delta_norm_xy = tmp_norm - tmp_shifted_norm_xy
	tmp_delta_norm_xz = tmp_norm - tmp_shifted_norm_xz
	tmp_delta_norm_yz = tmp_norm - tmp_shifted_norm_yz
	tmp_delta_norm_xyz = tmp_norm - tmp_shifted_norm_xyz

	#fill deltas grids: positions
	#----------------------------
	tmp_delta_pos_x = tmp_pos - tmp_shifted_pos_x
	tmp_delta_pos_y = tmp_pos - tmp_shifted_pos_y
	tmp_delta_pos_z = tmp_pos - tmp_shifted_pos_z
	tmp_delta_pos_xy = tmp_pos - tmp_shifted_pos_xy
	tmp_delta_pos_xz = tmp_pos - tmp_shifted_pos_xz
	tmp_delta_pos_yz = tmp_pos - tmp_shifted_pos_yz
	tmp_delta_pos_xyz = tmp_pos - tmp_shifted_pos_xyz

	#find out which voxels are populated
	#-----------------------------------
	tmp_pop_vox = np.any(tmp_norm[:,:,:] != 0, 3)
	tmp_pop_x = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_x[:,:,:] != 0, 3))
	tmp_pop_y = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_y[:,:,:] != 0, 3))
	tmp_pop_z = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_z[:,:,:] != 0, 3))
	tmp_pop_xy = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_xy[:,:,:] != 0, 3))
	tmp_pop_xz = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_xz[:,:,:] != 0, 3))	
	tmp_pop_yz = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_yz[:,:,:] != 0, 3))
	tmp_pop_xyz = np.multiply(tmp_pop_vox, np.any(tmp_shifted_norm_xyz[:,:,:] != 0, 3))
	
	#calculate derivatives
	#---------------------
	
	#----------------------------------------------------------------------------------------------
	#NB: we could calculate a lot more (9 in total - each component of the norm in each direction)
	#    but we only really care about d(nx)/dx, d(ny)/dy and (d(nz)/dz (that last one is particularly telling)
	#----------------------------------------------------------------------------------------------
	
	#d(nx)/dx ([0] = nx and x, xy, xz and xyz contribute)
	if np.count_nonzero(tmp_pop_x) > 0:
		tmp_norm_d[tmp_pop_x, 0] += tmp_delta_norm_x[tmp_pop_x, 0] / tmp_delta_pos_x[tmp_pop_x, 0]
	if np.count_nonzero(tmp_pop_xy) > 0:
		tmp_norm_d[tmp_pop_xy, 0] += tmp_delta_norm_xy[tmp_pop_xy, 0] / tmp_delta_pos_xy[tmp_pop_xy, 0]
	if np.count_nonzero(tmp_pop_xz) > 0:
		tmp_norm_d[tmp_pop_xz, 0] += tmp_delta_norm_xz[tmp_pop_xz, 0] / tmp_delta_pos_xz[tmp_pop_xz, 0]
	if np.count_nonzero(tmp_pop_xyz) > 0:
		tmp_norm_d[tmp_pop_xyz, 0] += tmp_delta_norm_xyz[tmp_pop_xyz, 0] / tmp_delta_pos_xyz[tmp_pop_xyz, 0]
	#d(ny)/dy ([1] = ny and y, xy, yz and xyz contribute)
	if np.count_nonzero(tmp_pop_y) > 0:
		tmp_norm_d[tmp_pop_y, 1] += tmp_delta_norm_y[tmp_pop_y, 1] / tmp_delta_pos_y[tmp_pop_y, 1]
	if np.count_nonzero(tmp_pop_xy) > 0:
		tmp_norm_d[tmp_pop_xy, 1] += tmp_delta_norm_xy[tmp_pop_xy, 1] / tmp_delta_pos_xy[tmp_pop_xy, 1]
	if np.count_nonzero(tmp_pop_yz) > 0:
		tmp_norm_d[tmp_pop_yz, 1] += tmp_delta_norm_yz[tmp_pop_yz, 1] / tmp_delta_pos_yz[tmp_pop_yz, 1]
	if np.count_nonzero(tmp_pop_xyz) > 0:
		tmp_norm_d[tmp_pop_xyz, 1] += tmp_delta_norm_xyz[tmp_pop_xyz, 1] / tmp_delta_pos_xyz[tmp_pop_xyz, 1]
	#d(nz)/dz ([2] = nz and z, xz, yz and xyz contribute)
	if np.count_nonzero(tmp_pop_z) > 0:
		tmp_norm_d[tmp_pop_z, 2] += tmp_delta_norm_z[tmp_pop_z, 2] / tmp_delta_pos_z[tmp_pop_z, 2]
	if np.count_nonzero(tmp_pop_xz) > 0:
		tmp_norm_d[tmp_pop_xz, 2] += tmp_delta_norm_xz[tmp_pop_xz, 2] / tmp_delta_pos_xz[tmp_pop_xz, 2]
	if np.count_nonzero(tmp_pop_yz) > 0:
		tmp_norm_d[tmp_pop_yz, 2] += tmp_delta_norm_yz[tmp_pop_yz, 2] / tmp_delta_pos_yz[tmp_pop_yz, 2]
	if np.count_nonzero(tmp_pop_xyz) > 0:
		tmp_norm_d[tmp_pop_xyz, 2] += tmp_delta_norm_xyz[tmp_pop_xyz, 2] / tmp_delta_pos_xyz[tmp_pop_xyz, 2]

	#calculate avg derivatives
	if np.count_nonzero(tmp_pop_x + tmp_pop_xy + tmp_pop_xz + tmp_pop_xyz) > 0:
		norm_dnxdx_avg[f_nb] = np.average(tmp_norm_d[tmp_pop_x + tmp_pop_xy + tmp_pop_xz + tmp_pop_xyz, 0])
		norm_dnxdx_std[f_nb] = np.std(tmp_norm_d[tmp_pop_x + tmp_pop_xy + tmp_pop_xz + tmp_pop_xyz, 0])
	if np.count_nonzero(tmp_pop_y + tmp_pop_xy + tmp_pop_yz + tmp_pop_xyz) > 0:
		norm_dnydy_avg[f_nb] = np.average(tmp_norm_d[tmp_pop_y + tmp_pop_xy + tmp_pop_yz + tmp_pop_xyz, 1])
		norm_dnydy_std[f_nb] = np.std(tmp_norm_d[tmp_pop_y + tmp_pop_xy + tmp_pop_yz + tmp_pop_xyz, 1])
	if np.count_nonzero(tmp_pop_z + tmp_pop_xz + tmp_pop_yz + tmp_pop_xyz) > 0:
		norm_dnzdz_avg[f_nb] = np.average(tmp_norm_d[tmp_pop_z + tmp_pop_yz + tmp_pop_xz + tmp_pop_xyz, 2])
		norm_dnzdz_std[f_nb] = np.std(tmp_norm_d[tmp_pop_z + tmp_pop_yz + tmp_pop_xz + tmp_pop_xyz, 2])

	#calculate angle average
	#=======================
	if len(theta_tot_tmp) != 0:
		theta_deg_avg[f_nb] = np.average(theta_tot_tmp)
		theta_deg_std[f_nb] = np.std(theta_tot_tmp)
	else:
		theta_deg_avg[f_nb] = np.nan
		theta_deg_std[f_nb] = np.nan

	return
def calculate_stats():													#DONE
	
	global z_upper, z_lower
	global max_density_particles_pc
	max_density_particles_pc  = float("-inf")
	if args.chargesfilename != "no":
		global max_density_charges
		global min_density_charges
		max_density_charges = float("-inf")
		min_density_charges = float("inf")
		
	#leafelts coords
	#---------------
	z_upper /= float(nb_voxel_processed)
	z_lower /= float(nb_voxel_processed)
	
	#calculate grid coverage statistics
	#----------------------------------
	if args.normal != 'z':
		#average nb of times beads are sampled per frame
		global avg_sampling_upper_avg, avg_sampling_upper_std
		global avg_sampling_lower_avg, avg_sampling_lower_std
		avg_sampling_upper_each_frame = np.average(grid_statistics_upper_coverage, axis = 1)
		avg_sampling_lower_each_frame = np.average(grid_statistics_lower_coverage, axis = 1)
		avg_sampling_upper_avg = np.average(avg_sampling_upper_each_frame)
		avg_sampling_upper_std = np.std(avg_sampling_upper_each_frame)
		avg_sampling_lower_avg = np.average(avg_sampling_lower_each_frame)
		avg_sampling_lower_std = np.std(avg_sampling_lower_each_frame)
		#average % of beads sampled 0 times per frame
		global upper_pc_zeros_avg, upper_pc_zeros_std
		global lower_pc_zeros_avg, lower_pc_zeros_std
		upper_pc_zeros_each_frame = np.sum(grid_statistics_upper_coverage == 0, axis = 1) / float(leaflet_nb_beads_upper) * 100
		lower_pc_zeros_each_frame = np.sum(grid_statistics_lower_coverage == 0, axis = 1) / float(leaflet_nb_beads_lower) * 100
		upper_pc_zeros_avg = np.average(upper_pc_zeros_each_frame)
		upper_pc_zeros_std = np.std(upper_pc_zeros_each_frame)
		lower_pc_zeros_avg = np.average(lower_pc_zeros_each_frame)
		lower_pc_zeros_std = np.std(lower_pc_zeros_each_frame)
		#average % of beads sampled more than 1 time per frame
		global upper_pc_mult_avg, upper_pc_mult_std
		global lower_pc_mult_avg, lower_pc_mult_std
		upper_pc_mult_each_frame = np.sum(grid_statistics_upper_coverage > 1, axis = 1) / float(leaflet_nb_beads_upper) * 100
		lower_pc_mult_each_frame = np.sum(grid_statistics_lower_coverage > 1, axis = 1) / float(leaflet_nb_beads_lower) * 100
		upper_pc_mult_avg = np.average(upper_pc_mult_each_frame)
		upper_pc_mult_std = np.std(upper_pc_mult_each_frame)
		lower_pc_mult_avg = np.average(lower_pc_mult_each_frame)
		lower_pc_mult_std = np.std(lower_pc_mult_each_frame)
	
		#average nb of upper beads found around each voxel center
		global grid_statistics_upper_nb_beads_avg_avg, grid_statistics_upper_nb_beads_avg_std
		global grid_statistics_upper_nb_beads_std_avg, grid_statistics_upper_nb_beads_std_std
		global grid_statistics_lower_nb_beads_avg_avg, grid_statistics_lower_nb_beads_avg_std
		global grid_statistics_lower_nb_beads_std_avg, grid_statistics_lower_nb_beads_std_std
		grid_statistics_upper_nb_beads_avg_avg = np.average(grid_statistics_upper_nb_beads[:,0])
		grid_statistics_upper_nb_beads_avg_std = np.std(grid_statistics_upper_nb_beads[:,0])
		grid_statistics_upper_nb_beads_std_avg = np.average(grid_statistics_upper_nb_beads[:,1])
		grid_statistics_upper_nb_beads_std_std = np.std(grid_statistics_upper_nb_beads[:,1])
		grid_statistics_lower_nb_beads_avg_avg = np.average(grid_statistics_lower_nb_beads[:,0])
		grid_statistics_lower_nb_beads_avg_std = np.std(grid_statistics_lower_nb_beads[:,0])
		grid_statistics_lower_nb_beads_std_avg = np.average(grid_statistics_lower_nb_beads[:,1])
		grid_statistics_lower_nb_beads_std_std = np.std(grid_statistics_lower_nb_beads[:,1])
	
	#calculate normalisation factor for each group for each size
	#------------------------------------------------------------
	tmp_normalisation = {}
	for part_g in particles_groups.keys():
		tmp_normalisation[part_g] = 0
		for part in particles_groups[part_g]:
			try:
				tmp_normalisation[part_g] += np.sum(density_particles_nb[part]) * args.slices_thick
			except:
				pass		
	
	#density profile: particles
	#--------------------------
	for part in particles_def["labels"]:
		if particles_def_pres[part]:
			#relative density
			if tmp_normalisation[particles_def["group"][part]] > 0:
				density_particles_pc[part] = density_particles_nb[part] / float(tmp_normalisation[particles_def["group"][part]])
			
			#update scale
			if part != "Na+" and part != "Cl-":
				max_density_particles_pc = max(max_density_particles_pc, max(density_particles_pc[part]))
	
	#density profile: charges
	#------------------------
	if args.chargesfilename != "no":
		tmp_normalisation = nb_voxel_processed * slice_volume
		for charge_g in charges_groups.keys():
			if charges_groups_pres[charge_g]:
				#absolute density
				density_charges[charge_g] /= float(tmp_normalisation)
				
				#update scale
				max_density_charges = max(max_density_charges, max(density_charges[charge_g]))
				min_density_charges = min(min_density_charges, min(density_charges[charge_g]))
		
		#total charge: absolute density
		density_charges["total"] /= float(tmp_normalisation)
		
		#total charge: update scale
		max_density_charges = max(max_density_charges, max(density_charges["total"]))
		min_density_charges = min(min_density_charges, min(density_charges["total"]))

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_grid_statistics():
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/grid_statistics.txt'
	output_txt = open(filename_txt, 'w')
	
	#general header
	output_txt.write("[bilayer sampling statistics - written by membrane_prop v" + str(version_nb) + "]\n")
	
	#description
	output_txt.write("\n")
	output_txt.write("This file summarises sampling statistics of the lipid bilayer given the voxel grid parameters the lipids beads chosen.\n")
	output_txt.write("\n")
	output_txt.write("- lipid beads selection: " + str(leaflet_sele_string_original) + "\n")
	output_txt.write("- number beads in upper: " + str(leaflet_nb_beads_upper) + "\n")
	output_txt.write("- number beads in lower: " + str(leaflet_nb_beads_lower) + "\n")
	output_txt.write("- voxel x,y,z dimensions: " + str(args.voxel_x) + ", " + str(args.voxel_y) + ", " + str(args.voxel_z) + "\n")
	output_txt.write("- voxel min nb of points: " + str(args.voxel_nb) + "\n")
	output_txt.write("- search radius from voxel center to select beads used for local normal identification: " + str(args.normal_d) + " Angstrom\n")

	#data
	output_txt.write("\n")
	output_txt.write("Average number of times each bead is sampled per frame:\n")
	output_txt.write("-upper: " + str(round(avg_sampling_upper_avg,2)) + " (" + str(round(avg_sampling_upper_std,2)) + ")\n")
	output_txt.write("-lower: " + str(round(avg_sampling_lower_avg,2)) + " (" + str(round(avg_sampling_lower_std,2)) + ")\n")
	output_txt.write("\n")
	output_txt.write("Average % of beads not sampled at all per frame:\n")
	output_txt.write("-upper: " + str(round(upper_pc_zeros_avg,2)) + " (" + str(round(upper_pc_zeros_std,2)) + ")\n")
	output_txt.write("-lower: " + str(round(lower_pc_zeros_avg,2)) + " (" + str(round(lower_pc_zeros_std,2)) + ")\n")
	output_txt.write("\n")
	output_txt.write("Average % of beads sampled more than once per frame:\n")
	output_txt.write("-upper: " + str(round(upper_pc_mult_avg,2)) + " (" + str(round(upper_pc_mult_std,2)) + ")\n")
	output_txt.write("-lower: " + str(round(lower_pc_mult_avg,2)) + " (" + str(round(lower_pc_mult_std,2)) + ")\n")
	output_txt.write("\n")
	output_txt.write("Average number of beads used for normal identification:\n")
	output_txt.write("-upper: " + str(round(grid_statistics_upper_nb_beads_avg_avg,2)) + " (" + str(round(grid_statistics_upper_nb_beads_avg_std,2)) + ") +/-" + str(round(grid_statistics_upper_nb_beads_std_avg,2)) + " (" + str(round(grid_statistics_upper_nb_beads_std_std,2)) + ")\n")
	output_txt.write("-lower: " + str(round(grid_statistics_lower_nb_beads_avg_avg,2)) + " (" + str(round(grid_statistics_lower_nb_beads_avg_std,2)) + ") +/-" + str(round(grid_statistics_lower_nb_beads_std_avg,2)) + " (" + str(round(grid_statistics_lower_nb_beads_std_std,2)) + ")\n")
	output_txt.close()

	return

def density_write_particles():											#DONE

	tmp_file = 'membrane_prop_TM_density_particles'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/density/' + str(tmp_file) + '.xvg'
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/density/' + str(tmp_file) + '.txt'
	output_txt = open(filename_txt, 'w')
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_txt.write("@[relative particles frequency profile - written by membrane_prop v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_file) + ".xvg.\n")
	output_xvg.write("# [relative particles frequency profile - written by membrane_prop v" + str(version_nb) + "]\n")
 	output_xvg.write("# volume properties:\n")
	output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
	output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
	output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Relative particles frequency along z\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"relative particles frequency (Angstrom-3)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(particles_def["labels"])) + "\n")
	for part_index in range(0,len(particles_def["labels"])):
		part = particles_def["labels"][part_index]
		output_xvg.write("@ s" + str(part_index) + " legend \"" + str(part) + "\"\n")
		output_txt.write(str(tmp_file) + "," + str(part_index+1) + "," + str(part) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"][part])) + "\n")
	output_txt.close()
	
	#data
	for n in range(0,2*bins_nb):
		results = str(bins_labels[n])
		for part_index in range(0,len(particles_def["labels"])):
			part = particles_def["labels"][part_index]
			if particles_def_pres[part]:
				results += "	" + "{:.6e}".format(density_particles_pc[part][n])
			else:
				results += "	0"
		output_xvg.write(results + "\n")	
	output_xvg.close()
				
	return
def density_write_charges():											#DONE
	
	tmp_file = 'membrane_prop_TM_density_charges'
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/charge/' + str(tmp_file) + '.xvg'
	filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/charge/' + str(tmp_file) + '.txt'
	output_txt = open(filename_txt, 'w')
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_txt.write("@[charge density profile - written by membrane_prop v" + str(version_nb) + "]\n")
	output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_thickness_species.xvg.\n")
	output_xvg.write("# [charge density profile - written by membrane_prop v" + str(version_nb) + "]\n")
	output_xvg.write("# volume properties:\n")
	output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
	output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
	output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Charge density profile along z\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"charge density (e.Angstrom-3)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length " + str(len(charges_groups.keys()) + 1) + "\n")
	charges_groups_keys = charges_groups.keys()
	for q_index in range(0,len(charges_groups_keys)):
		charge_group = charges_groups_keys[q_index]
		output_xvg.write("@ s" + str(q_index) + " legend \"" + str(charge_group) + "\"\n")
		output_txt.write(str(tmp_file) + "," + str(q_index) + "," + str(charge_group) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(charges_colours[charge_group])) + "\n")
	output_xvg.write("@ s" + str(len(charges_groups.keys())) + " legend \"total\"\n")
	output_txt.write(str(tmp_file) + "," + str(len(charges_groups.keys())+1) + ",total," + mcolors.rgb2hex(mcolorconv.to_rgb(charges_colours["total"])) + "\n")
	output_txt.close()
	
	#data
	for n in range(0,2*bins_nb):
		results = str(bins_labels[n])
		for q_index in range(0,len(charges_groups_keys)):
			charge_g = charges_groups_keys[q_index]
			if charges_groups_pres[charge_g]:
				results += "	" + "{:.6e}".format(density_charges[charge_g][n])
			else:
				results += "	0"
		results += "	" + "{:.6e}".format(density_charges["total"][n])
		output_xvg.write(results + "\n")	
	output_xvg.close()
	
	return

def density_graph_particles():											#DONE
			
	#filenames
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/density/membrane_prop_TM_density_particles.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/density/membrane_prop_TM_density_particles.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Relative particles frequency along z")

	#plot data
	ax = fig.add_subplot(111)
	for part in particles_def["labels"]:
		if particles_def_pres[part]:
			plt.plot(np.arange(-args.slices_range,args.slices_range, args.slices_thick), density_particles_pc[part], color = particles_def["colour"][part], label = str(part))
	plt.vlines(z_upper, 0, max_density_particles_pc, linestyles = 'dashed')
	plt.vlines(z_lower, 0, max_density_particles_pc, linestyles = 'dashed')
	plt.vlines(0, 0, max_density_particles_pc, linestyles = 'dashdot')
	fontP.set_size("small")
	ax.legend(prop=fontP)
	plt.xlabel('z distance to bilayer center [$\AA$]')
	plt.ylabel('relative particles frequency [$\AA^{-3}$]')
	
	#save figure
	ax.set_xlim(-args.slices_range, args.slices_range)
	ax.set_ylim(0, max_density_particles_pc)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	
	return
def density_graph_charges():											#DONE
			
	#filenames
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/charge/membrane_prop_TM_density_charges.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/charge/membrane_prop_TM_density_charges.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Charge density profile in TM clusters")

	#plot data
	ax = fig.add_subplot(111)
	for charge_g in charges_groups.keys() + ["total"]:
		if charge_g == "total":
			plt.plot(np.arange(-args.slices_range,args.slices_range, args.slices_thick), density_charges[charge_g], color = charges_colours[charge_g], label = str(charge_g), linewidth = 2)
		elif charges_groups_pres[charge_g]:
			plt.plot(np.arange(-args.slices_range,args.slices_range, args.slices_thick), density_charges[charge_g], color = charges_colours[charge_g], label = str(charge_g))
	plt.vlines(z_upper, min_density_charges, max_density_charges, linestyles = 'dashed')
	plt.vlines(z_lower, min_density_charges, max_density_charges, linestyles = 'dashed')
	plt.vlines(0, min_density_charges, max_density_charges, linestyles = 'dashdot')
	plt.hlines(0,-args.slices_range, args.slices_range)
	fontP.set_size("small")
	ax.legend(prop=fontP)
	plt.xlabel('z distance to bilayer center [$\AA$]')
	plt.ylabel('average charge density [$e.\AA^{-3}$]')
	
	#save figure
	ax.set_xlim(-args.slices_range, args.slices_range)
	ax.set_ylim(min_density_charges, max_density_charges)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return

def angle_write_extent():

	#filename
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/angle/membrane_prop_angle_extent.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average local normal angle vs z axis - written by membrane_prop v" + str(version_nb) + "]\n")
	output_xvg.write("# nb of frames processed = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Evolution of average angle between local normal and z axis\"\n")
	output_xvg.write("@ xaxis label \"time (ns)\"\n")
	output_xvg.write("@ yaxis label \"average angle vs z axis (degrees)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 2\n")
	output_xvg.write("@ s0 legend \"angle (avg)\"\n")
	output_xvg.write("@ s1 legend \"angle (std)\"\n")
	
	#data
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index]) + "	" + "{:.6e}".format(theta_deg_avg[f_index]) + "	" + "{:.6e}".format(theta_deg_std[f_index])
		output_xvg.write(results + "\n")	
	output_xvg.close()

	return
def angle_graph_extent():
	
	#filenames
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/angle/membrane_prop_angle_extent.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/angle/membrane_prop_angle_extent.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of average angle between local normal and z axis")

	#plot data
	ax = fig.add_subplot(111)
	plt.plot(frames_time, theta_deg_avg, color = 'k', label = "avg", linewidth = 2)
	plt.fill_between(frames_time, theta_deg_avg - theta_deg_std, theta_deg_avg + theta_deg_std, color = '#A4A4A4', edgecolor = '#A4A4A4', linewidth = 0, alpha = 0.2)	
	#plt.hlines(0, min(frames_time), max(frames_time))
	fontP.set_size("small")
	ax.legend(prop=fontP)
	plt.xlabel('time (ns)')
	plt.ylabel('average angle vs z axis (degrees)')
	
	#save figure
	ax.set_xlim(0, max(frames_time))
	ax.set_ylim(0, 90)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	return

def angle_write_derivative():
	#filename
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/angle/membrane_prop_angle_derivatives.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [evolution of normal vector derivatives - written by membrane_prop v" + str(version_nb) + "]\n")
	output_xvg.write("# nb of frames processed = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Evolution of normal vector derivatives\"\n")
	output_xvg.write("@ xaxis label \"time (ns)\"\n")
	output_xvg.write("@ yaxis label \"d(component)/d(axis)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 3\n")
	output_xvg.write("@ s0 legend \"dnx/dx\"\n")
	output_xvg.write("@ s0 legend \"dny/dy\"\n")
	output_xvg.write("@ s0 legend \"dnz/dz\"\n")

	#data
	for f_index in range(0,nb_frames_to_process):
		results = str(frames_time[f_index]) + "	" + "{:.6e}".format(norm_dnxdx_avg[f_index]) + "	" + "{:.6e}".format(norm_dnydy_avg[f_index]) + "	" + "{:.6e}".format(norm_dnzdz_avg[f_index])
		output_xvg.write(results + "\n")	
	output_xvg.close()
	return
def angle_graph_derivative():
	
	#filenames
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/angle/membrane_prop_angle_derivatives.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/angle/membrane_prop_angle_derivatives.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Evolution of normal vector derivatives")

	#plot data: dnx/dx
	#-----------------
	ax1 = fig.add_subplot(311)
	plt.plot(frames_time, norm_dnxdx_avg, color = 'k', linewidth = 2)
	#plt.fill_between(frames_time, norm_dnxdx_avg - norm_dnxdx_std, norm_dnxdx_avg + norm_dnxdx_std, color = '#A4A4A4', edgecolor = '#A4A4A4', linewidth = 0, alpha = 0.2)
	plt.hlines(0, min(frames_time), max(frames_time), linestyle = "dashed")
	fontP.set_size("small")
	plt.xlabel('time (ns)')
	plt.ylabel('d(nx)/dx')	
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax1.xaxis.set_ticks_position('bottom')
	ax1.yaxis.set_ticks_position('left')
	ax1.set_xlim(0, max(frames_time))
	if args.nxmax != -1:
		ax1.set_ylim(top = args.nxmax)
	if args.nxmin != -1:
		ax1.set_ylim(bottom = args.nxmin)
	if args.nbx != -1:
		ax1.xaxis.set_major_locator(MaxNLocator(nbins=args.nbx))
	if args.nby != -1:
		ax1.yaxis.set_major_locator(MaxNLocator(nbins=args.nby))
	ax1.xaxis.labelpad = 20
	ax1.yaxis.labelpad = 20
	plt.setp(ax1.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax1.yaxis.get_majorticklabels(), fontsize = "small")
	
	#plot data: dny/dy
	#-----------------
	ax2 = fig.add_subplot(312)
	plt.plot(frames_time, norm_dnydy_avg, color = 'k', linewidth = 2)
	#plt.fill_between(frames_time, norm_dnydy_avg - norm_dnydy_std, norm_dnydy_avg + norm_dnydy_std, color = '#A4A4A4', edgecolor = '#A4A4A4', linewidth = 0, alpha = 0.2)
	plt.hlines(0, min(frames_time), max(frames_time), linestyle = "dashed")
	fontP.set_size("small")
	plt.xlabel('time (ns)')
	plt.ylabel('d(ny)/dy')	
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax2.xaxis.set_ticks_position('bottom')
	ax2.yaxis.set_ticks_position('left')
	ax2.set_xlim(0, max(frames_time))
	if args.nymax != -1:
		ax2.set_ylim(top = args.nymax)
	if args.nymin != -1:
		ax2.set_ylim(bottom = args.nymin)
	if args.nbx != -1:
		ax2.xaxis.set_major_locator(MaxNLocator(nbins=args.nbx))
	if args.nby != -1:
		ax2.yaxis.set_major_locator(MaxNLocator(nbins=args.nby))
	ax2.xaxis.labelpad = 20
	ax2.yaxis.labelpad = 20
	plt.setp(ax2.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax2.yaxis.get_majorticklabels(), fontsize = "small")
	
	#plot data: dnz/dz
	#-----------------
	ax3 = fig.add_subplot(313)
	plt.plot(frames_time, norm_dnzdz_avg, color = 'k', linewidth = 2)
	#plt.fill_between(frames_time, norm_dnzdz_avg - norm_dnzdz_avg, norm_dnzdz_avg + norm_dnzdz_std, color = '#A4A4A4', edgecolor = '#A4A4A4', linewidth = 0, alpha = 0.2)
	plt.hlines(0, min(frames_time), max(frames_time), linestyle = "dashed")
	fontP.set_size("small")
	plt.xlabel('time (ns)')
	plt.ylabel('d(nz)/dz')	
	ax3.spines['top'].set_visible(False)
	ax3.spines['right'].set_visible(False)
	ax3.xaxis.set_ticks_position('bottom')
	ax3.yaxis.set_ticks_position('left')
	ax3.set_xlim(0, max(frames_time))
	if args.nzmax != -1:
		ax3.set_ylim(top = args.nzmax)
	if args.nzmin != -1:
		ax3.set_ylim(bottom = args.nzmin)
	if args.nbx != -1:
		ax3.xaxis.set_major_locator(MaxNLocator(nbins=args.nbx))
	if args.nby != -1:
		ax3.yaxis.set_major_locator(MaxNLocator(nbins=args.nby))
	ax3.xaxis.labelpad = 20
	ax3.yaxis.labelpad = 20
	plt.setp(ax3.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax3.yaxis.get_majorticklabels(), fontsize = "small")

	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
#process inputs
#=========================================================================================

#data loading
set_lipids_beads()
set_particles()
if args.chargesfilename != "no":
	set_charges()
load_MDA_universe()
if args.selection_file_ff != "no":
	identify_ff()
identify_leaflets()

#create data structures
struct_time()
struct_data()

#=========================================================================================
# generate data
#=========================================================================================
print "\nCalculating density profiles..."

#case: structure only
if args.xtcfilename=="no":
	calculate_properties(U.trajectory.ts.dimensions, 0)
	print ""
	
#case: browse xtc frames
else:
	for f_index in range(0,nb_frames_to_process):
		#frame properties
		ts = U.trajectory[frames_to_process[f_index]]
		f_time = ts.time/float(1000)
		f_nb = ts.frame
		frames_nb[f_index] = f_nb
		frames_time[f_index] = f_time
		
		#calculate densities
		calculate_properties(U.trajectory.ts.dimensions, f_index)
	print ""

#=========================================================================================
# process data
#=========================================================================================
calculate_stats()

#=========================================================================================
# produce outputs
#=========================================================================================

print "\nWriting outputs..."
if args.normal != 'z':
	write_grid_statistics()
density_write_particles()
density_graph_particles()
angle_write_extent()
angle_graph_extent()
angle_write_derivative()
angle_graph_derivative()
if args.chargesfilename != "no":
	density_write_charges()
	density_graph_charges()
	
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
