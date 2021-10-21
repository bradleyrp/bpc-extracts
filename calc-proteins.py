#!/usr/bin/python

import time
import os
from numpy import *
import argparse
import MDAnalysis
from joblib import Parallel,delayed
from looptools import status,framelooper
from timer import checktime
from mesh import *
from tools import store

def protein_abstractor(grofile,trajfile,**kwargs):

	"""
	PROTEIN ABSTRACTOR
	Compute the centroids of proteins in a simulation.
	"""

	#---unpack
	# sn = kwargs['sn']
	# work = kwargs['workspace']
	parallel = kwargs.get('parallel',False)
	#---MDAnalysis uses Angstroms not nm
	lenscale = 10.
	
	#---get protein coms here
	uni = MDAnalysis.Universe(grofile,trajfile)
	#---! cgmd removed here sel = uni.select_atoms(work.vars['selectors']['protein_selection'])
	sel = uni.select_atoms('protein')
	# nprots = work.meta.get(sn,{}).get('nprots',1)
	# ...!!! HARDCODED 
	nprots = 1
	beads_per_protein = len(sel.resids)/nprots
	nframes = len(uni.trajectory)
	inds = [arange(i*beads_per_protein,(i+1)*beads_per_protein) for i in range(nprots)]
	trajectory,trajectory_all,vecs = [],[],[]
	start = time.time()
	for fr in range(nframes):
		status('collecting protein centroids',i=fr,looplen=nframes,start=start,tag='compute')
		uni.trajectory[fr]
		#---center of geometry not centroid because masses are all 72 in martini
		pts = sel.positions[array(inds).astype(int)]/lenscale
		pts_mean = pts.mean(axis=0)
		trajectory.append(pts_mean)
		trajectory_all.append(pts)
		vecs.append(sel.dimensions[:3])

	#---pack
	attrs,result = {},{}
	result['resnames'] = array(sel.residues.resnames)
	result['names'] = array(sel.atoms.names)
	result['vecs'] = array(vecs)/lenscale
	result['nframes'] = array(nframes)
	result['points'] = array(trajectory)
	result['points_all'] = array(trajectory_all)
	return result,attrs	

# iterative reexection
import sys
__me__ = sys.argv[0]
assert __me__.endswith('.py')
def go(): exec(open(__me__).read(),globals())

if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		epilog='Extract protein positions from a GROMACS simulation.')
	parser.add_argument('-i',dest='input',required=True,
		help='Filename prefix for GROMACS trajectories (gro,xtc).')
	parser.add_argument('-o',dest='output',required=True,
		help='Filename for the protein_abstractor data.')
	args = parser.parse_args()
	fn_prefix = args.input
	fn_out = args.output
	
	if os.path.isfile(fn_out):
		raise Exception('file exists: %s'%fn_out)
	if not os.path.isdir(os.path.dirname(fn_out)):
		raise Exception('directory for the output file does not exist')

	gro_fn = fn_prefix+'.gro'
	xtc_fn = fn_prefix+'.xtc'
	if not os.path.isfile(gro_fn):
		raise Exception('cannot find %s'%gro_fn)
	if not os.path.isfile(xtc_fn):
		raise Exception('cannot find %s'%xtc_fn)

	result,attrs = protein_abstractor(grofile=gro_fn,trajfile=xtc_fn)
	# write the resulting data
	store(
		name=os.path.basename(fn_out),
		path=os.path.dirname(fn_out),
		obj=result,attrs=attrs)
