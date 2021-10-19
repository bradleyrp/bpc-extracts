#!/usr/bin/env python

import time
import argparse
import os
import numpy as np
from joblib import Parallel,delayed

from tools import load,store
from mesh import makemesh_regular
from looptools import status,framelooper
from tools import datmerge

def undulations(data,grid_spacing=0.5,nprocs=4):

	"""
	Compute bilayer midplane structures for studying undulations.
	"""
	vecs = data['vecs']
	nframes = data['nframes']
	trajectory = data['points']
	attrs,result = {},{}
	monolayer_indices = data['monolayer_indices']
	# choose grid dimensions
	grid = np.array([round(i) for i in np.mean(vecs,axis=0)/grid_spacing])[:2]
	# parallel
	mesh = [[],[]]
	for mn in range(2):
		start = time.time()
		mesh[mn] = Parallel(n_jobs=nprocs,verbose=0,require='sharedmem')(
			delayed(makemesh_regular)(
				trajectory[fr][np.where(monolayer_indices==mn)],vecs[fr],grid)
			for fr in framelooper(nframes,start=start,text='monolayer %d, frame'%mn))

	# pack
	result['mesh'] = np.array(mesh)
	result['grid'] = np.array(grid)
	result['nframes'] = np.array(nframes)
	result['vecs'] = vecs
	attrs['grid_spacing'] = grid_spacing
	return result,attrs	

# iterative reexection
import sys
__me__ = sys.argv[0]
assert __me__.endswith('.py')
def go(): exec(open(__me__).read(),globals())

if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		epilog='Convert lipid_abstractor data to undulations data.')
	parser.add_argument('-i',dest='input',required=True,
		help='Filename for the lipid abstractor (dat) file.')
	parser.add_argument('-o',dest='output',required=True,
		help='Filename for the undulations data output.')
	args = parser.parse_args()
	fn = args.input
	fn_out = args.output
	
	if os.path.isfile(fn_out):
		raise Exception('file exists: %s'%fn_out)
	if not os.path.isfile(fn):
		raise Exception('cannot find')
	if not os.path.isdir(os.path.dirname(fn_out)):
		raise Exception('directory for the output file does not exist')

	# load the lipid_abstractor data
	data = load(fn)
	# run the calculation
	result,attrs = undulations(data=data,grid_spacing=0.5,nprocs=4)
	# write the resulting data
	store(
		name=os.path.basename(fn_out),
		path=os.path.dirname(fn_out),
		obj=result,attrs=attrs)

