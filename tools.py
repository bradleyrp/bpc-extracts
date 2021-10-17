#!/usr/bin/env python

import os
import sys
import json
import numpy as np
import h5py

def load(name,cwd=None,verbose=False,exclude_slice_source=False,filename=False):
	"""
	Get binary data from a computation.
	"""
	if not cwd: cwd,name = os.path.dirname(name),os.path.basename(name)
	cwd = os.path.abspath(os.path.expanduser(cwd))
	fn = os.path.join(cwd,name)
	if not os.path.isfile(fn): raise Exception('[ERROR] failed to load %s'%fn)
	data = {}
	rawdat = h5py.File(fn,'r')
	for key in [i for i in rawdat if i!='meta']: 
		if verbose:
			print('[READ] '+key)
			print('[READ] object = '+str(rawdat[key]))
		data[key] = np.array(rawdat[key])
	if 'meta' in rawdat: 
		if sys.version_info<(3,0): out_text = rawdat['meta'].value
		else: out_text = rawdat['meta'][()].decode()
		attrs = json.loads(out_text)
	else: 
		print('[WARNING] no meta in this pickle')
		attrs = {}
	if exclude_slice_source:
		for key in ['grofile','trajfile']:
			if key in attrs: del attrs[key]
	for key in attrs: data[key] = attrs[key]
	if filename: data['filename'] = fn
	rawdat.close()
	return data
