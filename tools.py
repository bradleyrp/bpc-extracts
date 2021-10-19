#!/usr/bin/env python

import os
import sys
import re
import json
import numpy as np
import h5py
import glob
import collections
from PIL import Image
from PIL import PngImagePlugin
from looptools import status
str_types = [str]

# extractted from omni/base/store.py

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

def store(obj,name,path,attrs=None,print_types=False,verbose=True):
	"""
	Use h5py to store a dictionary of data.
	"""
	import h5py
	if type(obj) != dict: raise Exception('except: only dictionaries can be stored')
	if os.path.isfile(path+'/'+name): raise Exception('except: file already exists: '+path+'/'+name)
	path = os.path.abspath(os.path.expanduser(path))
	if not os.path.isdir(path): os.mkdir(path)
	fobj = h5py.File(path+'/'+name,'w')
	for key in obj.keys(): 
		if print_types: 
			print('[WRITING] '+key+' type='+str(type(obj[key])))
			print('[WRITING] '+key+' dtype='+str(obj[key].dtype))
		# python3 cannot do unicode so we double check the type
		if (type(obj[key])==np.ndarray and re.match('^str|^unicode',obj[key].dtype.name) 
			and 'U' in obj[key].dtype.str):
			obj[key] = obj[key].astype('S')
		try: dset = fobj.create_dataset(key,data=obj[key])
		except: 
			# multidimensional scipy ndarray must be promoted to a proper numpy list
			try: dset = fobj.create_dataset(key,data=obj[key].tolist())
			except: raise Exception("failed to write this object so it's probably not numpy"+
				"\n"+key+' type='+str(type(obj[key]))+' dtype='+str(obj[key].dtype))
	if attrs != None: 
		try: fobj.create_dataset('meta',data=np.string_(json.dumps(attrs)))
		except Exception as e: raise Exception('failed to serialize attributes: %s'%e)
	if verbose: status('[WRITING] '+path+'/'+name)
	fobj.close()

# extracted from omni/datapack.py

def delve(o,*k): 
	"""
	Return items from a nested dict.
	"""
	return delve(o[k[0]],*k[1:]) if len(k)>1 else o[k[0]]

# extracted from omni/base/hypothesis.py

import copy

def hypothesis(sweep,default=None):
	"""
	Code for sweeping an arbitrarily deep dictionary over many dimensions in combinations.
	Adapted from hypothesize so that no default is necessary. Only changed the name slightly for readability.
	"""
	# create the default hypothesis
	if default==None: default = {}
	for pathway in sweep:
		if pathway['route'][0] not in default: default[pathway['route'][0]] = {}
		for i in range(1,len(pathway['route'])-1):
			level = delve(default,*pathway['route'][:i])
			if pathway['route'][i] not in level: level[pathway['route'][i]] = {}
		if len(pathway['route'])>1: 
			delve(default,*pathway['route'][:-1])[pathway['route'][-1]] = list(pathway['values'])[0]
	for i in default: default[i] = None if default[i] == {} else default[i]

	# extract a list of lists of parameters to sweep over
	t = [i['values'] for i in sweep]
	# note that this non-numpythonic way of doing this has not been rigorously tested
	# note that the previous meshgrid method did not work on all types
	allcombos = list([[i] for i in t[0]])
	for s in t[1:]:
		for bi in range(len(allcombos)):
			b = allcombos.pop(0)
			for r in list(s): allcombos.append(b + [r])

	# assemble a list of hypotheses from all possible combinations of the sweep values
	# note that this code is general, and works for an arbitrarily deep dictionary
	hypotheses = []
	# for each combo generate a new hypothesis
	for combo in allcombos:
		# start with the default hypothesis
		newhypo = copy.deepcopy(default)
		# each combo has a value and a route which is a sequence of dictionary keys
		#   we loop over each route to set each final value for the sweep
		for routenum in range(len(sweep)):
			# to get to the deepest part of that route we use tmp as a pointer
			#   and iteratively traverse one level until the second to last level
			tmp = newhypo[sweep[routenum]['route'][0]]
			# the following checks if we are already at the end of the dictionary 
			if type(newhypo[sweep[routenum]['route'][0]]) != dict:
				newhypo[sweep[routenum]['route'][0]] = combo[routenum]
			else:
				for i in sweep[routenum]['route'][1:-1]: tmp = tmp[i]
				# at the final level, we now have a pointer to the lowest dictionary to set the value
				tmp[sweep[routenum]['route'][-1]] = combo[routenum]
		# once we set all the values, the hypothesis is ready
		hypotheses.append(newhypo)	
	return hypotheses

def sweeper(**kwargs):
	"""
	A simple method for sweeping many parameters. Provides a simple interface to the hypothesis function
	but does not allow arbitrary nested dictionaries. Instead of supplying explicit routes and values, the
	keys in kwargs are the routes (hence only one is possible).
	"""
	return hypothesis([dict(route=[key],values=val) for key,val in kwargs.items()])

# extracted from omni/plotter/panels.py

import matplotlib as mpl 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import itertools
import numpy as np

"""
EXAMPLE LAYOUT 

lay = {
	'out':{'grid':[2,1],'hratios':[1,1.2]},
	'ins':[
		{'grid':[1,3],'wspace':0.5,'hspace':0.2},
		{'grid':[3,4],'wspace':0.75,'hspace':0.75},
		],}

Layout must be a dictionary with 'out' and 'ins' keys.
The 'out' object should be a dictionary with a 'grid' item which is a tuple specifying columns and rows.
The 'ins' object should be a list of dictionaries that each have their own 'grid' item (rows by columns).
Use wspace and hspace to set the spacing between axes within a subplot.
Use hratios and wratios to specify the ratio of the outside object.
All grids of axes are placed in order by column, then row, so you don't need to index the exact position.
The returned axes object only requires two indices: the outside index and the inside one.
If there is only one outer object (e.g. 'out':{'grid':[2,1]}) then you do not need the outside index number.
"""

def panelplot(layout=None,figsize=(8,8),explicit_indexing=False):
	"""
	Nested panel plots.
	"""
	# default is a single plot
	lay = layout if layout != None else {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
	axes,axpos = [],[]
	fig = plt.figure(figsize=figsize)
	onrows,oncols = lay['out']['grid']
	outer_hspace = 0.45 if 'hspace' not in lay['out'] else lay['out']['hspace']
	outer_wspace = 0.45 if 'wspace' not in lay['out'] else lay['out']['wspace']
	outer_grid = gridspec.GridSpec(onrows,oncols,wspace=outer_wspace,hspace=outer_hspace)
	if 'hratios' in lay['out']: outer_grid.set_height_ratios(lay['out']['hratios'])
	if 'wratios' in lay['out']: outer_grid.set_width_ratios(lay['out']['wratios'])
	if type(lay['ins'])==dict: lay_ins = [lay['ins'] for i in range(np.product(lay['out']['grid']))]
	else: lay_ins = lay['ins']
	for ii,spot in enumerate(list(itertools.product(*[np.arange(i) for i in lay['out']['grid']]))):
		if ii>len(lay_ins)-1: raise Exception('looks like you have too few ins')
		hspace = lay_ins[ii]['hspace'] if 'hspace' in lay_ins[ii] else None
		wspace = lay_ins[ii]['wspace'] if 'wspace' in lay_ins[ii] else None
		inner_grid = gridspec.GridSpecFromSubplotSpec(*lay_ins[ii]['grid'],
			wspace=wspace,hspace=hspace,subplot_spec=outer_grid[ii])
		if 'hratios' in lay_ins[ii]: inner_grid.set_height_ratios(lay_ins[ii]['hratios'])
		if 'wratios' in lay_ins[ii]: inner_grid.set_width_ratios(lay_ins[ii]['wratios'])
		inaxs = [fig.add_subplot(j) for j in inner_grid]
		axpos.append(list(itertools.product(*[np.arange(i) for i in lay_ins[ii]['grid']])))
		axes.append(inaxs)
	return (axes[0] if len(axes)==1 and not explicit_indexing else axes,fig)

def square_tiles(ntiles,figsize,favor_rows=False,wspace=None,hspace=None):
	"""
	Create a grid of tiles with sequential order.
	"""
	nrows = ncols = int(np.ceil(np.sqrt(ntiles)))
	nrows -= int(1*(ntiles<=(nrows-1)*ncols))
	if not favor_rows: nrows,ncols = ncols,nrows
	layout = {'out':{'grid':[1,1]},'ins':{'grid':[nrows,ncols]}}
	if wspace: layout['ins']['wspace'] = wspace
	if hspace: layout['ins']['hspace'] = hspace
	# send a single number and we will make the plot proportional (however this might 
	#   not be perfect if tiles are not square)
	if type(figsize) not in [tuple,list]:
		figsize = tuple([figsize*ncols/max([nrows,ncols]),figsize*nrows/max([nrows,ncols])])
	axes,fig = 	panelplot(figsize=figsize,layout=layout)
	for i in range(nrows*ncols-ntiles): fig.delaxes(axes[-1*(i+1)])
	return axes[:ntiles],fig

# extracted from omni/base/store.py

def picturesave(savename,directory='./',meta=None,extras=[],backup=False,
	dpi=300,form='png',version=False,pdf=False,tight=True,pad_inches=0,figure_held=None,loud=True,
	redacted=False):
	"""
	Function which saves the global matplotlib figure without overwriting.
	!Note that saving tuples get converted to lists in the metadata so if you notice that your plotter is not 
	overwriting then this is probably why.
	"""
	#! amazing bug: if you keep a comma after meta it makes it a tuple and then there must be a 
	#!   one-way conversion to dict when it is written to the metadata of the image and this causes
	#!   the figure counts to keep increasing no matter what. a very subtle error! corrected below
	if type(meta)==tuple:
		if len(meta)!=1 or type(meta[0])!=dict: raise Exception('meta must be a dict')
		else: meta = meta[0]
	# automatically share images with group members (note that you could move this to config)
	os.umask(0o002)
	# earlier import allows users to set Agg so we import here, later
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	# intervene here to check the wordspace for picture-saving "hooks" that apply to all new pictures
	#! is it necessary to pass the workspace here?
	if 'work' in globals() and 'picture_hooks' in work.metadata.variables:
		extra_meta = work.metadata.variables['picture_hooks']
		# redundant keys are not allowed: either they are in picture_hooks or passed to picturesave
		redundant_extras = [i for i in extra_meta if i in meta]
		if any(redundant_extras):
			raise Exception(
				'keys "%r" are incoming via meta but are already part of picture_hooks'
				%redundant_extras)
	# redacted figures have blurred labels
	if redacted:
		directory_redacted = os.path.join(directory,'REDACTED')
		if not os.path.isdir(directory_redacted): os.mkdir(directory_redacted)
		directory = directory_redacted
		status('you have requested redacted figures, so they are saved to %s'%directory,tag='warning')
		import random
		color_back = work.metadata.director.get('redacted_background_color','')
		color_fore = work.metadata.director.get('redacted_foreground_color','k')
		if 'redacted_scrambler' in work.metadata.director:
			scrambler_code = work.metadata.director['redacted_scrambler']
			try: 
				scrambler = eval(scrambler_code)
				scrambler('test text')
			except: raise Exception(
				'failed to evaluate your `redacted_scrambler` from the director: `%s`'%scrambler_code)
		else: 
			#! method below is deprecated because it looks silly. best to use hashes
			if False: scrambler = lambda x,max_len=12:''.join([
				chr(ord('a')+random.randint(0,25)) for i in x][:max_len])
			scrambler = lambda x,max_len=10:('#'*len(x))[:max_len]
		num_format = re.compile("^[\-]?[1-9][0-9]*\.?[0-9]+$")
		isnumber = lambda x:re.match(num_format,x)
		for obj in [i for i in plt.findobj() if type(i)==mpl.text.Text]:
		    text_this = obj.get_text()
		    if text_this!='' and not isnumber(text_this):
		        obj.set_text(scrambler(text_this))
		        if color_back: obj.set_backgroundcolor(color_back)
		        obj.set_color(color_fore)
	# if version then we choose savename based on the next available index
	if version:
		# check for this meta
		search = picturefind(savename,directory=directory,meta=meta,loud=loud)
		if not search:
			if meta == None: raise Exception('[ERROR] versioned image saving requires meta')
			fns = glob.glob(os.path.join(directory,savename+'.v*'))
			nums = [int(re.findall('^.+\.v([0-9]+)\.png',fn)[0]) for fn in fns 
				if re.match('^.+\.v[0-9]+\.png',fn)]
			ind = max(nums)+1 if nums!=[] else 1
			savename += '.v%d'%ind
		else: savename = re.findall('(.+)\.[a-z]+',os.path.basename(search))[0]
	# backup if necessary
	savename += '.'+form
	base_fn = os.path.join(directory,savename)
	if loud: status('saving picture to %s'%savename,tag='store')
	if os.path.isfile(base_fn) and backup:
		for i in range(1,100):
			latestfile = '.'.join(base_fn.split('.')[:-1])+'.bak'+('%02d'%i)+'.'+base_fn.split('.')[-1]
			if not os.path.isfile(latestfile): break
		if i == 99 and os.path.isfile(latestfile):
			raise Exception('except: too many copies')
		else: 
			if loud: status('backing up '+base_fn+' to '+latestfile,tag='store')
			os.rename(base_fn,latestfile)
	# intervene to use the PDF backend if desired
	#   this is particularly useful for the hatch-width hack 
	#   (search self.output(0.1, Op.setlinewidth) in 
	#   python2.7/site-packages/matplotlib/backends/backend_pdf.py and raise it to e.g. 3.0)
	if pdf and form!='png': raise Exception('can only use PDF conversion when writing png')
	elif pdf:
		alt_name = re.sub('.png$','.pdf',savename)
		# holding the figure allows other programs e.g. ipython notebooks to show and save the figure
		(figure_held if figure_held else plt).savefig(alt_name,dpi=dpi,bbox_extra_artists=extras,
			bbox_inches='tight' if tight else None,pad_inches=pad_inches if pad_inches else None,format=form)
		# convert pdf to png
		os.system('convert -density %d %s %s'%(dpi,alt_name,base_fn))
		os.remove(alt_name)
	else: 
		(figure_held if figure_held else plt).savefig(base_fn,dpi=dpi,bbox_extra_artists=extras,
			bbox_inches='tight' if tight else None,pad_inches=pad_inches if pad_inches else None,format=form)
	plt.close()
	# add metadata to png
	if form=='png' and meta!=None:
		im = Image.open(base_fn)
		imgmeta = PngImagePlugin.PngInfo()
		imgmeta.add_text('meta',json.dumps(meta))
		im.save(base_fn,form,pnginfo=imgmeta)
	else: print('[WARNING] you are saving as %s and only png allows metadata-versioned pictures'%form)
	return base_fn

def picturesave_redacted(*args,**kwargs):
	"""Wrap picturesave with redacted plots."""
	return picturesave(*args,redacted=True,**kwargs)

def picturedat(savename,directory='./',bank=False):
	"""
	Read metadata from figures with identical names.
	"""
	directory = os.path.join(directory,'')
	if not bank: 
		if os.path.isfile(directory+savename): 
			try: return json.loads(Image.open(directory+savename).info['meta'])
			except: raise Exception('failed to load metadata from (possibly corrupted) image %s'%(
				directory+savename))
		else: return
	else:
		dicts = {}
		if os.path.isfile(directory+savename):
			dicts[directory+savename] = Image.open(directory+savename).info
		for i in range(1,100):
			base = directory+savename
			latestfile = '.'.join(base.split('.')[:-1])+'.bak'+('%02d'%i)+'.'+base.split('.')[-1]
			if os.path.isfile(latestfile): dicts[latestfile] = json.loads(Image.open(latestfile).info)
		return dicts

def lowest_common_dict_denominator(data):
	"""..."""
	if sys.version_info>=(3,0) and isinstance(data,str): return str(data)
	elif sys.version_info<(3,0) and isinstance(data,basestring): return str(data)
	elif isinstance(data,collections.Mapping): 
		if sys.version_info<(3,0): return dict(map(lowest_common_dict_denominator,data.iteritems()))
		else: return dict(map(lowest_common_dict_denominator,data.items()))
	elif isinstance(data,collections.Iterable): 
		return type(data)(map(lowest_common_dict_denominator,data))
	else: return data

def compare_dicts(a,b):
	"""Compare dictionaries with unicode strings."""
	return lowest_common_dict_denominator(a)==lowest_common_dict_denominator(b)

def picturefind(savename,directory='./',meta=None,loud=True):
	"""
	Find a picture in the plot repository.
	"""
	if loud: status('searching pictures',tag='store')
	regex = '^.+\.v([0-9]+)\.png'
	fns = glob.glob(directory+'/'+savename+'.v*')
	nums = map(lambda y:(y,int(re.findall(regex,y)[0])),filter(lambda x:re.match(regex,x),fns))
	matches = [fn for fn,num in nums if 
		compare_dicts(meta,picturedat(os.path.basename(fn),directory=directory))]
	if len(matches)>1 and meta!=None: 
		print('[ERROR] multiple matches found for %s'%savename)
		raise Exception('???')
	if matches==[] and meta==None:
		return dict([(os.path.basename(fn),
			picturedat(os.path.basename(fn),directory=directory)) for fn,num in nums]) 
	return matches if not matches else matches[0]

def datmerge(kwargs,name,key,same=False):
	"""
	Incoming upstream data are sometimes taken from multiple pickles.
	This function stitches together the key from many of these pickles.
	"""
	#---if there is only one upstream object with no index we simply lookup the key we want
	if name in kwargs.get('upstream',[]): return kwargs['upstream'][name][key]
	else:
		#---! this function seems to require upstream data so we except here
		if 'upstream' not in kwargs: raise Exception('datmerge needs upstream pointers')
		#---get indices for the upstream object added by computer
		inds = [int(re.findall(name+'(\d+)',i)[0]) for i in kwargs['upstream'] if re.match(name,i)]
		collected = [kwargs['upstream']['%s%d'%(name,i)][key] for i in sorted(inds)]
		if not same:
			if collected==[]: raise Exception('collected is empty, argument to datmerge is wrong')
			if type(collected[0])==list: return [i for j in collected for i in j]
			elif type(collected[0])==numpy.ndarray: 
				#---sometimes single numbers are saved as 0-dimensional arrays
				if numpy.shape(collected[0])==(): return numpy.array([float(i) for i in collected])
				else: return numpy.concatenate(collected)
			else: raise Exception('\n[ERROR] not sure how to concatenate')
		else:
			#---the same flag forces a check that the key is the same in each item of the collected data
			if any([any(collected[0]!=c) for c in collected[1:]]): 
				raise Exception('\n[ERROR] objects not same')
			else: return collected[0]

