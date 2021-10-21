#!/usr/bin/env python

"""
Plot undulation spectra and height profiles.
"""

import os,sys,re
import argparse
import glob
import json
from looptools import basic_compute_loop,status
from undulate import calculate_undulations
from undulate_plot import undulation_panel,add_undulation_labels,add_axgrid,add_std_legend
from tools import load,sweeper,panelplot,picturesave,picturefind
from tools_plot import plothull
# height correlation requires scipy
import scipy
import scipy.spatial
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
str_types = [str]

def undulation_spectra(figname,style='all_on_one'):
	"""
	Manage different plot combinations.
	This plot covers the method used in Bradley-2016 and the same content in the thesis.
	Later, undulation spectra were generated for the dextran project at  
	plot-analyze_dextran_undulations.py. It would be useful to compare these to be
	sure the method is exactly the same.
	"""
	global sns
	# external settings
	# dev: define this globally: plotspecs = work.plots[plotname].get('specs',{})
	wavevector_limits = plotspecs.get('wavevector_limits',[1.0])
	fit_styles = [
		'band,perfect,simple',
		'band,perfect,basic',
		'band,perfect,curvefit',
		'band,blurry,curvefit',
		'band,perfect,curvefit-crossover',
		'band,perfect,fit',][:-2]
	# top level sweep over plot settings and styles
	# see below for an override that just runs the standard spectra
	# the following sweep is overwritten below and retained here only so you
	#   could look at the alternative fitting methods later. note that almost all
	#   of the fitting methods that do not cause an optimizer malfunction almost 
	#   almost always produce nearly identical kappa values. the ones that do not 
	#   (namely the band,perfect,curvefit-crossover) are probably not working because
	#   they try to incorporate surface tension or protrusion without enough
	#   constraints on the optimizer 
	sweep_specs = {
		'wavevector_limit':wavevector_limits,
		'style':['all_on_one'],
		# if you use average_normal you must run an additional calculation
		#   to compute the normal vectors at each grid point. this method
		#   is really only necessary if you are way outside of Monge gauge
		# the standard method is "flat" and the average method is not 
		#   strictly within the Monge gauge. the average method subtracts the
		#   average structure, and was used for comparison purposes only
		'midplane_method':['flat','average','average_normal'][:2],
		'fit_style':fit_styles,
		'residual_form':['log'],}
	# override the excessively large parameter sweep with the standard method
	sweep_specs = {
		'wavevector_limit':wavevector_limits,
		'style':['all_on_one'],
		'midplane_method':['flat'],
		'fit_style':['band,perfect,curvefit'],
		'residual_form':['log'],}
	# get colors first from art which always
	colors_reqs = ['binned','fitted','line']
	if colors!=None: 
		# must match the receiver in the plot function
		sweep_specs['color'] = [dict([(sn,dict([(key,colors[sn]) 
			for key in colors_reqs])) for sn in colors])]
	# we could also get the colors systematically from the metadata here?
	# fallback randomly generates random colors
	else: 
		def random_colormaker(x,n,name='jet'): 
			return mpl.cm.__dict__[name](np.linspace(0.,1.0,len(sns))[x]/n)
		colors_random = dict([(sn,dict([(key,random_colormaker(snum,len(sns)))
			for key in colors_reqs])) for snum,sn in enumerate(sns)])
		sweep_specs['color'] = [colors_random]
	plots_this = sweeper(**sweep_specs)
	# settings
	art = {'fs':{'legend':12}}
	# loop over all plot settings and styles
	for pnum,plotspec in enumerate(plots_this):
		# unpack settings
		style = plotspec['style']
		# router over plot layouts
		if style=='all_on_one':
			figsize = (5,5)
			layout = {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
			axes,fig = panelplot(layout,figsize=figsize)
			ax = axes[0]
			for snum,sn in enumerate(sns):
				# the crossover requires tension
				if plotspec['fit_style']=='band,perfect,curvefit-crossover':
					plotspec['fit_tension'] = True
				# log residuals with the more complicated fits causes weird results
				# note that these fits are not recommended so we raise
				if plotspec['fit_style'] in [
					'band,perfect,curvefit-crossover','band,perfect,fit']:
					plotspec['residual_form'] = 'linear'
					raise Exception('these fits are not tested')
				plot_undulation_spectrum(ax,sn,**plotspec)
			decorate_undulation_plot(ax=ax,art=art)
		else: raise Exception('invalid plot style %s'%style)
		# save the plot, with relevant keys
		meta_reqs = ['midplane_method','wavevector_limit','style','fit_style','residual_form']
		meta = dict([(key,plotspec[key]) for key in meta_reqs])
		# separate the filename and directory
		figname_abs = os.path.abspath(os.path.expanduser(figname))
		plotdir = os.path.dirname(figname)
		basename = os.path.basename(figname) 
		# we embed the metadata with the figure
		# if you are making multiple figures with the same name (and different metadata), you can 
		#   set version=True to avoid overwriting them with a version number on each file
		picturesave(basename,plotdir,backup=False,version=False,meta=meta)

def plot_undulation_spectrum(ax,sn,**kwargs):
	"""
	Plot a single undulation spectrum.
	"""
	global data
	mesh = data[sn]['data']['mesh']
	surf = mesh.mean(axis=0)
	surf = (surf - np.tile(surf.reshape(len(surf),-1).mean(axis=1),
		(surf.shape[1],surf.shape[2],1)).transpose(2,0,1))
	vecs = data[sn]['data']['vecs']
	# the standard fitting method is band,perfect,curvefit however
	#   the code has many different now-discarded options for performing 
	#   the fit. these methods were for reference only and the approved 
	#   method uses "perfect" binning which averages all q_x,q_y with the same
	#   magnitude onto a single point. the blurry version takes finite wavevector
	#   windows and produces similar results. the remaining fitting parameters
	#   refer only to the way the program performs the fit
	fit_style_default = 'band,perfect,curvefit'
	fit_style = kwargs.get('fit_style',fit_style_default)
	lims = kwargs.get('lims',[0.,kwargs.get('wavevector_limit',1.0)])
	colors = kwargs.get('color','k')
	midplane_method = kwargs.get('midplane_method','flat')
	colors_reqs = ['binned','fitted','line']
	if type(colors) in str_types: 
		colors = dict([(key,colors) for key in colors_reqs])
	elif not all([key in colors.get(sn,{}) for key in colors_reqs]): 
		raise Exception('need keys in colors %s'%colors_reqs)
	uspec = calculate_undulations(surf,vecs,fit_style=fit_style,lims=lims,
		midplane_method=midplane_method,fit_tension=kwargs.get('fit_tension',False))
	label = meta.get(sn,{}).get('label',sn)+'\n'+\
		r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'
	q_binned,energy_binned = uspec['q_binned'][1:],uspec['energy_binned'][1:]
	ax.plot(q_binned,energy_binned,'.',lw=0,markersize=10,markeredgewidth=0,
		label=None,alpha=0.2,color=colors[sn]['binned'])
	q_fit,energy_fit = np.transpose(uspec['points'])
	ax.plot(q_fit,energy_fit,'.',lw=0,markersize=4,markeredgewidth=0,
		label=label,alpha=1.,zorder=4,color=colors[sn]['fitted'])
	def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
		return 1.0/(area/1.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
	ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
		area=uspec['area']),lw=1,zorder=3,color=colors[sn]['line'])

def decorate_undulation_plot(ax,art):
	"""
	Uniform appearance for all undulation spectra.
	"""
	add_undulation_labels(ax,art=art)
	add_std_legend(ax,loc='upper right',art=art)
	add_axgrid(ax,art=art)

def compute_protein_proximity_height_correlation(fr,pbc=False):
	"""Calculation for basic_compute_loop to correlate membrane 
	height with proximity to proteins."""
	global surf,vecs,protein_pts
	# gather points
	surf = mesh[fr]
	ngrid = surf.shape
	if pbc: surf = np.tile(surf,(3,3))
	vec = vecs[fr]
	# unpack pts per the structure of points_all
	pts_this = np.concatenate(protein_pts[fr])[...,:2]
	# get XY points for the grid
	if not pbc:
		xypts = np.concatenate(np.transpose(np.meshgrid(*[np.linspace(0,v,n) 
			for v,n in zip(vec[:2],ngrid)])))
	else:
		xypts = np.concatenate(np.transpose(np.meshgrid(*[np.linspace(0,3*v,3*n) 
			for v,n in zip(vec[:2],ngrid)])))
		pts_this += vec[:2]
	# round because the cKDTree is finnicky
	pts_back = pts_this
	pts_fore = np.floor(xypts*10.**3)/10.**3
	# make the tree
	tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=vec[:2]*(3 if pbc else 1))
	close,nns = tree.query(pts_fore,k=1)
	# return the minimum distance to the protein and the bilayer height
	return np.array([close,surf.reshape(-1)]).T

def plot_height_proximity_correlation(**kwargs):
	"""
	Plot the instantaneous membrane height vs proximity to protein points.
	"""
	global sns
	import seaborn as sb
	# stash to globals to iterate the plot aesthetics
	if 'post' not in globals():
		global post,mesh,vecs,protein_pts
		post = {}
		sample_rate = 1
		for sn in sns:
			# points_all for the dimer simulations has dimensions frames, monomer, points, xyz
			protein_pts = data_prot[sn]['data']['points_all']
			vecs = data[sn]['data']['vecs']
			nframes = len(vecs)
			mesh = data[sn]['data']['mesh'].mean(axis=0)
			ngrid = mesh.shape[-2:]
			mesh -= np.tile(mesh.reshape(nframes,-1).mean(axis=1),
				(ngrid[0],ngrid[1],1)).transpose((2,0,1))
			incoming = basic_compute_loop(compute_protein_proximity_height_correlation,
				looper=[dict(fr=fr) for fr in range(0,nframes,sample_rate)])
			post[sn] = dict(sizes=[len(i) for i in incoming],incoming=np.concatenate(incoming))

	# regular plot
	binw = 1.0
	axes,fig = square_tiles(1,figsize=(8,8))
	ax = axes[0]
	pbc_spacing = min([min(data[sn]['data']['vecs'].mean(axis=0)[:2]) for sn in sns])
	colors = sb.color_palette("hls",len(sns))
	for snum,sn in enumerate(sns):
		rmax,zmax = [max([np.abs(v['incoming'][:,i]).max() for v in post.values()]) for i in range(2)]
		bins = np.arange(0,rmax+binw,binw)
		rate = 1
		sample = post[sn]['incoming'][::rate]
		binned = [sample[np.all((sample.T[0]>=bins[ii],
			sample.T[0]<=bins[ii+1]),axis=0)][:,1] for ii,i in enumerate(bins[:-1])]
		means = np.array([np.mean(i) for i in binned])
		stds = np.array([np.std(i) for i in binned])
		ax.plot(bins[:-1],means,label=sn,color=colors[snum])
		ax.fill_between(bins[:-1],means-stds,means+stds,alpha=0.1,color=colors[snum])
	# there is very little difference between doing the expensive PBC
	ax.set_xlim((0.,pbc_spacing/2.))
	ax.axhline(0,c='k',lw=1)
	plt.legend()
	plt.savefig(os.path.join(plotdir,'fig.height_proximity.png'))

def plot_height_profiles(figname_prefix):
	"""
	Plot the bilayer height profile with protein positions.
	"""
	global sns
	# one plot per simulation
	for sn in sns:
		mesh = data[sn]['data']['mesh']
		surf = mesh.mean(axis=0).mean(axis=0)
		surf -= surf.mean()
		hmax = np.abs(surf).max()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		vecs = data[sn]['data']['vecs'].mean(axis=0)
		im = ax.imshow(surf.T,origin='lower',
			interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
			extent=[0,vecs[0],0,vecs[1]],vmax=hmax,vmin=-1*hmax)
		if sn in data_prot:
			# try:
			if 1:
				points_all = data_prot[sn]['data']['points_all'] 
				points_all_mean_time = points_all.mean(axis=0)
				plothull(ax,points_all_mean_time[...,:2],griddims=surf.shape,vecs=vecs,c='k',lw=0)	
			#! except: status('failed to get protein points for %s'%sn,tag='warning')
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		axins = inset_axes(ax,width="5%",height="100%",loc=3,
			bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
		cbar = plt.colorbar(im,cax=axins,orientation="vertical")
		cbar.set_label(r"$\mathrm{\langle z \rangle (nm)}$",labelpad=-40,y=1.05,rotation=0)
		axins.tick_params(axis='y',which='both',left='off',right='off',labelright='on')
		ax.set_title('average bilayer height')
		ax.set_xlabel('x (nm)')
		ax.set_ylabel('y (nm)')
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		# separate the filename and directory
		figname_abs = os.path.abspath(os.path.expanduser(figname_prefix))
		plotdir = os.path.dirname(figname_abs)
		basename = os.path.basename(figname_abs) 
		# use version and meta to embed a dictionary of metadata and write versioned png files
		# by default we append the simulation code to the basename
		picturesave(basename+'-'+sn,plotdir,backup=False,version=False,meta={})

def load_data(calcnames,post_dn):
	"""
	Find and load data for upstream calculations (calcnames) from a directory (post_dn).
	Data must follow the BioPhysCode naming convention.
	"""
	# parse specs
	fns = [i for i in glob.glob(os.path.join(post_dn,'*.spec'))]
	toc = {}
	for fn in fns:
		with open(fn) as fp:
			spec = json.load(fp)
		sn = spec['meta']['sn']
		calcname = spec['calc']['name']
		if calcname not in calcnames: continue
		if (sn,calcname) not in toc:
			toc[(sn,calcname)] = []
		spec['fn'] = fn
		toc[(sn,calcname)].append(spec)

	def get_single_spec(toc,sn,calc_name):
		candidates = toc.get((sn,calc_name),[])
		if not candidates:
			raise Exception('cannot find sn=%s, calc_name=%s'%(sn,calc_name))
		elif len(candidates)>1:
			raise Exception('redundant calculations for sn=%s, calc_name=%s'%(sn,calc_name))
		else: return candidates[0]

	# load data
	data = {}
	for calcname in calcnames:
		for sn in sns:
			key = (sn,calcname)
			spec = get_single_spec(toc,sn,calcname)
			fn = spec['fn']
			fn_dat = re.sub('.spec','.dat',fn)
			# load the h5 data
			dat = load(fn_dat)
			# nest the data when we need multiple upstream calculations
			# note that the extra 'data' key is left in for backwards compatibility
			if len(calcnames)>1:
				if 'data' not in data[sn]: data[sn] = {'data':{}}
				data[sn]['data'][calcname] = dat
			# no nesting if omly one calculation is required
			else:
				data[sn] = {'data':dat}
			print('[STATUS] loaded %s'%os.path.basename(fn_dat))
	return data

# iterative reexection
import sys
__me__ = sys.argv[0]
assert __me__.endswith('.py')
def go(): exec(open(__me__).read(),globals())

if __name__=='__main__': 

	# SETTINGS		
	# target upstream calculation name
	calcnames = ['undulations']
	# specifications for the upstream data
	plotspecs = {'grid_spacing': 0.5}
	# define metadata including labels for the plots (this is an example)
	meta = {'v1021':{'label':'v1021'}}
	# optional color definitions
	colors = None
	
	# arguments
	modes = ['spectra','height-profiles']
	parser = argparse.ArgumentParser(
		epilog='Generate undulation spectra.')
	parser.add_argument('-s',dest='sources',required=True,
		help='List of undulations dat files.',nargs='+')
	parser.add_argument('-o',dest='output',required=True,
		help='Output plot name.')
	parser.add_argument('-m',dest='mode',required=True,
		help='Mode (one of: %s).'%(', '.join(modes)))
	args = parser.parse_args()
	if args.mode not in modes:
		raise Exception('mode must be one of: %s'%(', '.join(modes)))
	if args.mode=='spectra':
		if not args.output.endswith('.png'):
			raise Exception('the output must be a png')
		args.output = args.output.rstrip('.png')
	if args.mode=='height-profiles':
		if args.output.endswith('.png'):
			raise Exception('output for height-profiles mode '
				'should be a file prefix (we will add png)')
	for fn in args.sources:
		if not os.path.isfile(fn):
			raise Exception('missing file in the sources list: %s'%fn)

	def infer_name(fn):
		"""
		The name of the simulation precedes the first dot in the filename.
		This is a placeholder function to infer the simulation name from e.g. a dat file.
		"""
		return re.match('^(.*?)\.',os.path.basename(fn)).group(1)

	def infer_protein_file(fn):
		"""Infer the associated `protein_abstractor` dat file from the existing `undulations` dat file
		whenever we want the protein positions (i.e. plotting the height profiles)."""
		fn_alt = re.sub('.undulations.','.protein_abstractor.',fn)		
		if not os.path.isfile(fn_alt):
			raise Exception('failed to infer protein_abstractor for: %s'%fn)
		return fn_alt

	# load data once if interactive
	if 'data' not in globals():

		global data
		data = {}
		for fn in args.sources:
			sn = infer_name(fn)
			data[sn] = {'data':load(fn)}
		# define the simulation names from the existing dataset
		sns = sorted(data.keys())
		
		# if we are making height profiles, and the simulations have proteines,
		#   then we need protein data
		if args.mode=='height-profiles':
			data_prot = {}
			for fn in args.sources:
				fn_prot = infer_protein_file(fn)
				sn = infer_name(fn_prot)
				data_prot[sn] = {'data':load(fn_prot)}
		
	if args.mode=='spectra':
		undulation_spectra(figname=args.output)	
	if args.mode=='height-profiles':
		plot_height_profiles(figname_prefix=args.output)
