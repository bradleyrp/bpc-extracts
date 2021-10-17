#!/usr/bin/env python

"""
Looping tools.
"""

import sys,time
import joblib
# some functions extracted manually below: from base.tools import framelooper

def status(string,i=0,looplen=None,bar_character=None,width=None,spacer='.',
	bar_width=25,tag='status',start=None,pad=None,refresh=True):
	"""
	Show a status bar and counter for a fixed-length operation.
	Taken from AUTOMACS to work in python 2 and 3.
	!NOTE need to fix the thing where lines get shorter and garbage is left behind...
	"""
	#---! it would be useful to receive a signal here to suppress the status bar from 
	#---! ...printing to the log file on backrun.
	#---use unicode if not piping to a log file
	logfile = (not hasattr(sys.stdout,'isatty')) or sys.stdout.isatty()==False
	#---use of equals sign below is deprecated when we suppress status bars in the log file below
	if not logfile: 
		left,right,bb = u'\u2590',u'\u258C',(u'\u2592' if bar_character==None else bar_character)
	else: left,right,bb = '|','|','='
	string = '[%s] '%tag.upper()+string if tag != '' else string
	if width: string = string.ljust(width,spacer)[:width]
	if pad: string = ('%-'+str(int(pad))+'s')%string
	if not looplen:
		if not logfile: sys.stdout.write(string+'\n')
		else: sys.stdout.write(string+'\n')
	elif looplen and logfile and i==0: sys.stdout.write('[STATUS] running a loop ')
	#---suppress progress bar in the log file except on the last item
	elif looplen and logfile and i>0 and i<looplen-1: sys.stdout.write('.')
	else:
		if start != None:
			esttime = (time.time()-start)/(float(i+1)/looplen)
			timestring = ' %s minutes'%str(abs(round((esttime-(time.time()-start))/60.,1)))
			bar_width = 15
		else: timestring = ''
		countstring = str(i+1)+'/'+str(looplen)
		bar = ' %s%s%s '%(left,int(bar_width*(i+1)/looplen)*bb+' '*\
			(bar_width-int(bar_width*(i+1)/looplen)),right)
		if not logfile: 
			output = (u'\r' if refresh else '')+string+bar+countstring+timestring+' '
			if sys.version_info<(3,0): output = output.encode('utf-8')
			if refresh:
				sys.stdout.flush()
				sys.stdout.write(output)
			else: print(output)
		else: 
			#---suppressed progress bar in the logfile avoids using carriage return
			sys.stdout.write('[STATUSBAR] '+string+bar+countstring+timestring+' ')
		if i+1<looplen: sys.stdout.flush()
		else: sys.stdout.write('\n')

def framelooper(total,start=None,text='frame'):
	"""
	When performing parallel calculations with joblib we pass a generator to count the number of 
	tasks and report the time.
	"""
	for fr in range(total):
		status(text,i=fr,looplen=total,tag='parallel',start=start)
		yield fr

def basic_compute_loop(compute_function,looper,run_parallel=True,debug=None):
	"""
	Canonical form of the basic compute loop.
	!!! remove this from contacts.py when it works
	"""
	#---send the frame as the debug argument
	if debug!=None and debug!=False:
		fr = debug
		incoming = compute_function(**looper[fr])
		import ipdb;ipdb.set_trace()
		sys.quit()
	start = time.time()
	if run_parallel:
		import joblib
		from joblib import Parallel,delayed
		#! your joblib is probably fine
		if 0 and joblib.__version__<0.12:
			from joblib.pool import has_shareable_memory
			incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
				delayed(compute_function,has_shareable_memory)(**looper[ll]) 
				for ll in framelooper(len(looper),start=start))
		else:
			incoming = Parallel(n_jobs=8,verbose=10 if debug else 0,require='sharedmem')(
				delayed(compute_function)(**looper[ll]) 
				for ll in framelooper(len(looper),start=start))
		incoming = []
		for ll in framelooper(len(looper)):
			incoming.append(compute_function(**looper[ll]))
	return incoming
