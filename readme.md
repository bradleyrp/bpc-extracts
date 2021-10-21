# Notes: how to compute undulations without BioPhysCode

These notes will explain how to extract the undulations code from the BioPhysCode workflow so you can use it as a standalone script.

## Status and Objective

I have reproduced the undulations code for the "banana" proteins (IRSp53 a.k.a. I-BAR and endophilin) on the local cluster named Rockfish. Our current objective is to strip the undulations code from the data structures and code in which it is embedded (which we previously called "BioPhysCode", see https://github.com/biophyscode). The undulations code can then be reused elsewhere.

## Existing workflow

The existing workflow uses code from two git repositories:

1. The "legacy factory" at `http://github.com/biophyscode/factory`
2. The calculation codes at `https://github.com/bradleyrp/omni-single`

To set up the existing workflow, you can follow a complete set of instructions [packaged with the calculations](https://github.com/bradleyrp/omni-single/blob/master/readme.md). These instructions cover all of the code and configuration. 

While the instructions include all of the command you need to run to execute the workflow, you may also need to build an Anaconda environment, compile GROMACS or locate a module that provides it, and use Singularity to run ffmpeg (this last item is only required for rendering videos). The `readme.md` instructions provide an explanation for how to supply the right software. 

## Data

We have three kinds of data which we may want to ingest.

1. Source data: the simulations. These data were previously uploaded to the compbio cluster at `/compbio/share/2021.02.15-banana-set`. This dataset includes 8 simulations with a naming convention `v1000` that include the so-called banana-shaped proteins.

2. Simulation slices. These are samples of the raw data which are written in the GROMACS formats for structures (gro) or trajectories (xtc). Simulation slices follow a naming format that includes the simulation name, the start, end, and sampling time, the group, and the PBC method. For example, the slice `v1034.2000000-12000000-10000.lipids.pbcmol.xtc` comes from the `v1034` simulations, from 2-12 microseconds at a sampling rate of 0.1 microsecond, including the lipids only, using molecular periodic boundary conditions (`pbcmol` in the GROMACS manual).

2. Post-processing data. This data type includes any data extracted from the simulations in a more usable format, most notably the bilayer midplane height profiles. They are stored with either an `hdf5` fomat (the `dat` suffix) and bundled alongside metadata files written in JSON format (with the `spec` suffix), for example `v1026.2000000-12000000-10000.lipid_abstractor.n0.spec`. You can see that the naming includes the simulation name and time range. We also number these files (e.g. `n0`) in case we run the calculation more than once, with different parameters.

Simulation slices and post-processing data were uploaded to `compbio:/compbio/share/2021.10.11-banana-post`. For the remainder of this guide, we will use the filename and data structure used by BioPhysCode for compatibility reasons (so we can use the existing data). The data are stored in `hdf5` files with a simple naming scheme.

## Environment

We will use the following environment to perform these calculations. The primary requirements are `scipy`, `h5py`, and `joblib`.

```
# install or load anaconda
ml anaconda
# create an environment
conda create bpcex 
conda activate bpcex
pip install scipy
pip install h5py
pip install joblib
```

## CODE A: plotting undulations

We will start by extracting the code to plot undulations from the postprocessing data (item 2 above). This code will have the following features:

1. The input data are the `undulations` post-processing data described above, e.g. `v1021.2000000-12000000-10000.undulations.n0.spec`.
2. The code will produce an undulation plot.

### Extract the code

The following block of code includes a list of the copy commands required to remove this code from the orginal pipeline. 

```
# go to the cluster
ssh rf
# change to the calc directory in the existing pipeline
cd /data/rbradley/legacy-factory/calc/banana
# make a new directory to hold the extracts
mkdir bpc-extracts
cd bpc-extracts
# next we will copy the codes we need from the calc/banana folder
# copy the undulations code available at https://github.com/bradleyrp/omni-single/blob/master/plot-undulations.py
cp ../calcs/plot-undulations.py .
# copy the supporting codes required by the import statements at the top of this script
# then we change the import statements to find these codes
# each of the following files is located in https://github.com/bradleyrp/omni-single/tree/master/codes
cp ../calcs/codes/undulate.py .
cp ../calcs/codes/undulate_plot.py .
cp ../calcs/codes/looptools.py .
# next I stripped the @autoplot decorators. these are not necessary
# remove the load function which we will replace
```

Almost all of the undulation code exists in the `https://github.com/bradleyrp/omni-single/` repository. To extract the `plot-undulations.py` so it could be used in a standalone fashion, I copied the files above, and then added some code to load the data (see below) and collected a few supporting functions from BioPhysCode.

### Loading the data

BioPhysCode includes a method for linking different calculations and data together. For example, we compute the undulation spectra from the lipid centers-of-mass, which are computed in a separate step. Since our goal is to write a stand-alone script to compute the undulation spectra, we need to adapt the code used to load the data. The original `plot-undulations.py` script has a few references to this data:

- `data,calc = plotload('undulations')`
- `sns = work.sns()`
- `work.meta` and `work.plotdir`

This code would typically search the `post` folder for pairs of `dat` and `spec` files which match the so-called "upstream" calculation, and then load these in a loop over the relevant simulations.

Instead of using the BioPhysCode method to manage the data, we simply select a list of dat files containing the regular grid of the lipids (we call these the "undulations" dat files) and send these to the script with the `-s` flag. The code then infers the simulation name from the file name (the text before the first `.` i.e. `v1021`) and then later links the undulations dat files to the `protein_abstractor` dat files if we are plotting the height profiles as well.

After adding code to find and load the right data, we are ready to run the calculation. 

### Plotting the undulation spectra

I have exposed the input directory, output directory, and simulation names to the command-line. The following commands will activate the environment and load the data.

```
# activate your environment to get scipy etc
ml anaconda
conda activate bpcex
# go to the code, wherever you cloned it
cd /data/rbradley/legacy-factory/calc/banana/bpc-extracts
# infer the input files from a glob
DAT_FILES=$(ls /data/rbradley/legacy-factory/data/banana/post/v1*undulations*n*dat)
# plot the spectra
python -i plot-undulations.py -m spectra -s $DAT_FILES -o ./fig.undulations.png
# plot the height profiles (simulation name is appended to the output name)
python -i plot-undulations.py -m height-profiles -s $DAT_FILES -o ./fig.height-profiles
```

I have recommended the interactive flag (`-i`) in case you want to inspect the data after the calculation is complete. You can also run the `go()` command in the interactive session to rerun the script without loading data.

In the block above, I have provided example commands for each of two modes (the `-m` flag). One mode plots the spectra for all simulations, and the other plots the height profiles iteratively. The above commands produce `fig.undulations.png` and e.g. `fig.height-profile-v1021.png`.

## CODE B: Generate height profiles from lipid positions

In "code A" above we used the existing undulations dat files to plot the undulation spectra and height profiles for our simulations. In this section, I will explain how the undulation data were generated from another piece of upstream data. The undulation calculation requires the following steps:

1. Run the simulation.
2. Sample the trajectory.
3. Sample the lipid centers of mass.
4. Interpolate the lipid centers of mass onto a regular grid (code B, this guide).
5. Use the regular grid to generate height profiles and undulation spectra (code A).

The code in this section will have the following features:

1. It will accept a single "lipid abstractor" file.
2. It will run the "undulations" calculation copied verbatim from BioPhysCode
3. It will output one of the "undulations" dat files used in code A.

### Extract the code 

To copy this code from BioPhysCode, I used the following commands.

```
# go to the bpc-extracts folder
cd /data/rbradley/legacy-factory/calc/banana/bpc-extracts
# copy the code from the repo at https://github.com/bradleyrp/omni-single/tree/master/codes
cp ../calcs/undulations.py calc-undulations.py
cp ../calcs/codes/mesh.py ./
cp ../omni/base/compute_loop.py .
cp ../omni/base/timer.py .
# each of the files above required minor edits to the import statements since the relative paths changed when removed them
```

In the BioPhysCode code, the `undulations` function was found inside the `undulations.py` module. The purpose of the YAML "metadata" file is to standardize the input and outputs for functions like  this one. The code inside the `undulations` function is the workhorse of the pipeline. This code calls a `makemesh_regular` function to generate a regular grid from the lipid positions.

In this exercise, we are removing the boilerplate so this calculation stands alone. The input for this calculation is a data structure which contains the box vectors, the lipid centers of mass (`points`), and the indices for the monolayers (`monolayer_indices`). Since we have removed the `undulations` function from BioPhysCode, I changed the input data. In the resulting script, `calc-undulations.py`, we load the lipid_abstractor dat file and then produce and undulation dat file.

### Run the calculation

After extracting `codes.undulations.undulations` (i.e. the `undulations` function in `undulations.py`) from the original BioPhysCode pipeline, we can now run it with one command.

```
python -i calc-undulations.py -i /data/rbradley/legacy-factory/data/banana/post/v1021.2000000-12000000-10000.lipid_abstractor.n0.dat -o /data/rbradley/legacy-factory/data/banana/post/v1021.2000000-12000000-10000.undulations.dat
```

This command produces an `undulations.dat` file which has no numbering scheme. In BioPhysCode, we add a number (e.g. `n0`) to match the data up with the parameters used to generate it.

### File organization

In the "code A" and "code B" examples above, I have extracted and repackaged the code which converts the lipid positions into a regular mesh which represents the bilayer height profiles (code B), then plots the undulation spectra and height profiles from these regular meshes (code A).

Since these codes still use the existing data structure in the `lipid_abstractor` dat files, they still depend on other (upstream) parts of BioPhysCode. Removing the codes from BioPhysCode also removes any file management features.

## CODE C: Generate lipid positions from a sampled trajectory

In this example, we start with a timewise sub-sampled trajectory of our simulation:

```
v1021.2000000-12000000-10000.lipids.pbcmol.gro
v1021.2000000-12000000-10000.lipids.pbcmol.xtc
```

We require both a structure (gro) and trajectory file (xtc) which are both prepared from the GROMACS simulation. In this example, we wish to generate the `lipid_abstractor` dat file required as an input for code B (above).

```
python -i calc-lipids.py -i /data/rbradley/legacy-factory/data/banana/post/v1021.2000000-12000000-10000.lipids.pbcmol -o /data/rbradley/legacy-factory/data/banana/post/v1021.2000000-12000000-10000.lipid_abstractor.dat
```

We supply the file prefix (without the `.gro` suffix) to the input flag on the `calc-lipids.py` script. This script uses MDAnalysis to read the trajectory, compute the lipid centers of mass, and then write the data to an h5 file which follows the format expected by code b.

## CODE D: Generate protein positions

To complete this calculation, we need to extract protein positions from the GROMACS trajectories to create a `protein_abstractor` data file. This file must match the name of the `undulations` dat file in order for code A (i.e. `plot-undulations.py` to find the right data and match it up to the height profiles (which come from the `lipid_abstractor` data. We ported this script into the `calc-proteins.py`. 

```
python -i calc-proteins.py -i /data/rbradley/legacy-factory/data/banana/post/v1021.2000000-12000000-10000.lipids.pbcmol -o /data/rbradley/legacy-factory/data/banana/post/v1021.2000000-12000000-10000.protein_abstractor.dat
```

## Recap

If you start from an xtc and gro file, you can run these scripts in sequence to generate the undulation spectrum and height profile (with protein positions):

1. `calc-lipids.py` (i.e. code C)
2. `calc-proteins.py` (i.e. code D)
3. `calc-undulations.py` (i.e. code B)
4. `plot-undulations.py` (i.e. code A)

Each of these scripts now has an argument parser. You can see the inputs and outputs with: `python calc-lipids.py -h`. Each simulation needs to run through all of the `calc-` scripts to produce the data required by the `plot-undulations.py` script.
