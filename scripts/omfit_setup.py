# Created on Nov. 21, 2024 by qpratt
# Example degas2 run based on the micerscript.py from  A. Angulo and G. Wilkie
# 
# SETUP: This script calls degas2 executables which should be in the $DEGAS2_BIN dir. 
#        Also, make sure degas2/scripts is added to $PYTHONPATH for the python modules below.
# The user should run the following commands before executing this scipt,
# >> module purge
# >> module load degas2/1.0/run_gcc8.x
# This will set up the necessary env. vars and modify the path.

# Core python,
import os
import subprocess
import sys
import json
# From degas2/scripts,
import dg2d
import problem
import source
import defineback
import postprocess
# Other,
import numpy as np
import scipy.interpolate as interpolate
import netCDF4 as nc
import matplotlib.pyplot as plt

# ----------------
# Macroscopic things for this script,
# - input filenames (for this script),
profile_fname = "input_profiles.nc"
geqdsk_fname = "geqdsk"
setup_kwargs_fname = "omfit_setup_dict.json" # optional

def check_for_file(filename, fail=False):
    """ Helper function, used throughout"""
    file_exists = os.path.exists(filename)
    if file_exists:
        print(f"INFO (omfit_setup): found '{filename}'.")
    else:
        if fail:
                print(f"ERROR (omfit_setup): '{filename}' not found - ending.")
                exit(1)
        else:
            print(f"WARN (omfit_setup): '{filename}' not found.")
    return file_exists

# Check for necessary files,
necessary_files = [profile_fname, geqdsk_fname, 'tally.in', 'degas2.in']
for f in necessary_files:
    check_for_file(f, fail=True)
# Check for optional settings file,
setup_kwargs_exists = check_for_file(setup_kwargs_fname)
if setup_kwargs_exists:
    with open(setup_kwargs_fname) as f:
        setup_kwargs = json.load(f)
else:
    setup_kwargs = dict()

# ----------------
# Macroscopic DEGAS2 setup,
d2path = os.environ["DEGAS2_BIN"] # set with module load degas2/1.0/run_gcc8.x.lua
# start with inputs = ['degas2.in','tally.in']
print(f"INFO (omfit_setup): d2path={d2path}")

# ----------------
# Profile setup,
# This script follows the convention that psi is normalized (0 = magnetic axis, 1 = separatrix)
# Kinetic profiles are read from a .nc dataset,
profiles = nc.Dataset(profile_fname)
profiles.set_auto_mask(False) # makes variables come in as np.array rather than MaskedArrays
# Append something far away (rho = 2.0) so things are defined over the whole grid.
psi_data = list(profiles["rho"][:]) + [2.0]
ne_data = list(profiles["n_e"][:]) + [0.] # [/m3]
Te_data = list(profiles["T_e"][:]) + [0.] # [eV]
Ti_data = list(profiles["T_i"][:])  + [0.]# [eV]
# Toroidal angular rotation data, V_tor = R*omega
omega_data = list(profiles["omega"][:]) + [0.] # [rad/s]
profiles.close()

# ----------------
# Magnetic equilibrium,
# Must be put in working directory, 
geqdsk_file = geqdsk_fname

# ----------------
# DEGAS2 problemsetup,
run_problemsetup = setup_kwargs.get("run_problemsetup", True)
if run_problemsetup:
    # From the degas2/scripts/problem.py
    # Generates degas2 problem input files based on a predefined 'C-D' case.
    p = problem.genStdProblem("C-D")
    # >>>>
    # outputs = ['problem.in','problem.nc']
    print("Running problemsetup...\n")
    subprocess.run(d2path+"/problemsetup",shell=True)
    # <<<<
else:
    # We can skip this if a 'problem.nc' file is provided,
    print("Skipping problemsetup...\n")
    check_for_file("problem.nc", fail=True)

# ----------------
# DEGAS2 definegeometry2d,
# NOTE: we still need to run 'generateGeometryFromEFITfile'
#       to create the psifunc interpolant for defineback, even if
#       we don't run the actual definegeometry2d executable.
run_definegeometry2d = setup_kwargs.get("run_definegeometry2d", True)
dg2d_kwargs = setup_kwargs.get("definegeometry2d", {})

# Recycling coefficient
recyc = dg2d_kwargs.get("recyc", 0.98)
# Wall Material
material = dg2d_kwargs.get("material","C")
# Wall temperature in Kelvin
walltemp = dg2d_kwargs.get("walltemp",300.0)
# When refining mesh, this is the largest segment permitted along the wall [m]
# max distance along limiter for triangulation. Roughly sets the spatial res.
dlim_max = dg2d_kwargs.get("dlim_max", 0.02)

# From the degas2/scripts/dg2d.py
geo_kw = dict(recyc_coef=recyc, Twall=walltemp,dlim_max=dlim_max,
              clockwise=True, # Not sure why this is True, default is False.
              RZlim_start=[1.0128, 1.2053], # Upper HFS point on limiter.
         )
# [AA] George will change this to generate rho instead of psifunc
# This function uses the gEQDSK file to generate the \psi_n(R, Z) interpolant (psifunc).
# outputs = ['wallfile.txt','dg2d.in']
psifunc, nodes = dg2d.generateGeometryFromEFITfile(geqdsk_file, material, **geo_kw)
if run_definegeometry2d:
    # >>>>
    # outputs = ["geometry.nc","geomtestc.silo","polygons.nc"]
    print("Running definegeometry2d...\n")
    subprocess.run(d2path+"/definegeometry2d dg2d.in",shell=True)
    # <<<<
else:
    # We can skip this if a 'geometry.nc' file is provided.
    # I don't think we actually need the 'polygons.nc' file.
    print("Skipping definegeometry2d...\n")
    check_for_file("geometry.nc", fail=True)

# ----------------
# DEGAS2 defineback,
run_defineback = setup_kwargs.get("run_defineback", True)
if run_defineback:
    # From the degas2/scripts/defineback.py script,
    # Defines the background plasma,
    # outputs = ['plasmafile.txt']
    back_kw = dict(psi_data=psi_data, rot_data=omega_data)
    defineback.generate_plasma_file_through_psi(ne_data,Te_data,Ti_data, psifunc,**back_kw)

    # From the degas2/scripts/source.py script,
    # Source Option #1: Puff,
    #     This will produce a "puff" type source at a temperature of 300K uniformly around limiter at a strength of
    #     1.0e24 nuclei per m^2 per s. Using Nsample flights and treating as a 
    #sgroup = source.Source(Nsample,"puff","D",rootspecies="D",pufftemp=300.0,strength=1.0e24,stratum=3,segment="112", specify_flux=True)

    # Source Option #2: Plate,
    #     This treats the source as if it's recycling, same flux as above, energy distribution of produced neturals
    #     are determined by ions near the wall and the recycling properies of the PFC.
    
    db_kwargs = setup_kwargs.get("defineback", {})
    # The number of flight samples
    Nsample = db_kwargs.get("Nsample", 1000)
    source_strength = db_kwargs.get("source_strength", 1.e24) # [m2/s]
    stratum = db_kwargs.get("stratum", 3)   
    segment = db_kwargs.get("segment","*")

    source_kw = dict(rootspecies="D+", #  
                     strength=source_strength, # [m2/s]
                     stratum=stratum, # 
                     segment=segment, #
                    )
    # outputs = ['sourcefile.txt'] CURRENTLY ONLY ONE GROUP IS SUPPORTED,
    #   may only exist for large numbers of source segments.
    sgroup = source.Source(Nsample,"plate","D",**source_kw)
    # outputs = ['db.in']
    source.write_db_input([sgroup])
    # >>>>
    # outputs = ['background.nc','density*.txt','temperature*.txt']
    print("Running defineback...\n")
    subprocess.run(d2path+"/defineback db.in",shell=True)
    # <<<<
else:
    # We can skip this if a 'background.nc' file is provided.
    print("Skipping defineback...\n")
    check_for_file("background.nc", fail=True)

# ----------------
# DEGAS2 tallysetup,
run_tallysetup = setup_kwargs.get("run_tallysetup", True)
if run_tallysetup:
    # Use supplied tally.in file
    # >>>>
    # outputs = ['tally.nc']
    print("Running tallysetup...\n")
    subprocess.run(d2path+"/tallysetup",shell=True)
    # <<<<
else:
    # We can skip this if a 'tally.nc' file is provided.
    print("Skipping tallysetup...\n")
    check_for_file("tally.nc", fail=True)

print("DONE (omfit_setup)")
