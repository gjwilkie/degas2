# Created on Nov. 21, 2024 by Quinn Pratt
# DEGAS2 setup based on the micerscript.py from A. Angulo and G. Wilkie
# 
# SETUP: This script calls degas2 executables which should be in the $DEGAS2_BIN directory. 
#        Make sure degas2/scripts is added to $PYTHONPATH for the python modules below.
#
# The user should run the following commands (on the omega cluster at GA) before executing this scipt,
# >> module purge
# >> module load degas2
# This will set up the necessary env. vars and modify the path/pythonpath.
# ----------------
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
#import matplotlib.pyplot as plt

# ----------------
# General,
# - input files (for this script),
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
d2path = os.environ["DEGAS2_BIN"] # set with >> module load degas2
# start with inputs = ['degas2.in','tally.in']
print(f"INFO (omfit_setup): DEGAS2_BIN={d2path}")

# ----------------
# Magnetic equilibrium,
geqdsk_file = geqdsk_fname
from geomutils import read_geqdsk
g = read_geqdsk(geqdsk_file)
rgrid = g.rgrid
zgrid = g.zgrid
# Defines psi_n(R,Z) = (psi - min(psi))/(max(psi) - min(psi)), Normalized poloidal flux.
psin_rz = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag)
psifunc = lambda R,Z : interpolate.RectBivariateSpline(rgrid, zgrid, psin_rz.T)(R, Z)[0]

# ----------------
# Profile setup,
# We use the normalized poloidal magnetic flux, psi_n, as the radial coordinate.
# Kinetic profiles are read from a .nc dataset,
profiles = nc.Dataset(profile_fname)
profiles.set_auto_mask(False) # makes variables come in as np.array rather than MaskedArrays

# clip temperatures, 
min_temp = 300./11650 # 300 K --> eV
T_i = np.clip(profiles["T_i"][:],a_max=None, a_min=min_temp)
T_e = np.clip(profiles["T_e"][:],a_max=None, a_min=min_temp)
# clip densities,
min_dens = 0.0
n_e = np.clip(profiles["n_e"][:],a_max=None, a_min=min_dens)
n_i = np.clip(profiles["n_i"][:],a_max=None, a_min=min_dens)

# Append 0 at the grid boundary so the plasma is defined over the whole grid.
psi_data = list(profiles["psi_n"][:]) + [np.amax(psin_rz)]
ne_data = list(n_e) + [0.] # [/m3]
ni_data = list(n_i) + [0.] # [/m3]
Te_data = list(T_e) + [0.] # [eV]
Ti_data = list(T_i) + [0.]# [eV]
# Toroidal angular rotation data, V_tor = R*omega
omega_data = list(profiles["omega"][:]) + [0.] # [rad/s]
profiles.close()

# ----------------
# DEGAS2 problemsetup,
# From the degas2/scripts/problem.py
run_problemsetup = setup_kwargs.get("run_problemsetup", True)
if run_problemsetup:
    problemsetup_kwargs = setup_kwargs.get("problemsetup", {})
    # Toggle between 'genStdProblem' and the underlying 'generateProblemInput',
    custom_problem_input = problemsetup_kwargs.get("custom_problem_input", False)
    # 
    if custom_problem_input:
        print(f"INFO (omfit_setup): Using custom 'problemsetup' inputs - check values in 'problem.in'")
        # Gather args for 'generateProblemInput' from problemsetup_kwargs, use "C-D" as default,
        cd_problem_defaults = dict(test_species=["0","D","D2","D2+"],
                               background_species=["e","D+"],
                               reactions=["hionize5","dd_chargex","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec"],
                               materials=["C"],
                               pmis=["hdesorbc","h2desorbc","dreflc"],
                           )
        args = []
        for k, v in cd_problem_defaults.items():
            args += [problemsetup_kwargs.get(k, v)]
        # pass args to 'generateProblemInput()'
        p = problem.generateProblemInput(*args)
    else:
        # Generates degas2 problem input files based on a predefined cases.
        std_problem_label = problemsetup_kwargs.get("std_problem_label","C-D")
        print(f"INFO (omfit_setup): Using standard problem label='{std_problem_label}'")
        p = problem.genStdProblem(std_problem_label)
    # >>>>
    print("Running problemsetup...\n")
    subprocess.run(d2path+"/problemsetup",shell=True)
    # outputs = ['problem.in','problem.nc']
    # <<<<
else:
    # We can skip this if a 'problem.nc' file is provided,
    print("Skipping problemsetup...\n")
    check_for_file("problem.nc", fail=True)

# Check if recombination is present in the problem by parsing the problem.nc file.
# This will influence the output of defineback later. 
# i.e. if we need to edit the recombination source_n_flights.
rxn_dict = problem.get_reactions_from_problem("degas2.in", "problem.nc")
rxn_names = list(rxn_dict.keys())
recomb_included = any(["recomb" in s for s in rxn_names])

# ----------------
# DEGAS2 definegeometry2d,
# NOTE: defineback needs the psifunc interpolant to map the 1D profiles onto the
#       psifunc is now defined above from the gEQDSK file.
run_definegeometry2d = setup_kwargs.get("run_definegeometry2d", True)
if run_definegeometry2d:
    dg2d_kwargs = setup_kwargs.get("definegeometry2d", {}) # dict of options for dg2d scripts.
    # Whether or not we've been provided with a custom mesh, 
    custom_tri = dg2d_kwargs.get("custom_tri", False)
    # Recycling coefficient
    recyc = dg2d_kwargs.get("recyc", 0.98)
    # Wall Material
    material = dg2d_kwargs.get("material","C")
    # Wall temperature in Kelvin
    walltemp = dg2d_kwargs.get("walltemp",300.0)
    # When refining mesh, this is the largest segment permitted along the wall [m]
    # max distance along limiter for triangulation. Roughly sets the spatial res.
    dlim_max = dg2d_kwargs.get("dlim_max", 0.02)
    minarea = dg2d_kwargs.get("minarea",-1)
    # This is an optional point on/near the limiter to index as 0.
    # This can assist with defining distributed sources using the "start:end" method.
    # Set to None to disable.
    # Set to (1.0128, 1.2053) [m] for the upper HFS edge.
    RZlim_start = dg2d_kwargs.get("RZlim_start", None)

    # From the degas2/scripts/dg2d.py
    geo_kw = dict(recyc_coef=recyc, 
                  Twall=walltemp,
                  dlim_max=dlim_max,
                  clockwise=True, # Not sure why this is True, default is False.
                  RZlim_start=RZlim_start, # Upper HFS point on limiter.
                  minarea=minarea,
             )
    if custom_tri:
        tri_basename = dg2d_kwargs.get("custom_tri_basename", "flux_surfaces")
        # check for necessary files,
        for ext in [".ele", ".node"]:
            check_for_file(tri_basename+ext, fail=True)
        # generate the geometry files,
        # this script will write, 
        # 1. the dg2d.in file
        # 2. the wallfile. 
        from xgcpost import write_geometry_files
        write_geometry_files(material=material,
            recyc=recyc,
            walltemp=walltemp,
            use_xgc_mesh=True,
            polygonfilename="polygons.nc",
            trifile_base=tri_basename,
            make_plot=False,
        )
        # NOTE: psifunc is still defined from above.
    else:
        # Use the new DG2D class to deal with cases direct from a gEQDSK file...
        # dg2d.setup(material, recyc, gfile=geqdsk_file, walltemp=walltemp, run_dg2d=False)
	    # OLD...
        # This function uses the gEQDSK file to generate the \psi_n(R, Z) interpolant (psifunc).
        # the psifunc is used later when we run defineback.
        psifunc, nodes = dg2d.generateGeometryFromEFITfile(geqdsk_file, material, **geo_kw)
        # outputs = ['wallfile.txt','dg2d.in']
    # >>>>
    print("Running definegeometry2d...\n")
    subprocess.run(d2path+"/definegeometry2d dg2d.in",shell=True)
    # outputs = ["geometry.nc","geomtestc.silo","polygons.nc"]
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
    db_kwargs = setup_kwargs.get("defineback", {})
    
    # Define the background plasma,
    # From the degas2/scripts/defineback.py script,
    custom_plasmafile = db_kwargs.get("custom_plasmafile",False)
    if custom_plasmafile:
        print(f"INFO (omfit_setup): Using custom 'plasmafile.txt.")
        check_for_file("plasmafile.txt", fail=True)
    else:
        back_kw = dict(psi_data=psi_data, rot_data=omega_data, ni_data=ni_data)
        defineback.generate_plasma_file_through_psi(ne_data, Te_data, Ti_data, psifunc,**back_kw)
        # outputs = ['plasmafile.txt']

    n_sourcegroups = db_kwargs.get("n_sourcegroups", 1)
    sourcegroups = []
    for i in range(n_sourcegroups):
        sg_kwargs = db_kwargs.get(f"sourcegroup_{i}", {})
        
        # Unpack args for the source.Source class,
        n_flights = sg_kwargs.get("n_flights", 1000)
        source_type = sg_kwargs.get("source_type",'plate')
        source_species = sg_kwargs.get("source_species", "D")
        source_root_species = sg_kwargs.get("source_rootspecies", source_species+"+")
        if source_type in ['vol_source','puff']:
            source_root_species = source_species

        if source_type == "vol_source":
            # Force certain settings,
            sg_kwargs['custom_sourcefile'] = True # must use sourcefile, forces the check for sourcefile.txt later.
            
        # Unpack kwargs for the source.Source class,
        source_kw = dict(rootspecies=source_root_species,
                         specify_flux=sg_kwargs.get("specify_flux", True),
                         pufftemp=sg_kwargs.get("puff_temp",None),
                         puffexp=sg_kwargs.get("puff_exp", None),
                         strength=sg_kwargs.get("source_strength", 1.e23), # [/m2/s] when specify_flux in 'db.in'
                         stratum=sg_kwargs.get("stratum", 3),
                         segment=sg_kwargs.get("segment", "*"),
                         sourcefile=f"sourcefile_{i}.txt", 
                         )
        
        # Option to pass a sourcefile.txt directly...
        custom_sourcefile = sg_kwargs.get("custom_sourcefile",False)
        if custom_sourcefile:
            print(f"INFO (omfit_setup): Source Group {i} requires custom 'sourcefile_{i}.txt.")
            check_for_file(f"sourcefile_{i}.txt", fail=True)
            source_kw["strength"] = None # disables the source-by-segment methods in the Source class.
            
        # outputs = ['sourcefile_{i}.txt']
        #   may only exist for large numbers of source segments.
        sourcegroups += [ source.Source(n_flights,source_type,source_species,**source_kw) ]
    source.write_db_input(sourcegroups)
    # outputs = ['db.in']
    # >>>>
    print("Running defineback...\n")
    subprocess.run(d2path+"/defineback db.in",shell=True)
    # outputs = ['background.nc','density*.txt','temperature*.txt']
    # <<<<

    if recomb_included:
        # If a recombination reaction (e.g. 'hrecombine5') is included in the problem, 
        # DEGAS2 will automatically add an additional (final) source-group when defineback is executed above.
        # However, DEGAS2 only uses 100 flights for this source, here we update the 'background.nc' directly
        # to update the num_flights used for recombination,
        recomb_n_flights = db_kwargs.get("recomb_n_flights", 1000)
        print(f"INFO: (omfit_setup) Recombination included in problem, will set n_flights={recomb_n_flights} in background.nc")
        with nc.Dataset("background.nc", 'r+') as ds:
            ds.variables["source_num_flights"][-1] = recomb_n_flights
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
    print("Running tallysetup...\n")
    subprocess.run(d2path+"/tallysetup",shell=True)
    # outputs = ['tally.nc']
    # <<<<
else:
    # We can skip this if a 'tally.nc' file is provided.
    print("Skipping tallysetup...\n")
    check_for_file("tally.nc", fail=True)

print("DONE (omfit_setup)")
