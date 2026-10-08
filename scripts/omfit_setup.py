# Created on Nov. 21, 2024 by Quinn Pratt
# DEGAS2 setup inspired by 'micerscript.py' from A. Angulo and G. Wilkie
# 
# REQUIREMENTS
#   - DEGAS2 executables must be in $DEGAS2_BIN directory. 
#   - $DEGAS2_DIR/scripts must be added to $PYTHONPATH for certain imports below.
#
# PUBLIC DEGAS2 INSTALLs
#   - GA-OMEGA cluster,
#       >> module purge
#       >> module load degas2
#   - PPPL-FLUX cluster,
#       >> ?
#
# This script is designed to be called from a DEGAS2 work-directory (generally via OMFITx.executable)
# The directory must have: input_profiles.nc, geqdsk, tally.in, and degas2.in.
# 
# ----------------
# Core python,
import os
import glob
import subprocess
import sys
import json
from time import perf_counter
try:
    # From degas2/scripts,
    import dg2d
    import problem
    import source
    import defineback
    import postprocess
except ImportError:
    print("ERROR: (omfit_setup) Failed to import from Python API - is $DEGAS2_DIR/scripts in the PYTHONPATH?")
# Other,
import numpy as np
import scipy.interpolate as interpolate
import netCDF4 as nc
#import matplotlib.pyplot as plt

xid = "omfit_setup" # script ID.
verbose = True
all_programs = ['problemsetup', 'definegeometry2d', 'defineback', 'tallysetup']
performance_dict = dict(zip(all_programs, [None]*len(all_programs)))

d2path = os.environ.get('DEGAS2_BIN')
if d2path is None:
    print(f"ERROR: ({xid}) $DEGAS2_BIN not set! - load module or otherwise set env.")
    exit(1)
#d2path = os.environ["DEGAS2_BIN"] # set with >> module load degas2
print(f"INFO: ({xid}) DEGAS2_BIN={d2path}")
# --------------------------------
# Local lib of helper functions,

def check_for_file(filename, fail=False):
    """ Helper function, used throughout"""
    file_exists = os.path.exists(filename)
    if file_exists:
        print(f"INFO: ({xid}) found '{filename}'.")
    else:
        if fail:
                print(f"ERROR: ({xid}) '{filename}' not found - ending.")
                exit(1)
        else:
            print(f"WARN: ({xid}) '{filename}' not found.")
    return file_exists

def run_program(program, inputs=None, outputs=None, timeout_min=5):
    """ Helper wrapper for running degas2 programs with subprocess.run()"""
    spkw = dict(            # kwargs passed to subprocess.run() throughout.
        shell=False,        # generally recommended to avoid shell injection (although unlikely in this case)
        check=True,         # capture errors.
        timeout=int(timeout_min*60),
        )
    if program not in all_programs:
        raise ValueError(f"ERROR: ({xid}) 'run_program' wrapper encountered unknown program='{program}'.")
    print(f"EXEC: ({xid}) Running {program}...")
    command = [d2path+f"/{program}"]
    if inputs is None:
        inputs = []
    else:
        if not isinstance(inputs, list):
            printe(f"ERROR: ({xid}) inputs kwarg to run_program() must be a list of str.")
            exit(1)
    command += inputs # append args to list.

    try:
        t0 = perf_counter()
        subprocess.run(command, **spkw)
        performance_dict[program] = perf_counter() - t0 # capture runtime.
    except subprocess.CalledProcessError as e:
        # This block triggers automatically if the command returns a non-zero exit code
        message = f"""ERROR: DEGAS2 {program} failed!
- Command: {' '.join(e.cmd)}
- Exit code: {e.returncode}"""
        print(message)
        exit(1)
    except FileNotFoundError:
        # Triggers if the program/executable itself cannot be found on the system path
        print(f"ERROR: The executable '{program}' was not found.")
        exit(1)
    except subprocess.TimeoutExpired as e:
        print(f"ERROR: The executable '{program}' timed-out after {e.timeout} seconds.")
        exit(1)
    # if we made it to here, the program finished.
    # check for outputs if given before declaring success.
    if outputs is None:
        success = True
    else:
        output_check = []
        for f in outputs: 
            output_check += [ check_for_file(f) ]
        success = all(output_check)
    if success:
        print(f"EXEC: ({xid}) Completed {program}.\n")
    else:        
        print(f"WARN: ({xid}) Completed {program}, but some output files were not found...\n")
    return success

def print_performance():
    print("="*32)
    print(f"INFO: ({xid}) DEGAS2 performance report:")
    print("Timing...")
    for prog, time in performance_dict.items():
        if isinstance(time, float):
            print(f"\t- {prog} = {time:.3f} [sec]")
        else:
            print(f"\t- {prog} skipped")
    
    print("Memory use...")
    total = 0
    for ext in ['in','txt','nc']:
        # Filter files in the current directory only
        files = [f for f in glob.glob(f"*.{ext}") if os.path.isfile(f)]
        # Sum the sizes and convert to MB
        size = sum(os.path.getsize(f) for f in files) / (1024 * 1024)
        print(f"\t- size of *.{ext} = {size:.3f} [MB]")
        total += size
    print(f"Total = {total:.3f} [MB]")
    print("="*32)

# --------------------------------
# General,
# - input files (for this script),
profile_fname = "input_profiles.nc"
geqdsk_fname = "geqdsk"
setup_kwargs_fname = "omfit_setup_dict.json" # technically optional

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

# --------------------------------
# DEGAS2 problemsetup,
# From the degas2/scripts/problem.py
run_problemsetup = setup_kwargs.get("run_problemsetup", True)
if run_problemsetup:
    problemsetup_kwargs = setup_kwargs.get("problemsetup", {})
    # Toggle between 'genStdProblem' and the underlying 'generateProblemInput',
    custom_problem_input = problemsetup_kwargs.get("custom_problem_input", False)
    # 
    if custom_problem_input:
        print(f"INFO: ({xid}) Using custom 'problemsetup' inputs - verify values in 'problem.in'")
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
        print(f"INFO: ({xid}) Using standard problem label='{std_problem_label}'")
        p = problem.genStdProblem(std_problem_label)
    
    # >>>>
    # outputs = ['problem.in','problem.nc']
    run_program("problemsetup", outputs=['problem.nc'])
    # <<<<
else:
    # We can skip this if a 'problem.nc' file is provided,
    print(f"INFO: ({xid}) Skipping problemsetup...\n")
    check_for_file("problem.nc", fail=True)

# Check if recombination is present in the problem by parsing the problem.nc file.
# This will influence the output of defineback later. 
# i.e. if we need to edit the recombination source_n_flights.
rxn_dict = problem.get_reactions_from_problem("degas2.in", "problem.nc")
rxn_names = list(rxn_dict.keys())
recomb_included = any(["recomb" in s for s in rxn_names])
if recomb_included:
    print(f"INFO: ({xid}) Recombination reaction found in problem.nc - will check for 'recomb_n_flights' param.")

# --------------------------------
# DEGAS2 definegeometry2d,
# WARN: updated workflow expects a "limiterfile.txt" from OMFIT.
# TODO: 
#       - option to write basic limiterfile.txt in degas2/scripts.
#       - move code to read limiterfile.txt to degas2/scritpts.
#       - carefully check cases where custom_tri=True but no limiter found (!)     
run_definegeometry2d = setup_kwargs.get("run_definegeometry2d", True)
if run_definegeometry2d:
    # 1. Unpack dg2d settings,
    dg2d_kwargs = setup_kwargs.get("definegeometry2d", {}) # dict of options for dg2d scripts.

    # Whether or not we've been provided with a custom mesh in the form of Triangle .node/.ele files.
    custom_tri = dg2d_kwargs.get("custom_tri", False)

    # 2. Read limiter file, TO DO: make DEGAS2/scripts/utils.py method to handle this.
    limiterfile_basename = dg2d_kwargs.get("limiterfile_basename", "limiterfile")
    limiterfile_fname = limiterfile_basename + ".txt"
    check_for_file(limiterfile_fname, fail=True)
    seg_inds, Rlim, Zlim, material, walltemp, recyc, exitzone = np.genfromtxt(limiterfile_fname,
        delimiter=' ',dtype=object,encoding='utf-8').T
    # format arrays...
    seg_inds, exitzone = seg_inds.astype(int), exitzone.astype(int)
    Rlim, Zlim = Rlim.astype(float), Zlim.astype(float)
    material = material.astype("U8") # up to 8 chars.
    walltemp, recyc = walltemp.astype(float), recyc.astype(float)

    # 4. Call the *new* dg2d api workflow,
    dg2d_obj = dg2d.DG2D()
    dg2d_obj.polygonfile_name = "polygons.nc" # backward compatibility, now defaults to singular 'polygon.nc'
    if custom_tri:
        tri_basename = dg2d_kwargs.get("custom_tri_basename", "flux_surfaces")
        # check for necessary files,
        for ext in [".ele", ".node"]:
            check_for_file(tri_basename+ext, fail=True)
        # extract nodes and connectivity,
        coords, conn, wallnode_ids = dg2d_obj.load_from_triangle(tri_basename)
        # define mesh instead of defining limiter (see below),
        dg2d_obj.define_mesh(coords, conn, boundary_nodes=wallnode_ids, Rb=Rlim, Zb=Zlim)
    else:
        dg2d_obj.define_limiter(Rlim,Zlim)
    dg2d_obj.set_wallprops() # preallocates wall properties as arrays.
    dg2d_obj.set_wallprops(walltemp=walltemp, 
        material=material, 
        Rcoeff=recyc, 
        exitzone=exitzone, 
        wallidx=seg_inds) # sets all properties from seg_inds indexing.
    dg2d_obj.write_files(aux_thickness=0.005) # <-- this will write 'dg2d.in', and 'wallfile.txt'

    # >>>>
    # outputs = ["geometry.nc","geomtestc.silo","polygons.nc"]
    run_program("definegeometry2d", inputs=['dg2d.in'], outputs=['geometry.nc'])
    # <<<<
else:
    # We can skip this if a 'geometry.nc' file is provided.
    # I don't think we actually need the 'polygons.nc' file.
    print(f"INFO: ({xid}) Skipping definegeometry2d...\n")
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
        defineback.generate_plasma_file_through_psi(ne_data, Te_data, Ti_data, psifunc, **back_kw)
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
    source.write_db_input(sourcegroups) # outputs = ['db.in']
    # >>>>
    # outputs = ['background.nc','density*.txt','temperature*.txt']
    run_program("defineback", inputs=['db.in'], outputs=['background.nc'])
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
    print(f"INFO: ({xid}) Skipping defineback...\n")
    check_for_file("background.nc", fail=True)

# ----------------
# DEGAS2 tallysetup,
run_tallysetup = setup_kwargs.get("run_tallysetup", True)
if run_tallysetup:
    # Use supplied tally.in file
    # >>>>
    # outputs = ['tally.nc']
    run_program("tallysetup", outputs=['tally.nc'])
    # <<<<
else:
    # We can skip this if a 'tally.nc' file is provided.
    print(f"INFO: ({xid}) Skipping tallysetup...\n")
    check_for_file("tally.nc", fail=True)

if verbose:
    print_performance()

print(f"DONE: ({xid})")
