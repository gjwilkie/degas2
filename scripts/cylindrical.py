import dg2d
import defineback
import netCDF4 as nc
import numpy as np

# Outputs neutral density and temperature from netcdf output
# Assumes first len(rgrid) zones are the ones of interest
def process_output(rgrid,outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="tally.nc",plot=False,asciioutfilename=None):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)

    # Find the tally indicee to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    pres_idx = -1
    flux_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
        elif "neutral pressure" in str(nc.chartostring(tallynames[itally])):
            pres_idx = itally
        elif "neutral flux vector" in str(nc.chartostring(tallynames[itally])):
            flux_idx = itally

    # Find the indicee of the independent variables
    zone_idx = -1
    testsp_idx = -1
    sp_idx = -1
    rc_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar
        elif "test " in str(nc.chartostring(varnames[itally])):
            testsp_idx = ivar
        elif "problem_sp" in str(nc.chartostring(varnames[itally])):
            sp_idx = ivar
        elif "reaction" in str(nc.chartostring(varnames[itally])):
            rc_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    pres_base = tallydata["tally_base"][pres_idx]
    flux_base = tallydata["tally_base"][flux_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]
    NR=len(rgrid)
    Nspec = tally_indices[dens_idx,1]

    density = np.zeros((NR,Nspec-1))
    density_err = np.zeros((NR,Nspec-1))

    pressure = np.zeros((NR,Nspec-1))
    pressure_err = np.zeros((NR,Nspec-1))

    for izone in range(0,NR):
        for ispec in range(1,Nspec):
            density[izone,ispec-1] = outputdata["out_post_all"][dens_base+ispec*Nzone+izone,0]
            density_err[izone,ispec-1] = outputdata["out_post_all"][dens_base+ispec*Nzone+izone,1]
            pressure[izone,ispec-1] = outputdata["out_post_all"][pres_base+ispec*Nzone+izone,0]
            pressure_err[izone,ispec-1] = outputdata["out_post_all"][pres_base+ispec*Nzone+izone,1]

    if asciioutfilename:
        f = open(asciioutfilename,"w")
        f.write("#%14s %15s "%("izone","radius (m)"))
        for j in range(0,Nspec-1):
            f.write("%15s %15s "%("Density("+str(j+1)+")","rel. err."))
        for j in range(0,Nspec-1):
            f.write("%15s %15s "%("Pressure("+str(j+1)+")","rel. err."))
        f.write("\n")
        for i in range(0,NR):
            f.write("%15d %15e "%(i,rgrid[i]))
            for j in range(0,Nspec-1):
                f.write("%15e %15e "%(density[i,j],density_err[i,j]))
            for j in range(0,Nspec-1):
                f.write("%15e %15e "%(pressure[i,j],pressure_err[i,j]))
            f.write("\n")
        f.close()

    return density,pressure


def write_cylinder_input(R_tot,NR,source,ne,Te,Ti,material,Nflights=10000,walltemp=300.0,source_sp="H2",Rmin=1.0e-3,Zfac=10.0,t0=-1.0,tf=-1.0,init=False):
    rgrid, source_stratum = dg2d.write_cylinder_dg2d_input(R_tot,NR,material=material,walltemp=walltemp,Rmin=Rmin,Zfac=Zfac)

    defineback.write_cylinder_db_input(rgrid,ne,Te,Ti,source,Nflights,source_stratum,walltemp=walltemp,source_sp=source_sp,t0=t0,tf=tf,init=init)
    return rgrid

