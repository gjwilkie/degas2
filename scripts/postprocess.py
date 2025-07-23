import matplotlib.pyplot as plt
import numpy as np
import sys
import netCDF4 as nc
import matplotlib.tri as tri
from importlib import reload
from matplotlib import ticker
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator


def get_output(tallyname,sgroup=None,outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc",debug=False,with_err=False,Nzone=None):
    Nvar_max = 5

    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)
    geomdata = nc.Dataset(geometryfilename)

    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)

    tally_idx = -1
    found=False
    for itally in range(0,Ntally):
        if debug:
            print(str(nc.chartostring(tallynames[itally])).strip())
        if tallyname.strip() == str(nc.chartostring(tallynames[itally])).strip():
            found = True
            tally_idx = itally
    if not found:
        print("ERROR: Could not find tally named "+tallyname+"\n")

    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)

    base_idx = tallydata["tally_base"][tally_idx]

    tally_indices = np.array(tallydata["tally_tab_index"][tally_idx][:])

    if not Nzone == None:
        tally_indep_var = tallydata["tally_indep_var"][tally_idx][:]
        if tally_indep_var[0] == 1:
            tally_indices[0]=Nzone
        if tally_indep_var[1] == 1:
            tally_indices[1]=Nzone
        if tally_indep_var[2] == 1:
            tally_indices[2]=Nzone
        if tally_indep_var[3] == 1:
            tally_indices[3]=Nzone
        if tally_indep_var[4] == 1:
            tally_indices[4]=Nzone
    Ndat_tally = np.prod(tally_indices)

    Ngroup = 1
    if (sgroup == None):
        outdata_raw = np.array(outputdata["out_post_all"][base_idx:base_idx+Ndat_tally,0])
        error_raw = np.array(outputdata["out_post_all"][base_idx:base_idx+Ndat_tally,1])
    elif (sgroup == "all" or sgroup == "ALL" or sgroup == "All"):
        outdata_raw = np.array(outputdata["out_post_grp"][:,base_idx:base_idx+Ndat_tally,0])
        error_raw = np.array(outputdata["out_post_grp"][:,base_idx:base_idx+Ndat_tally,1])
        Ngroup = len(outdata_raw[:,0])
    else:
        outdata_raw = np.array(outputdata["out_post_grp"][sgroup,base_idx:base_idx+Ndat_tally,0])
        error_raw = np.array(outputdata["out_post_grp"][sgroup,base_idx:base_idx+Ndat_tally,1])
    outdata = np.squeeze(np.reshape(outdata_raw,np.append(Ngroup,tally_indices),order='F'))
    error = np.squeeze(np.reshape(error_raw,np.append(Ngroup,tally_indices),order='F'))

    outputdata.close()
    geomdata.close()
    tallydata.close()

    if with_err:
        return outdata, error
    else:
        return outdata

    

# Returns neutral density and detector signals 
def process_output_with_llama(outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc",plot=False,lymandatafile=None):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)

    geomdata = nc.Dataset(geometryfilename)
    zone_coords_3D = geomdata["zone_center"]
    zone_type = geomdata["zone_type"]
    x_zone = zone_coords_3D[:][0]
    z_zone = zone_coords_3D[:][2]
    zone_volumes = geomdata["zone_volume"]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    emission_idx = -1
    signal_idx = -1
    ioniz_idx = -1
    esource_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
#        elif "Lyman detector view" in str(nc.chartostring(tallynames[itally])):
        elif "LLAMA detector view" in str(nc.chartostring(tallynames[itally])):
            signal_idx = itally
        elif "Lyman emission rate" in str(nc.chartostring(tallynames[itally])):
            emission_idx = itally
        elif "ion source rate" in str(nc.chartostring(tallynames[itally])) and not "by reaction" in str(nc.chartostring(tallynames[itally]) ) and not "total" in str(nc.chartostring(tallynames[itally])):
            ioniz_idx = itally
        elif "ion energy source" in str(nc.chartostring(tallynames[itally])) and not "by reaction" in str(nc.chartostring(tallynames[itally])):
            esource_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    reac_idx = -1
    sp_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar
        elif "detector " in str(nc.chartostring(varnames[itally])):
            detector_idx = ivar
        elif "reaction " in str(nc.chartostring(varnames[itally])):
            reac_idx = ivar
        elif "problem_sp " in str(nc.chartostring(varnames[itally])):
            sp_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    signal_base = tallydata["tally_base"][signal_idx]
    emission_base = tallydata["tally_base"][emission_idx]
    ioniz_base = tallydata["tally_base"][ioniz_idx]
    esource_base = tallydata["tally_base"][esource_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    Ndetector = tally_indices[signal_idx,0]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    x = []
    z = []
    density = []
    density_err = []
    emission = []
    ioniz = []
    esource = []
    vols = []
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
           x.append(geomdata["zone_center"][izone,0])
           z.append(geomdata["zone_center"][izone,2])
           density.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,0])
           density_err.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,1])
           emission.append(outputdata["out_post_all"][emission_base+izone,0])
           ioniz.append(outputdata["out_post_all"][ioniz_base+izone,0])
           esource.append(outputdata["out_post_all"][esource_base+3*Nzone+izone,0])
           vols.append(geomdata["zone_volume"][izone])

    signal = []
    signal_err = []
    for idet in range(0,Ndetector):
        signal.append(outputdata["out_post_all"][signal_base+idet,0])
        signal_err.append(outputdata["out_post_all"][signal_base+idet,1])

    x = np.array(x)
    z = np.array(z)
    density = np.array(density)
    density_err = np.array(density_err)
    emission = np.array(emission)
    signal = np.array(signal)
    signal_err = np.array(signal_err)
    ioniz = np.array(ioniz)
    esource = np.array(esource)

    if plot:
        triang = tri.Triangulation(x,z)
        plt.gca().set_aspect("equal")
        plt.title("Neutral density (m^-3)")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,np.maximum(1e10,density),locator=ticker.LogLocator())
        plt.colorbar()
        plt.savefig("density.pdf",bbox_inches="tight")
        plt.close()

#        triang = tri.Triangulation(x,z)
#        plt.title("Neutral density (m^-3)")
#        plt.xlabel("x (m)")
#        plt.ylabel("z (m)")
#        plt.tricontourf(triang,np.log(density))
#        plt.colorbar()
#        plt.savefig("logdensity.pdf",bbox_inches="tight")
#        plt.close()


        plt.title("Relative error of neutral density")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,density_err)
        plt.colorbar()
        plt.savefig("error.pdf",bbox_inches="tight")
        plt.close()

        plt.title("Lyman-alpha emission (W / m^2)")
        plt.gca().set_aspect("equal")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,emission)
        plt.gca().set_aspect("equal")
        plt.colorbar()
        plt.savefig("emission.pdf",bbox_inches="tight")
        plt.close()

        if lymandatafile:
            xint_raw = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = open(lymandatafile,"rb")
            first = True
            signal_raw = []
            i=0
            for line in f:
                if not first:
                    data = line.split()
                    if data[1] == "nan":
                        xinit_raw.pop(i)
                    else:
                        signal_raw.append(float(data[1]))
                    i+=1
                first=False
            f.close()

            xint_raw = np.array(xint_raw)

            E_per_photon = 6.626e-34*3.0e8/1216.0e-10
            signal_raw = np.array(signal_raw) * E_per_photon
            xint = range(1,len(signal)+1)
            plt.plot(xint_raw,signal_raw,"-o")
            plt.errorbar(xint,signal,yerr=signal*signal_err,fmt="-o",capsize=3)
            plt.xticks(xint)
            plt.xlabel("Detector number")
            plt.ylabel("Predicted signal (W / m^2-sr)") 
            plt.legend(["Observation","DEGAS2"])
            plt.savefig("signals.pdf",bbox_inches="tight")
            plt.close()
        else:
            xint = range(1,len(signal)+1)
            plt.errorbar(xint,signal,yerr=signal*signal_err,fmt="-o",capsize=3)
            plt.xticks(xint)
            plt.xlabel("Detector number")
            plt.ylabel("Predicted signal (W / m^2-sr)") 
            plt.savefig("signals.pdf",bbox_inches="tight")
            plt.close()

            E_per_photon = 6.626e-34*3.0e8/1216.0e-10
            plt.figure(figsize=(4,3))
            plt.errorbar(xint,signal/E_per_photon,yerr=signal*signal_err,fmt="-o",capsize=3)
#            plt.xticks(xint)
            plt.xlabel(r"Detector number")
            plt.ylabel(r"Predicted brightness (ph.$/s m^2-\mathrm{sr}$)") 
            plt.savefig("photonsignal.png",bbox_inches="tight")
            plt.close()
    return x,z,density,density_err,emission,signal,signal_err, vols, ioniz, esource


# Returns neutral density and detector signals 
def process_output_with_detectors(outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc",plot=False,lymandatafile=None):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)

    geomdata = nc.Dataset(geometryfilename)
    zone_coords_3D = geomdata["zone_center"]
    zone_type = geomdata["zone_type"]
    x_zone = zone_coords_3D[:][0]
    z_zone = zone_coords_3D[:][2]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    emission_idx = -1
    signal_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
        elif "Lyman detector view" in str(nc.chartostring(tallynames[itally])):
            signal_idx = itally
        elif "Lyman emission rate" in str(nc.chartostring(tallynames[itally])):
            emission_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar
        elif "detector " in str(nc.chartostring(varnames[itally])):
            detector_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    signal_base = tallydata["tally_base"][signal_idx]
    emission_base = tallydata["tally_base"][emission_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    Ndetector = tally_indices[signal_idx,0]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    x = []
    z = []
    density = []
    density_err = []
    emission = []
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
           x.append(geomdata["zone_center"][izone,0])
           z.append(geomdata["zone_center"][izone,2])
           density.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,0])
           density_err.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,1])
           emission.append(outputdata["out_post_all"][emission_base+izone,0])

    signal = []
    signal_err = []
    for idet in range(0,Ndetector):
        signal.append(outputdata["out_post_all"][signal_base+idet,0])
        signal_err.append(outputdata["out_post_all"][signal_base+idet,1])

    x = np.array(x)
    z = np.array(z)
    density = np.array(density)
    density_err = np.array(density_err)
    emission = np.array(emission)
    signal = np.array(signal)
    signal_err = np.array(signal_err)

    if plot:
        triang = tri.Triangulation(x,z)
        plt.gca().set_aspect("equal")
        plt.title("Neutral density (m^-3)")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,density)
        plt.colorbar()
        plt.savefig("density.pdf",bbox_inches="tight")
        plt.close()

#        triang = tri.Triangulation(x,z)
#        plt.title("Neutral density (m^-3)")
#        plt.xlabel("x (m)")
#        plt.ylabel("z (m)")
#        plt.tricontourf(triang,np.log(density))
#        plt.colorbar()
#        plt.savefig("logdensity.pdf",bbox_inches="tight")
#        plt.close()


        plt.title("Relative error of neutral density")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,density_err)
        plt.colorbar()
        plt.savefig("error.pdf",bbox_inches="tight")
        plt.close()

        plt.title("Lyman-alpha emission (W / m^2)")
        plt.gca().set_aspect("equal")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,emission)
        plt.gca().set_aspect("equal")
        plt.colorbar()
        plt.savefig("emission.pdf",bbox_inches="tight")
        plt.close()

        if lymandatafile:
            xint_raw = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
            f = open(lymandatafile,"rb")
            first = True
            signal_raw = []
            i=0
            for line in f:
                if not first:
                    data = line.split()
                    if data[1] == "nan":
                        xinit_raw.pop(i)
                    else:
                        signal_raw.append(float(data[1]))
                    i+=1
                first=False
            f.close()

            xint_raw = np.array(xint_raw)

            E_per_photon = 6.626e-34*3.0e8/1216.0e-10
            signal_raw = np.array(signal_raw) * E_per_photon
            xint = range(1,len(signal)+1)
            plt.plot(xint_raw,signal_raw,"-o")
            plt.errorbar(xint,signal,yerr=signal*signal_err,fmt="-o",capsize=3)
            plt.xticks(xint)
            plt.xlabel("Detector number")
            plt.ylabel("Predicted signal (W / m^2-sr)") 
            plt.legend(["Observation","DEGAS2"])
            plt.savefig("signals.pdf",bbox_inches="tight")
            plt.close()
        else:
            xint = range(1,len(signal)+1)
            plt.errorbar(xint,signal,yerr=signal*signal_err,fmt="-o",capsize=3)
            plt.xticks(xint)
            plt.xlabel("Detector number")
            plt.ylabel("Predicted signal (W / m^2-sr)") 
            plt.savefig("signals.pdf",bbox_inches="tight")
            plt.close()

            E_per_photon = 6.626e-34*3.0e8/1216.0e-10
            plt.errorbar(xint,signal/E_per_photon,yerr=signal*signal_err,fmt="-o",capsize=3)
            plt.xticks(xint)
            plt.xlabel("Detector number")
            plt.ylabel("Predicted signal (photons / s m^2-sr)") 
            plt.savefig("photonsignal.pdf",bbox_inches="tight")
            plt.close()
    return x,z,density,density_err,emission,signal,signal_err

def process_output(wallnodes_ordered,outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc",plot=True):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)

    geomdata = nc.Dataset(geometryfilename)
    zone_coords_3D = geomdata["zone_center"]
    zone_type = geomdata["zone_type"]
    x_zone = zone_coords_3D[:][0]
    z_zone = zone_coords_3D[:][2]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    emission_idx = -1
    signal_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
        elif "Lyman emission rate" in str(nc.chartostring(tallynames[itally])):
            emission_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    emission_base = tallydata["tally_base"][emission_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    x = []
    z = []
    density = []
    density_err = []
    emission = []
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
           x.append(geomdata["zone_center"][izone,0])
           z.append(geomdata["zone_center"][izone,2])
           density.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,0])
           density_err.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,1])
           emission.append(outputdata["out_post_all"][emission_base+izone,0])

    x = np.array(x)
    z = np.array(z)
    density = np.array(density)
    density_err = np.array(density_err)

    Nwall = len(wallnodes_ordered)
    xwall = np.array(wallnodes_ordered)
    zwall = np.array(wallnodes_ordered)
    for i in range(0,Nwall):
        xwall[i] = wallnodes_ordered[i].coords[0]
        zwall[i] = wallnodes_ordered[i].coords[1]

    triang = tri.Triangulation(x,z)
    if plot:

        #TODO: Try tripcolor

        plt.title("Neutral density (m^-3)")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.plot(x,z,".")
        plt.plot(xwall,zwall,"-")
        plt.savefig("zones.pdf",bbox_inches="tight")
        plt.close()


        plt.title("Log10(Neutral density (m^-3))")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
#        plt.tricontourf(triang,np.maximum(density,1.0e10),20,locator=ticker.LogLocator(),cmap="Oranges")
        plt.tricontourf(triang,np.log10(np.maximum(density,1.0e10)),20)
        ax = plt.gca()
        ax.set_aspect("equal")
#        plt.tricontourf(triang,density)
        plt.plot(xwall,zwall,"-c")
        plt.colorbar()
        plt.savefig("density.pdf",bbox_inches="tight")
        plt.close()

#        triang = tri.Triangulation(x,z)
#        plt.title("Neutral density (m^-3)")
#        plt.xlabel("x (m)")
#        plt.ylabel("z (m)")
#        plt.tricontourf(triang,np.log(density))
#        plt.plot(xwall,zwall,"-")
#        plt.colorbar()
#        plt.savefig("logdensity.pdf",bbox_inches="tight")
#        plt.close()


        plt.title("Relative error of neutral density")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,density_err)
        plt.colorbar()
        plt.savefig("error.pdf",bbox_inches="tight")
        plt.close()

        plt.title("Lyman-alpha emission (W / m^2)")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,np.maximum(emission,1.0),locator=ticker.LogLocator())
        plt.plot(xwall,zwall,"-")
        plt.colorbar()
        plt.savefig("emission.pdf",bbox_inches="tight")
        plt.close()

    return x,z,density,emission

def get_emission_rates(outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc"):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)
    geomdata = nc.Dataset(geometryfilename)

    zone_volumes = geomdata["zone_volume"]

    x_zone = zone_coords_3D[:][0]
    z_zone = zone_coords_3D[:][2]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    emission_idx = -1
    signal_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
        elif "Lyman emission rate" in str(nc.chartostring(tallynames[itally])):
            emission_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    emission_base = tallydata["tally_base"][emission_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    x = []
    z = []
    density = []
    density_err = []
    emission = []
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
           x.append(geomdata["zone_center"][izone,0])
           z.append(geomdata["zone_center"][izone,2])
           density.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,0])
           density_err.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,1])
           emission.append(outputdata["out_post_all"][emission_base+izone,0])

    x = np.array(x)
    z = np.array(z)
    density = np.array(density)
    density_err = np.array(density_err)

    Nwall = len(wallnodes_ordered)
    xwall = np.array(wallnodes_ordered)
    zwall = np.array(wallnodes_ordered)
    for i in range(0,Nwall):
        xwall[i] = wallnodes_ordered[i].coords[0]
        zwall[i] = wallnodes_ordered[i].coords[1]

    triang = tri.Triangulation(x,z)
    if plot:

        #TODO: Try tripcolor

        plt.title("Neutral density (m^-3)")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.plot(x,z,".")
        plt.plot(xwall,zwall,"-")
        plt.savefig("zones.pdf",bbox_inches="tight")
        plt.close()


        plt.title("Log10(Neutral density (m^-3))")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
#        plt.tricontourf(triang,np.maximum(density,1.0e10),20,locator=ticker.LogLocator(),cmap="Oranges")
        plt.tricontourf(triang,np.log10(np.maximum(density,1.0e10)),20)
        ax = plt.gca()
        ax.set_aspect("equal")
#        plt.tricontourf(triang,density)
        plt.plot(xwall,zwall,"-c")
        plt.colorbar()
        plt.savefig("density.pdf",bbox_inches="tight")
        plt.close()

#        triang = tri.Triangulation(x,z)
#        plt.title("Neutral density (m^-3)")
#        plt.xlabel("x (m)")
#        plt.ylabel("z (m)")
#        plt.tricontourf(triang,np.log(density))
#        plt.plot(xwall,zwall,"-")
#        plt.colorbar()
#        plt.savefig("logdensity.pdf",bbox_inches="tight")
#        plt.close()


        plt.title("Relative error of neutral density")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,density_err)
        plt.colorbar()
        plt.savefig("error.pdf",bbox_inches="tight")
        plt.close()

        plt.title("Lyman-alpha emission (W / m^2)")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,np.maximum(emission,1.0),locator=ticker.LogLocator())
        plt.plot(xwall,zwall,"-")
        plt.colorbar()
        plt.savefig("emission.pdf",bbox_inches="tight")
        plt.close()

    return x,z,density,emission

def get_nup(outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc"):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)
    geomdata = nc.Dataset(geometryfilename)

    zone_volumes = geomdata["zone_volume"]
    zone_type = geomdata["zone_type"]
    zone_coords = geomdata["zone_center"]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    flux_idx = -1
    pres_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
        if "neutral flux vector" in str(nc.chartostring(tallynames[itally])):
            flux_idx = itally
        if "neutral pressure" in str(nc.chartostring(tallynames[itally])):
            pres_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    flux_base = tallydata["tally_base"][flux_idx]
    pres_base = tallydata["tally_base"][pres_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    Nzone_p = 0
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
            Nzone_p += 1
    
    Nsp = tally_indices[dens_idx,1]

    density = np.zeros([Nzone_p,Nsp-1])
    density_err = np.zeros([Nzone_p,Nsp-1])
    flux = np.zeros([3,Nzone_p,Nsp-1])
    pressure = np.zeros([Nzone_p,Nsp-1])
    r = np.zeros([Nzone_p])
    z = np.zeros([Nzone_p])
    vols = np.zeros([Nzone_p,Nsp])
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
            for isp in range(0,Nsp-1):
                density[izone,isp] = outputdata["out_post_all"][dens_base+(1+isp)*Nzone+izone,0]
                density_err[izone,isp] = outputdata["out_post_all"][dens_base+(1+isp)*Nzone+izone,1]
                pressure[izone,isp] = outputdata["out_post_all"][pres_base+(1+isp)*Nzone+izone,0]
                for j in range(0,3):
                    flux[j,izone,isp] = outputdata["out_post_all"][flux_base+3*Nzone*(1+isp)+3*izone+j,0]
                vols[izone] = zone_volumes[izone]
                r[izone] = zone_coords[izone,0]
                z[izone] = zone_coords[izone,2]

    return density,density_err,flux,pressure,vols


def get_density(outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc"):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)
    geomdata = nc.Dataset(geometryfilename)

    zone_volumes = geomdata["zone_volume"]
    zone_type = geomdata["zone_type"]
    zone_coords = geomdata["zone_center"]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    pres_idx = -1
    psource_idx = -1
    msource_idx = -1
    esource_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    density = []
    density_err = []
    r = []
    z = []
    vols = []
    # The following implicitly assumes ions are the second background
    # and the relevant neutrals are the second test species.
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
           density.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,0])
           density_err.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,1])
           vols.append(zone_volumes[izone])
           r.append(zone_coords[izone,0])
           z.append(zone_coords[izone,2])
    density = np.array(density)
    density_err = np.array(density_err)
    r = np.array(r)
    z = np.array(z)
    vols = np.array(vols)

    return r,z,density,density_err, vols


def get_density_and_sources(outputfilename="output.nc",tallyfilename="tally.nc",geometryfilename="geometry.nc"):
    outputdata=nc.Dataset(outputfilename)
    tallydata=nc.Dataset(tallyfilename)
    geomdata = nc.Dataset(geometryfilename)

    zone_volumes = geomdata["zone_volume"]
    zone_type = geomdata["zone_type"]

    # Find the tally indices to use
    tallynames = tallydata["tally_name"]
    Ntally=len(tallynames)
    dens_idx = -1
    pres_idx = -1
    psource_idx = -1
    msource_idx = -1
    esource_idx = -1
    for itally in range(0,Ntally):
        if "neutral density" in str(nc.chartostring(tallynames[itally])):
            dens_idx = itally
        elif "neutral pressure" in str(nc.chartostring(tallynames[itally])):
            pres_idx = itally
        elif "ion source rate" in str(nc.chartostring(tallynames[itally])) and not "by reaction" in str(nc.chartostring(tallynames[itally])):
            psource_idx = itally
        elif "ion momentum source vector" in str(nc.chartostring(tallynames[itally])) and not "by reaction" in str(nc.chartostring(tallynames[itally])):
            msource_idx = itally
        elif "ion energy source" in str(nc.chartostring(tallynames[itally])) and not "by reaction" in str(nc.chartostring(tallynames[itally])):
            esource_idx = itally

    # Find the indices of the independent variables
    zone_idx = -1
    det_idx = -1
    varnames = tallydata["tally_var_list"]
    Nvar = len(varnames)
    for ivar in range(0,Nvar):
        if "zone " in str(nc.chartostring(varnames[itally])):
            zone_idx = ivar

    dens_base = tallydata["tally_base"][dens_idx]
    pres_base = tallydata["tally_base"][pres_idx]
    psource_base = tallydata["tally_base"][psource_idx]
    msource_base = tallydata["tally_base"][msource_idx]
    esource_base = tallydata["tally_base"][esource_idx]

    # tally_tab_index holds the dimensionality of each tally (Ntally x tally_rank_ind)
    tally_indices = tallydata["tally_tab_index"]

    max_tally_rank = len(tally_indices[0,:])

    Nzone = tally_indices[dens_idx,0]

    density = []
    density_err = []
    psource = []
    msource = []
    esource_i = []
    esource_e = []
    vols = []
    # The following implicitly assumes ions are the second background
    # and the relevant neutrals are the second test species.
    for izone in range(0,tally_indices[dens_idx,0]):
        if zone_type[izone] == 2:
           density.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,0])
           density_err.append(outputdata["out_post_all"][dens_base+1*Nzone+izone,1])
           psource.append(outputdata["out_post_all"][psource_base+1*Nzone+izone,0])
           esource_e.append(outputdata["out_post_all"][esource_base+0*Nzone+izone,0])
           esource_i.append(outputdata["out_post_all"][esource_base+1*Nzone+izone,0])
           msource.append(outputdata["out_post_all"][msource_base+3*Nzone+3*izone+1,0])
           vols.append(zone_volumes[izone])
    density = np.array(density)
    density_err = np.array(density_err)
    psource = np.array(psource)
    msource = np.array(msource)
    esource_i = np.array(esource_i)
    esource_e = np.array(esource_e)

    return density,density_err, psource, msource, esource_i, esource_e, vols


def append_tri_file(origfilename,zone_map,ndensity,psource,esource_i,esource_e,vols,newformat=False):
    newfilename = origfilename+".new"
    origfile = open(origfilename,'r')
    newfile = open(newfilename,'w')

    # Read original file line by line. 
    # Repeat line appended by neutral properties in additional columns
    done = False
    izone = 0
    iele = 0
    start = False
    nele = 9999999
    while not done:
        origline = origfile.readline()
        if start:
            if newformat:
#                validflag = int(origline.strip().split(",")[12])
                validflag = int(origline.strip().split(",")[14])
            else:
                validflag = int(origline.strip().split(",")[11])
            iele = int(origline.strip().split(",")[0])
            izone = zone_map[iele]
            if izone == -1:
                newline = origline.strip()+",nan,nan,nan,nan"+"\n"
            else:
                newline = origline.strip()+","+str(ndensity[izone])+","+str(psource[izone]/vols[izone])+","+str(esource_i[izone]/vols[izone])+","+str(esource_e[izone]/vols[izone])+"\n"
        else:
            newline = origline

        newfile.write(newline)
        if origline[0:2] != "//" and not start:
            nele = int(origline.strip().split(",")[0])
            start = True
        if iele >= nele-1:
            done = True

    newfile.close()
    origfile.close()

 

