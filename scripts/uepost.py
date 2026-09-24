import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import h5py
from problem import genStdProblem
from dg2d import DG2D, write_dg2d_header
import subprocess
import os
import shutil
from glob import glob
from tqdm import tqdm
import defineback
from matplotlib.colors import LogNorm
from source import Source, write_db_input
from postprocess import get_output

# Meant to traverse a directory structure set up the following way:
# - KSTAR_geo_parameters: 
#   - Ip300, ...
#     - Brbf2bt_Ip300kA, etc.
# - Ip300,...
#   - ML2d_[..params..]
#     - savedt.hdf5_[....]_successful
#
# All cases share the same mesh. Save geometry.nc in KSTAR_geo_parameters
# Each current has different magnetic field data and indeed has its own background.nc. 
# Globally shared nc files: problem, geometry, tally
# Individual-case nc files: background.nc

class CommonData:

    def __init__(self,Ip,casename,Nflight=10000,Mx=14,My=12,Lx=18):
        h5file = "ML2d_"+casename+"/savedt.hdf5_"+casename+"_successful"
        print(h5file)
        f=h5py.File(h5file,"r")
        Nx = f["com"]["nx"][()]
        Ny = f["com"]["ny"][()]
        rm=f['com']["rm"][()]
        zm=f['com']["zm"][()]
        f.close()

        self.Nx=Nx
        self.Ny=Ny
        self.Mx=Mx
        self.Lx=Lx
        self.My=My

        self.Nflight= Nflight

        rgrid = np.zeros([Nx+1,Ny+1])
        zgrid = np.zeros([Nx+1,Ny+1])

        rgrid[1:Nx+1,1:Ny+1] = rm[1:Nx+1,1:Ny+1,4]
        zgrid[1:Nx+1,1:Ny+1] = zm[1:Nx+1,1:Ny+1,4]
        rgrid[0,1:Ny+1] = rm[1,1:Ny+1,3]
        zgrid[0,1:Ny+1] = zm[1,1:Ny+1,3]
        rgrid[1:Nx+1,0] = rm[1:Nx+1,1,2]
        zgrid[1:Nx+1,0] = zm[1:Nx+1,1,2]
        rgrid[0,0] = rm[1,1,1]
        zgrid[0,0] = zm[1,1,1]

        self.rgrid = rgrid
        self.zgrid = zgrid

        f = open("Brbf2bt_Ip%3dkA"%int(Ip))
        self.Bratio = np.zeros([Nx,Ny])
        for i in range(0,Nx):
            line = f.readline().split()
            for j in range(0,Ny):
                self.Bratio[i,j] = float(line[j+1])
        f.close()
        self.Bratio = np.sqrt(1.0 - self.Bratio**2)

        self.area_inner = np.zeros(Ny)
        self.area_outer = np.zeros(Ny)
        self.area_inner = np.sqrt( (rgrid[0,1:Ny+1]-rgrid[0,0:Ny])**2 + (zgrid[0,1:Ny+1]-zgrid[0,0:Ny])**2)
        self.area_outer = np.sqrt( (rgrid[Nx,1:Ny+1]-rgrid[Nx,0:Ny])**2 + (zgrid[Nx,1:Ny+1]-zgrid[Nx,0:Ny])**2)
        self.area_inner *= 0.5*(rgrid[0,1:Ny+1]+rgrid[0,0:Ny])*2.0*np.pi
        self.area_outer *= 0.5*(rgrid[Nx,1:Ny+1]+rgrid[Nx,0:Ny])*2.0*np.pi

def setup_case(c,casename,Ip):
    Nx = c.Nx
    Ny = c.Ny

    os.chdir("ML2d_"+casename)

    h5file = "savedt.hdf5_"+casename+"_successful"
    f=h5py.File(h5file,"r")
    ni_zone = f["bbb"]["ni"][1:Nx+1,1:Ny+1,:]
    Te_zone = f["bbb"]["te"][1:Nx+1,1:Ny+1]
    Ti_zone = f["bbb"]["ti"][1:Nx+1,1:Ny+1]
    ng_zone = f["bbb"]["ng"][1:Nx+1,1:Ny+1]
    up_zone = f["bbb"]["up"][1:Nx+1,1:Ny+1,:]
    f.close()

    source_inner = ni_zone[0,:,0]*up_zone[0,:,0]*c.area_inner*c.Bratio[0,:]
    source_outer = ni_zone[-1,:,0]*up_zone[-1,:,0]*c.area_outer*c.Bratio[-1,:]

    source_stratum = np.zeros(2*Ny,dtype=int)
    source_segment = np.zeros(2*Ny,dtype=int)
    source_strength = np.zeros(2*Ny)

    source_strength[0:Ny] = -source_inner[::-1]
    source_strength[Ny:2*Ny] = source_outer[:]
    source_stratum[0:Ny] = [Nx*Ny+1]*Ny
    source_stratum[Ny:2*Ny] = [Nx*Ny+2]*Ny
    source_segment[0:Ny] = range(0,Ny)
    source_segment[Ny:2*Ny] = range(0,Ny)

    defineback.generate_sourcefile(source_stratum,source_segment,source_strength,sourcefilename="sourcefile")
    defineback.write_plasmafile(ni_zone[:,:,0].flatten(),Te_zone[:,:].flatten()/1.602e-19,Ti_zone[:,:].flatten()/1.602e-19,uy_zone=up_zone[:,:,0].flatten())
    sgroup = Source(c.Nflight,"plate","D",rootspecies="D+",specify_flux=False)
    write_db_input([sgroup])

    shutil.copy("../../degas2.in",".")
#    shutil.copy("../../tally.nc",".")

    subprocess.run("defineback db.in",shell=True)

    os.chdir("..")

def loop_through(Ip,Nflight):
    os.chdir("Ip%3d"%int(Ip))
    dirs = glob("ML2d_*")
    Ndir = len(dirs)
    common=CommonData(Ip,dirs[0][5:],Nflight=Nflight,Mx=14,My=12,Lx=18)


    for i in tqdm(range(0,Ndir)):
        setup_case(common,dirs[i][5:],Ip)
        
    os.chdir("..")
