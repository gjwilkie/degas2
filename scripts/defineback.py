import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy import interpolate
import numpy as np
import dg2d
import sys
from shapely.geometry import Point as sPoint
from shapely.geometry.polygon import Polygon as sPolygon
from importlib import reload 
import matplotlib.tri as tri
import source

class DB:
    """ 
    Class representing all input for defineback

    Attributes:
        sgroups: list of Source type
    """


    def __init__(self):
        """
        Constructor for defineback object. Generates for minimal defaults.
        """
        self.sgroups = []
        self.ne_zone = []
        self.Te_zone = []
        self.Ti_zone = []
        self.ui_zone = None
    
    def write_files(self):
        write_plasmafile(self.ne_zone,self.Te_zone,self.Ti_zone,self.ui_zone)


# Populates plasma density and temperature by zone. Depending on if the point (zone center) is inside the separatrix,
# this will either interpolate based on psi, or interpolate on a 2D table based on more general R and Z data
def get_zone_plasma_data_through_psi_inside_sep(zone_coords,R_outside,Z_outside,ne_outside,Te_outside,psifunc,R_inside,ne_inside,Te_inside,R_sep,Z_sep,hfs_R_lim=-1,hfs_fac=1.0,lfs_R_lim=-1,hfs_ne=None,hfs_Te=None,plot=False):
    psi_data = []

    for R in R_inside:
        # Get average value on psi on each surface. Nominally all points should have equal psi
        psi_data.append(psifunc(R,0.0))

#    ne_func_inside = interpolate.interp1d(psi_data,ne_inside,fill_value=(ne_inside[0],ne_inside[-1]),bounds_error=False)
#    Te_func_inside = interpolate.interp1d(psi_data,Te_inside,fill_value=(Te_inside[0],Te_inside[-1]),bounds_error=False)
    ne_func_inside = interpolate.interp1d(psi_data,ne_inside,fill_value="extrapolate",bounds_error=False)
    Te_func_inside = interpolate.interp1d(psi_data,Te_inside,fill_value="extrapolate",bounds_error=False)

    ne_zone = []
    Te_zone = []

    sep = []
    for i in range(0,len(R_sep)):
        sep.append( (R_sep[i], Z_sep[i]) )

    sep_poly = sPolygon(sep)

    R_out = []
    Z_out = []
    ne_out = []
    Te_out = []

    for point in zone_coords:
        point_temp = sPoint(point[0],point[1])

        if sep_poly.contains(point_temp):
            psi = psifunc(point[0],point[1])
            ne_zone.append(ne_func_inside(psi))
            Te_zone.append(Te_func_inside(psi))
        elif hfs_R_lim > 0.0 and point[0] < hfs_R_lim:
            ne_zone.append(hfs_ne)
            Te_zone.append(hfs_Te)
        elif lfs_R_lim > 0.0 and point[0] > lfs_R_lim:
            R_lcfs = np.max(R_sep)
            psi = psifunc(R_lcfs,0.0)
            ne_zone.append(hfs_fac*ne_func_inside(psi))
            Te_zone.append(Te_func_inside(psi))
        else:
#            ne = interpolate.griddata(np.vstack((R_outside,Z_outside)).transpose(),ne_outside,(point[0],point[1]),method="nearest")
#            Te = interpolate.griddata(np.vstack((R_outside,Z_outside)).transpose(),Te_outside,(point[0],point[1]),method="nearest")
#            ne_zone.append(ne)
#            Te_zone.append(Te)
             ne = interpolate.griddata((R_outside, Z_outside), ne_outside, (point[0], point[1]), method='nearest')
             Te = interpolate.griddata((R_outside, Z_outside), Te_outside, (point[0], point[1]), method='nearest')
             ne_zone.append(ne)
             Te_zone.append(Te)

        R_out.append(point[0])
        Z_out.append(point[1])
        ne_out.append(ne_zone[-1])
        Te_out.append(Te_zone[-1])

    R_out = np.array(R_out)
    Z_out = np.array(Z_out)
    ne_out = np.array(ne_out)
    Te_out = np.array(Te_out)
    Te_zone = np.array(Te_zone)
    ne_zone = np.array(ne_zone)

    if plot:
        triang = tri.Triangulation(R_out,Z_out)
        plt.title("Electron density")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
    #    plt.tricontourf(triang,ne_out,levels=np.linspace(0,4.0e18,8))
        plt.tricontourf(triang,ne_out)
        plt.colorbar()
        plt.savefig("ne.pdf",bbox_inches="tight")
        plt.close()
    
    
        ax = plt.axes(projection='3d')
        ax.scatter3D(R_out,Z_out,ne_out,c=ne_out,cmap="Blues")
        ax.view_init(azim=255,elev=10)
    #    ax.scatter3D(R_outside,Z_outside,ne_outside,c=ne_outside,cmap="Greens")
    #    ax.set_xlim([0.1,0.2])
    #    ax.set_ylim([-0.3,0.3])
        plt.savefig("ne3d.pdf",bbox_inches="tight")
        plt.cla()
        plt.clf()
        plt.close()
    
        plt.title("Electron temperature")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.tricontourf(triang,Te_out)
        plt.colorbar()
        plt.savefig("Te.pdf",bbox_inches="tight")
        plt.close()
    #
        plt.title("Zone centers")
        plt.xlabel("x (m)")
        plt.ylabel("z (m)")
        plt.plot(R_out,Z_out,"+")
        plt.savefig("zones.pdf",bbox_inches="tight")
        plt.close()


    return ne_zone, Te_zone


def get_zone_plasma_data_through_psi(zone_coords,ne_data,Te_data,psifunc,psi_data=None,R_data=None):

    if R_data != None:
        psi_data = []
        for R in R_data:
            # Get average value on psi on each surface. Nominally all points should have equal psi
            psi_data.append(psifunc(R,0.0))
    else:
        if psi_data == None:
            print("ERROR: Must specify either R_data or psi_data in get_zone_plasma_data_through_psi")

    ne_func = interpolate.interp1d(psi_data,ne_data,fill_value=(ne_data[0],ne_data[-1]))
    Te_func = interpolate.interp1d(psi_data,Te_data,fill_value=(Te_data[0],Te_data[-1]))

    ne_zone = []
    Te_zone = []

    for point in zone_coords:
        # scipy.interpolate.RectBivariateSpline returns an ndim=2 array for scalar input. 
        psi = np.squeeze(psifunc(point[0],point[1]))
        ne_zone.append(ne_func(psi))
        Te_zone.append(Te_func(psi))

    return ne_zone, Te_zone

def write_plasmafile(ne_zone,Te_zone,ni_zone,Ti_zone,ui_zone=None,b=None,plasmafilename="plasmafile.txt",ux_zone=None,uy_zone=None,uz_zone=None):
    """ Function responsible for actually writing the degas2 plasmafile.txt.
    """
    pfile = open(plasmafilename,'w')
    pfile.write("zone      T(1)         N(1)        T(2)        N(2)")
    if (not ui_zone is None) or (not ux_zone is None):
        pfile.write("      V1(2)        V2(2)       V3(2)       ") 
    pfile.write("\n")

    Nzone = len(ne_zone)
    for idx in range(0,Nzone):
        pfile.write(str(idx+1)+"  "+str(Te_zone[idx])+"  "+str(ne_zone[idx])
                +"  "+str(Ti_zone[idx])+"  "+str(ni_zone[idx]))
        if not ux_zone is None:
            ux = ux_zone[idx]
            uy = uy_zone[idx]
            uz = uz_zone[idx]
            pfile.write("  "+str(ux)+"  "+str(uy)+"  "+str(uz))
        elif not ui_zone is None:
            ux = ui_zone[idx]*b[idx,0]
            uy = ui_zone[idx]*b[idx,1]
            uz = ui_zone[idx]*b[idx,2]
            pfile.write("  "+str(ux)+"  "+str(uy)+"  "+str(uz))
        pfile.write("\n")
    pfile.close()

def write_sourcefile(sourcefilename,ne_zone,vpar_zone,area_zone,plasma_sector,sector_strata_segment,sector_zone,strata,exitstratum=-1):
    sfile = open(sourcefilename,'w')

    def write_array(label,data):
      sfile.write("#\n"+label+"\n#\n")

      N = len(data)
      for idx in range(0,N):
        sfile.write(str(data[idx])+"  ")
        if (idx+1)%10 == 0 or (idx == (N-1)):
           sfile.write("\n")

    stratum = []
    segment = []
    dens = []
    vpar = []
    area = []

    for iplasma in plasma_sector[1:]:
        if not strata[iplasma] == exitstratum:
          izone = sector_zone[iplasma]-1
          stratum.append(strata[iplasma])
          segment.append(sector_strata_segment[iplasma])
          dens.append(ne_zone[izone])
          area.append(area_zone[izone])
          vpar.append(vpar_zone[izone])
      
    write_array("stratum",stratum) 
    write_array("segment",segment) 
    write_array("N(2)",dens) 
    write_array("V_PAR",vpar) 
    write_array("AREA",area) 

    sfile.close()

    source=np.array(dens)*np.array(vpar)*np.array(area)

    segment = np.array(segment)
    dens = np.array(dens)
    vpar = np.array(vpar)
    source = np.array(source)
    area = np.array(area)

    idx = np.argsort(segment)
    segment = segment[idx]
    dens = dens[idx]
    vpar = vpar[idx]
    source = source[idx]
    area = area[idx]

    plt.plot(np.array(segment),source,"o")
    plt.savefig("source.pdf",bbox_inches="tight")
    plt.close()

    plt.plot(np.array(segment),np.array(dens)*np.array(vpar),"o")
    plt.savefig("sourceflux.pdf",bbox_inches="tight")
    plt.close()

    plt.plot(np.array(segment),np.array(dens),"-")
    plt.savefig("walldens.pdf",bbox_inches="tight")
    plt.close()
    plt.plot(np.array(segment),np.array(vpar),"-")
    plt.savefig("wallvpar.pdf",bbox_inches="tight")
    plt.close()
    plt.plot(np.array(segment),np.array(area),"-")
    plt.savefig("wallarea.pdf",bbox_inches="tight")
    plt.close()


def generate_plasma_files_from_zone_data(ne_zone,Te_zone,Ti_zone,vpar_zone,geomfilename="geometry.nc",ionmass=1.67e-27,trapped_fraction=0.0,plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt"):

    ncdata = nc.Dataset(geomfilename)
    zone_coords_3D = ncdata["zone_center"]
    zone_type = ncdata["zone_type"]
    plasma_sector = ncdata["plasma_sector"]
    sector_zone = ncdata["sector_zone"]
    strata = ncdata["strata"]
    sector_strata_segment = ncdata["sector_strata_segment"]
    sector_points = ncdata["sector_points"]

    zone_coords = []
    zone_idx = []
    for point in range(0,len(zone_coords_3D)):
        # Store the plasma zones in zone_idx and their center locations in zone_coords
        if zone_type[point] == 2:
            zone_coords.append([zone_coords_3D[point,0],zone_coords_3D[point,2]])
            zone_idx.append(point)

    area_zone = np.zeros(np.size(zone_idx))
    vpar_zone = np.zeros(np.size(zone_idx))
    wallzonecenter_r = np.zeros(np.size(zone_idx))
    wallzonecenter_z = np.zeros(np.size(zone_idx))

    with open(solfile_name,"rb") as f:
        lines = f.readlines()
    f.close()

    file = open(bfieldfilename,"r")
    first = True
    r_data = []
    z_data = []
    Br_data = []
    Bz_data = []
    Bt_data = []
    for line in file:
        if not first:
            first=False
            data = line.split()
            r_data.append(float(data[0]))
            z_data.append(float(data[1]))
            Br_data.append(float(data[2]))
            Bt_data.append(float(data[3]))
            Bz_data.append(float(data[4]))
        first=False
    file.close()

    # Find the zones corresponding to each plasma sector 
    for isector in range(1,len(plasma_sector)):
        psector = plasma_sector[isector]
        izone = sector_zone[psector]
        localidx = zone_idx.index(izone)-1
        vpar_zone[localidx] = np.sqrt(1.602e-19*Te_zone[localidx]/ionmass)

        # The two points that define the sector line segment
        point1 = np.array([sector_points[psector,0,0],sector_points[psector,0,2]])
        point2 = np.array([sector_points[psector,1,0],sector_points[psector,1,2]])

        center = 0.5*(point1+point2)

        # Get the unit vector normal to this surface, a_unit
        diff = point2-point1
        normal = [-diff[1],0.0,diff[0]]
        a_unit = normal/np.linalg.norm(normal)
        fullarea = 2.0*np.pi*center[0]*np.linalg.norm(diff)

        # Get the magnetic field unit vector in the poloidal plane, b_unit
        # Use nearest data point:
        data_idx = np.argmin( np.square(center[0]-r_data) + np.square(center[1]-z_data))

        Br=Br_data[data_idx]
        Bz=Bz_data[data_idx]
        Bt=Bt_data[data_idx]

        b_unit = [Br,Bt,Bz]/np.linalg.norm([Br,Bt,Bz])

        wallzonecenter_r[localidx] = center[0]
        wallzonecenter_z[localidx] = center[1]

        area_zone[localidx] = fullarea*np.abs(np.dot(b_unit,a_unit))

    write_plasmafile(plasmafilename,ne_zone,Te_zone,Ti_zone)

# Generates a plasma file for use in defineback from two data sources: n(psi), with psi(r,z) inside separatrix
# and as a general function n(r,z) for outside separatrix.
# Separatrix given as a polygon of points: arrays r_sep, z_sep
def generate_plasma_file_with_psi_and_rz(solfile_name,psifunc,tsfile_name,R_sep,Z_sep,TiTe_ratio=1.0,geomfilename="geometry.nc",bfieldfilename="gs_fields.dat",ionmass=1.67e-27,hfs_R_lim=-1.0,hfs_fac=1.0,lfs_R_lim=-1.0,trapped_fraction=0.0,hfs_ne=None,hfs_Te=None,plot=False,plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt"):

    ncdata = nc.Dataset(geomfilename)
    zone_coords_3D = ncdata["zone_center"]
    zone_type = ncdata["zone_type"]
    plasma_sector = ncdata["plasma_sector"]
    sector_zone = ncdata["sector_zone"]
    strata = ncdata["strata"]
    sector_strata_segment = ncdata["sector_strata_segment"]
    sector_points = ncdata["sector_points"]

    zone_coords = []
    zone_idx = []
    for point in range(0,len(zone_coords_3D)):
        # Store the plasma zones in zone_idx and their center locations in zone_coords
        if zone_type[point] == 2:
            zone_coords.append([zone_coords_3D[point,0],zone_coords_3D[point,2]])
            zone_idx.append(point)

    area_zone = np.zeros(np.size(zone_idx))
    vpar_zone = np.zeros(np.size(zone_idx))
    wallzonecenter_r = np.zeros(np.size(zone_idx))
    wallzonecenter_z = np.zeros(np.size(zone_idx))

    with open(solfile_name,"rb") as f:
        lines = f.readlines()
    f.close()

    R_outside = []
    Z_outside = []
    ne_outside = []
    Te_outside = []

    for line in lines:
        if float(line.split()[2]) > 0.0 and float(line.split()[3]) > 0.0:
            R_outside.append(float(line.split()[0]))
            Z_outside.append(float(line.split()[1]))
            ne_outside.append(float(line.split()[2]))
            Te_outside.append(float(line.split()[3]))
    R_outside = np.array(R_outside)
    Z_outside = np.array(Z_outside)
    ne_outside = np.array(ne_outside)
    Te_outside = np.array(Te_outside)

    with open(tsfile_name,"rb") as f:
        lines = f.readlines()
    f.close()

    R_inside = []
    ne_inside = []
    Te_inside = []

    first = True
    for line in lines:
        if not first:
            R_inside.append(0.01*float(line.split()[0]))
            ne_inside.append(float(line.split()[2]))
            Te_inside.append(float(line.split()[3]))
        first = False
    R_inside = np.array(R_inside)
    ne_inside = np.array(ne_inside)
    Te_inside = np.array(Te_inside)
        
    ne_zone, Te_zone = get_zone_plasma_data_through_psi_inside_sep(zone_coords,R_outside,Z_outside,ne_outside,Te_outside,psifunc,R_inside,ne_inside,Te_inside,R_sep,Z_sep,hfs_R_lim,hfs_fac,lfs_R_lim,hfs_ne,hfs_Te,plot=plot)

    file = open(bfieldfilename,"r")
    first = True
    r_data = []
    z_data = []
    Br_data = []
    Bz_data = []
    Bt_data = []
    for line in file:
        if not first:
            first=False
            data = line.split()
            r_data.append(float(data[0]))
            z_data.append(float(data[1]))
            Br_data.append(float(data[2]))
            Bt_data.append(float(data[3]))
            Bz_data.append(float(data[4]))
        first=False
    file.close()

    # Find the zones corresponding to each plasma sector 
    for isector in range(1,len(plasma_sector)):
        psector = plasma_sector[isector]
        izone = sector_zone[psector]
        localidx = zone_idx.index(izone)-1
        vpar_zone[localidx] = np.sqrt(1.602e-19*Te_zone[localidx]/ionmass)

        # The two points that define the sector line segment
        point1 = np.array([sector_points[psector,0,0],sector_points[psector,0,2]])
        point2 = np.array([sector_points[psector,1,0],sector_points[psector,1,2]])

        center = 0.5*(point1+point2)

        # Get the unit vector normal to this surface, a_unit
        diff = point2-point1
        normal = [-diff[1],0.0,diff[0]]
        a_unit = normal/np.linalg.norm(normal)
        fullarea = 2.0*np.pi*center[0]*np.linalg.norm(diff)

        # Get the magnetic field unit vector in the poloidal plane, b_unit
        # Use nearest data point:
        data_idx = np.argmin(np.sqrt( np.square(center[0]-r_data) + np.square(center[1]-z_data)) )

        Br=Br_data[data_idx]
        Bz=Bz_data[data_idx]
        Bt=Bt_data[data_idx]

        b_unit = [Br,Bt,Bz]/np.linalg.norm([Br,Bt,Bz])

        wallzonecenter_r[localidx] = center[0]
        wallzonecenter_z[localidx] = center[1]

        area_zone[localidx] = fullarea*np.abs(np.dot(b_unit,a_unit))

    write_plasmafile(ne_zone,Te_zone,TiTe_ratio*Te_zone)

    plt.plot(wallzonecenter_r,wallzonecenter_z,'.')
    plt.savefig("wallzonecenter.pdf")
    plt.clf()

    write_sourcefile(sourcefilename,(1.0-trapped_fraction)*np.array(ne_zone),vpar_zone,area_zone,plasma_sector,sector_strata_segment,sector_zone,strata,2)

def generate_sourcefile(strata,segments,source_strength,sourcefilename="sourcefile.txt"):
   sfile = open(sourcefilename,"w") 
   def write_array(label,data):
       sfile.write("#\n"+label+"\n#\n")

       N = len(data)
       for idx in range(0,N):
         sfile.write(str(data[idx])+"  ")
         if (idx+1)%10 == 0 or (idx == (N-1)):
           sfile.write("\n")

   write_array("stratum",strata) 
   write_array("segment",segments)
   write_array("F",source_strength) 
   sfile.close()


# Generates a plasma file for use in defineback
# Arguments:
#   R_data: array of the R coordinates at which ne and Te are given at Z=0
#   Z_data: array of corresponding the Z coordinates 
#   ne_data, Te_data: the arrays of plasma density and temperature data defined at the (R_data,Z_data). ne_data should be in units of m^-3 and Te_data in units of eV
#   TiTe_ratio: a float value that provides Ti/Te, used to infer Ti from Te uniformly. To eventually replace with a separate array.
#   geomfilename: the name of the geometry .nc file that contains the zone data
#   plasmafilename (optional): the name and/or path of the plasma file to write
def generate_plasma_file(R_data,Z_data,ne_data,Te_data,TiTe_ratio,psifunc,geomfilename,bfieldfilename="gs_fields.dat",ionmass=1.66e-27,plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt"):

    ncdata = nc.Dataset(geomfilename)
    zone_coords_3D = ncdata["zone_center"]
    zone_type = ncdata["zone_type"]
    plasma_sector = ncdata["plasma_sector"]
    sector_zone = ncdata["sector_zone"]
    strata = ncdata["strata"]
    sector_strata_segment = ncdata["sector_strata_segment"]
    sector_points = ncdata["sector_points"]

    zone_coords = []
    zone_idx = []
    for point in range(0,len(zone_coords_3D)):
        # Store the plasma zones in zone_idx and their center locations in zone_coords
        if zone_type[point] == 2:
            zone_coords.append([zone_coords_3D[point,0],zone_coords_3D[point,2]])
            zone_idx.append(point)

    area_zone = np.zeros(np.size(zone_idx))
    vpar_zone = np.zeros(np.size(zone_idx))
   
    ne_zone, Te_zone = get_zone_plasma_data(zone_coords,R_data,Z_data,ne_data,Te_data)

    #diag_r = np.zeros(len(plasma_sector))
    #diag_z = np.zeros(len(plasma_sector))

    file = open(bfieldfilename,"r")
    first = True
    r_data = []
    z_data = []
    Br_data = []
    Bt_data = []
    Bz_data = []
    for line in file:
        if not first:
            first=False
            data = line.split()
            r_data.append(float(data[0]))
            z_data.append(float(data[1]))
            Br_data.append(float(data[2]))
            Bt_data.append(float(data[3]))
            Bz_data.append(float(data[4]))
        first=False
    file.close()

    # Find the zones corresponding to each plasma sector 
    for isector in range(1,len(plasma_sector)):
        psector = plasma_sector[isector]
        izone = sector_zone[psector]-1
        localidx = zone_idx.index(izone)
        vpar_zone[localidx] = np.sqrt(1.602e-19*Te_zone[localidx]/ionmass)

        # The two points that define the sector line segment
        point1 = np.array([sector_points[psector,0,0],sector_points[psector,0,2]])
        point2 = np.array([sector_points[psector,1,0],sector_points[psector,1,2]])

        center = 0.5*(point1+point2)

        # Get the unit vector normal to this surface, a_unit
        diff = point2-point1
        normal = [-diff[1],0.0,diff[0]]
        a_unit = normal/np.linalg.norm(normal)
        fullarea = 2.0*np.pi*center[0]*np.linalg.norm(diff)

        # Get the magnetic field unit vector in the poloidal plane, b_unit
        # Use nearest data point:
        data_idx = np.argmin(np.sqrt( np.square(center[0]-r_data) + np.square(center[1]-z_data)) )

        Br=Br_data[data_idx]
        Bt=Bt_data[data_idx]
        Bz=Bz_data[data_idx]

        b_unit = [Br,Bt,Bz]/np.linalg.norm([Br,Bt,Bz])

        area_zone[localidx] = fullarea*np.abs(np.dot(b_unit,a_unit))

    write_plasmafile(plasmafilename,ne_zone,Te_zone,TiTe_ratio*Te_zone)

    write_sourcefile(sourcefilename,ne_zone,vpar_zone,area_zone,plasma_sector,sector_strata_segment,sector_zone,strata)


# Generates a plasma file for use in defineback
# Arguments:
#   R_data: array of the R coordinates at which ne and Te are given at Z=0
#   ne_data, Te_data: the arrays of plasma density and temperature data defined at the R_data locations. ne_data should be in units of m^-3 and Te_data in units of eV
#   TiTe_ratio: a float value that provides Ti/Te, used to infer Ti from Te uniformly. To eventually replace with a separate array.
#   psifunc: a function passed as an argument. This function should take R,Z as arguments and return psi
#   plasmafilename (optional): the name and/or path of the plasma file to write
def generate_plasma_file_through_psi(ne_data,Te_data,Ti_data,psifunc,geomfilename="geometry.nc",bfieldfilename=None,
                                     ionmass=1.66e-27,plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt",
                                     R_data=None,psi_data=None,rot_data=None,ni_data=None):
    """ Function to generate a plasmafile for use in defineback.
        This is the proimary function to create a background plasma given 1D profiles.
    Args:
        ne_data: (1D array-like)
        Te_data: (1D array-like)
        Ti_data: (1D array-like)
        psifunc: (callable)
    Kwargs:
        geomfilename: (str)
        bfieldfilename: (str or None)
        ionmass: (float)
        plasmafilename: (str)
        sourcefilename: (str)
        R_data: (1D array-like or None)
        psi_data: (1D array-like or None)
        rot_data: (1D array-like or None)
	ni_data: (1D array-like or None) if None, ni = ne.
    """
    if ni_data is None:
        ni_data = ne_data
    # ----------------
    # Load information from the DEGAS2 geometry.nc file
    ncdata = nc.Dataset(geomfilename)
    zone_coords_3D = ncdata["zone_center"] # (N_zone, 3) ~ [x,y,z] of each zone center.
    zone_type = ncdata["zone_type"]
    plasma_sector = ncdata["plasma_sector"]
    sector_zone = ncdata["sector_zone"]
    strata = ncdata["strata"]
    sector_strata_segment = ncdata["sector_strata_segment"]
    sector_points = ncdata["sector_points"]

    # ----------------
    # Unpack 
    zone_coords = []
    zone_idx = []
    for point in range(0,len(zone_coords_3D)):
        # Store the plasma zones in zone_idx and their center locations in zone_coords
        if zone_type[point] == 2:
            zone_coords.append([zone_coords_3D[point,0],zone_coords_3D[point,2]])
            zone_idx.append(point)
    Nzone = np.size(zone_idx)

    area_zone = np.zeros(Nzone)
    vpar_zone = np.zeros(Nzone)
    ux_zone = np.zeros(Nzone)
    uy_zone = np.zeros(Nzone)
    uz_zone = np.zeros(Nzone)
   
    if R_data == None:
        ne_zone, Te_zone = get_zone_plasma_data_through_psi(zone_coords,ne_data,Te_data,psifunc,psi_data=psi_data)
        ni_zone, Ti_zone = get_zone_plasma_data_through_psi(zone_coords,ni_data,Ti_data,psifunc,psi_data=psi_data)
    else:
        ne_zone, Te_zone = get_zone_plasma_data_through_psi(zone_coords,R_data,ne_data,Te_data,psifunc,R_data=R_data)
        ni_zone, Ti_zone = get_zone_plasma_data_through_psi(zone_coords,R_data,ni_data,Ti_data,psifunc,R_data=R_data)

    if rot_data != None:
        rot_func = interpolate.interp1d(psi_data,rot_data,fill_value=(rot_data[0],rot_data[-1]))
        for izone in range(0,Nzone):
            # v_phi = R*omega
            uy_zone[izone] = zone_coords[izone][0]*rot_func(psifunc(zone_coords[izone][0],zone_coords[izone][1]))

    #diag_r = np.zeros(len(plasma_sector))
    #diag_z = np.zeros(len(plasma_sector))

    if bfieldfilename != None:
        file = open(bfieldfilename,"r")
        first = True
        r_data = []
        z_data = []
        Br_data = []
        Bz_data = []
        for line in file:
            if not first:
                first=False
                data = line.split()
                r_data.append(float(data[0]))
                z_data.append(float(data[1]))
                Br_data.append(float(data[2]))
                Bz_data.append(float(data[4]))
            first=False
        file.close()

        # Find the zones corresponding to each plasma sector 
        for isector in range(1,len(plasma_sector)):
            psector = plasma_sector[isector]
            izone = sector_zone[psector]-1
            localidx = zone_idx.index(izone)
            vpar_zone[localidx] = np.sqrt(1.602e-19*Te_zone[localidx]/ionmass)
    
            # The two points that define the sector line segment
            point1 = np.array([sector_points[psector,0,0],sector_points[psector,0,2]])
            point2 = np.array([sector_points[psector,1,0],sector_points[psector,1,2]])
    
            center = 0.5*(point1+point2)
    
            # Get the unit vector normal to this surface, a_unit
            diff = point2-point1
            normal = [-diff[1],diff[0]]
            a_unit = normal/np.linalg.norm(normal)
            fullarea = 2.0*np.pi*center[0]*np.linalg.norm(diff)
    
            # Get the magnetic field unit vector in the poloidal plane, b_unit
            # Use nearest data point:
            data_idx = np.argmin(np.sqrt( np.square(center[0]-r_data) + np.square(center[1]-z_data)) )
    
            Br=Br_data[data_idx]
            Bz=Bz_data[data_idx]
    
            b_unit = [Br,Bz]/np.linalg.norm([Br,Bz])
    
            area_zone[localidx] = fullarea*np.abs(np.dot(b_unit,a_unit))

        write_sourcefile(sourcefilename,ne_zone,vpar_zone,area_zone,plasma_sector,sector_strata_segment,sector_zone,strata)

    write_plasmafile(ne_zone,Te_zone,ni_zone,Ti_zone,ux_zone=ux_zone,uy_zone=uy_zone,uz_zone=uz_zone,plasmafilename=plasmafilename)


def write_cylindrical_polygon_input(rgrid,ne,Te,TiTe_ratio,S0,R_tot,NR,Nflights,walltemp=300.0,source_sp="H2",plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt"):

    # Write the plasmafile
    pfile = open(plasmafilename,"w")
    pfile.write("zone   T(1)   N(1)   T(2)   N(2)\n")
    for i in range(0,NR):
        r = rgrid[i]
        pfile.write("%d %f %e %f %e\n"%(i+1,Te(r),ne(r),Te(r)*TiTe_ratio,ne(r)))
    pfile.close()

    generate_db_input_nosourcefile(S0,Nflights,[2*NR],dbfilename="db.in",walltemp=walltemp,source_sp=source_sp)


def write_cylinder_db_input(rgrid,ne,Te,Ti,S0,Nflights,source_stratum,walltemp=300.0,source_sp="H2",plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt",t0=-1.0,tf=-1.0,init=False):

    NR=len(rgrid)

    # Write the plasmafile
    pfile = open(plasmafilename,"w")
    pfile.write("zone   T(1)   N(1)   T(2)   N(2)\n")
    for i in range(0,NR):
        r = rgrid[i]
        pfile.write("%d %f %e %f %e\n"%(i+1,Te(r),ne(r),Ti(r),ne(r)))
    pfile.close()

    sourcegroup = source.Source(Nflights,"puff",source_sp,pufftemp=walltemp,strength=S0,stratum=source_stratum,segment=0)
    
    source.write_db_input([sourcegroup],t0=t0,tf=tf,init=init)


def generate_db_input_nosourcefile(source_strength,Nflights,strata,dbfilename="db.in",walltemp=300.0,source_sp="H2",plasmafilename="plasmafile.txt"):

    f = open(dbfilename,"w")
    f.write("plasma_file "+plasmafilename+"\n")
    f.write("new_source_group\n")
    f.write("  source_type puff\n")
    f.write("  source_geom surface\n")
    f.write("  source_species "+source_sp+"\n")
    f.write("  source_root_sp "+source_sp+"\n")
    f.write("  source_puff_temp "+str(walltemp)+"\n")
    f.write("  specify_flux\n")
    f.write("  source_nflights "+str(Nflights)+"\n")
    f.write("  source_strength %e \n"%(source_strength))
    f.write("  source_stratum ")
    for stratum in strata:
        f.write("%d "%stratum)
    f.write("\n")
    f.write("  source_segment ")
    for stratum in strata:
        f.write("* ")
    f.write("\n")
    f.write("end_source_group\n \n")
    f.close()

# Generates an input file for defineback
# Arguments:
#   Nflights: array of integers specifying the total number of flights from each strata
#   dbfilename: name of the defineback input file
# TODO:
# - Figure out how to have a strictly recycling source
# - Scale Nflights with the length of the strata

def generate_db_input(Nflights,dbfilename="db.in",sourcesp="H",specify_flux=False,plasmafilename="plasmafile.txt",sourcefilename="sourcefile.txt"):

    f = open(dbfilename,"w")
    f.write("plasma_file "+plasmafilename+"\n")
    f.write("new_source_group\n")
    f.write("  source_type plate\n")
    f.write("  source_geom surface\n")
    f.write("  source_species "+sourcesp+"\n")
    f.write("  source_root_sp "+sourcesp+"+"+"\n")
    if specify_flux:
        f.write("  specify_flux\n")
    else:
        f.write("  specify_current\n")
    f.write("  source_nflights "+str(Nflights)+"\n")
    f.write("  source_file "+sourcefilename+" row\n")
    f.write("end_source_group\n \n")

#    for i in range(0,len(wallstrata)):
#        f.write("new_source_group\n")
#        f.write("  source_type plate\n")
#        f.write("  source_geom surface\n")
#        f.write("  source_species H2\n")
#        f.write("  source_root_sp H+\n")
#        f.write("  specify_current\n")
#        f.write("  source_nflights "+str(Nflights[i])+"\n")
#        f.write("  source_stratum "+str(wallstrata[i])+"\n")
#        f.write("  source_segment *\n")
#        f.write("  source_strength "+str(source_strength[i])+"\n")
#        f.write("end_source_group\n \n")

def generatePlasmaFileFromFunctions(psi_func,ne_func,Te_func,Ti_func,gfile_name="geometry.nc",pfile_name="plasmafile.txt"):
    ncdata = nc.Dataset(gfile_name)
    zone_coords_3D = ncdata["zone_center"]
    zone_type = ncdata["zone_type"]

    r_zone = []
    z_zone = []
    ne_zone = []
    Te_zone = []
    Ti_zone = []
    zone_idx = []

    for i in range(0,len(zone_coords_3D)):
        if zone_type[i] == 2:
            r_zone.append(zone_coords_3D[i,0])
            z_zone.append(zone_coords_3D[i,2])
            zone_idx.append(i)
            psi = psi_func(r_zone[i],z_zone[i])
            ne_zone.append(ne_func(psi))
            Te_zone.append(Te_func(psi))
            Ti_zone.append(Ti_func(psi))
    r_zone = np.array(r_zone)
    z_zone = np.array(z_zone)
    ne_zone = np.array(ne_zone)
    Te_zone = np.array(Te_zone)
    Ti_zone = np.array(Ti_zone)
    zone_idx = np.array(zone_idx,dtype=int)

    write_plasmafile(ne_zone,Te_zone,Ti_zone,plasmafilename=pfile_name)

def generateSourceFileFromFunction(sfunc,wallnodes,R0,filename="sourcefile.txt",strata=None,integrated_source=-1.0):
    s_eps = 1.0

    sfile = open(filename,"w") 
    def write_array(label,data):
        sfile.write("#\n"+label+"\n#\n")
 
        N = len(data)
        for idx in range(0,N):
            sfile.write(str(data[idx])+"  ")
            if (idx+1)%10 == 0 or (idx == (N-1)):
                sfile.write("\n")

    Nwall = len(wallnodes)
    segments = []
    source_strength = []
    sumsource = 0.0
    for i in range(0,Nwall):

        l = np.sqrt((wallnodes[(i+1)%Nwall].coords[0]-wallnodes[i].coords[0])**2 + (wallnodes[(i+1)%Nwall].coords[1]-wallnodes[i].coords[1])**2)

        rmid = 0.5*(wallnodes[(i+1)%Nwall].coords[0] + wallnodes[i].coords[0])
        zmid = 0.5*(wallnodes[(i+1)%Nwall].coords[1] + wallnodes[i].coords[1])
        theta = np.arctan2(zmid,rmid-R0)

        if sfunc(theta) > s_eps:
            # Segments are actually in reverse order from the wallfile. So some convoluted index algebra is needed
            iseg = np.abs(i-Nwall-1)%Nwall
            segments.append(iseg)
            tempsource = sfunc(theta)
            tempsource = max(tempsource,sfunc(theta+2.0*np.pi))
            tempsource = max(tempsource,sfunc(theta-2.0*np.pi))
            source_strength.append(tempsource)
            sumsource+=tempsource*l*2.0*np.pi*rmid

    source_strength = np.array(source_strength)

    if integrated_source > 0.0:
        source_strength = source_strength*integrated_source/sumsource

    if not strata:
        strata = [2]*len(segments)
    strata = np.array(strata,dtype=int)
    segments = np.array(segments,dtype=int)
    source_strength = np.array(source_strength)

    write_array("stratum",strata) 
    write_array("segment",segments)
    write_array("F",source_strength) 
    sfile.close()



     

    
