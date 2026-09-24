import adios2
import dg2d
import defineback
from polygon import *
import numpy as np
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
import scipy.special as sp
import netCDF4 as nc
import sys
import os
import f90nml
import problem
import source
import subprocess
import postprocess

def setup_xgc_case(meshbase=None,d2_dir=None,wall_material="C",wall_temperature=300.0,runtest=False,skipgeo=False,aux_thickness=0.002):
    """
    Given XGC inputs, constructs DEGAS2 datafiles in working directory.
    Requires the executables: problemsetup, definegeometry2d, defineback, tallysetup
    Args:
        meshbase (partially required): string for the basename of the geometry files, shared by .ele and .node triangulation files. Include full relative path (e.g. "input_dir"). If not provided, will look for xgc.mesh.bp locally. 
        d2_dir (partially required): string for the absolute path of the degas2 installation. Assumes the build directory is populated with appropriate executables equipped with the required synthetic diagnostics. If this is not given, the user is assumed to have copied degas2.in and tally.in from the degas2/data/templates directory to the local working directory AND have degas2/scripts in their PYTHONPATH.
        wall_material (optional): string for the wall material, as specified in degas2 problem inputs. Defaults to "C" for graphite.
        wall_temperature (optional): float for the wall temperature in Kelvin; determines the energy of desorbed products. Defaults to 300.
    """

    adios2meshfile = "xgc.mesh.bp"
    if os.path.isdir(adios2meshfile):
        rz, conn,wallnodes = get_bp_mesh(filename=adios2meshfile)
    elif meshbase == None:
        raise Exception("If xgc.mesh.bp is not present, meshbase must be specified to find the triangulation files in input_dir.")
    else:
        rz, conn, wallnodes = get_tri_mesh(meshbase)

    Nwall = len(wallnodes)
    Ntri = len(conn[:,0])
    Nnode = len(rz[:,0])

    nml = f90nml.read('input')
    ionmass = nml["ptl_param"]["ptl_mass_au"][1]
    Rcoeff = nml["neu_param"]["neu_recycle_rate"]
    try:
        ebin_min = nml["neu_param"]["neu_ebin_min"]
    except:
        ebin_min = 0.1
    try:
        ebin_max = nml["neu_param"]["neu_ebin_max"]
    except:
        ebin_max = 100.0
    try: 
        ebin_num = nml["neu_param"]["neu_ebin_num"]
    except:
        raise Exception("neu_ebin_num is a required input to override default of 1.")
    try: 
        ebin_log = nml["neu_param"]["neu_ebin_log"]
    except:
        ebin_log = True

    spec = gen_problem_for_xgc(wall_material,ionmass)

    if not skipgeo:
        print("Generating geometry input...")
        g = dg2d.DG2D()
        g.define_mesh(rz,conn,progress=True)
        g.set_wallprops(walltemp=wall_temperature,Rcoeff=Rcoeff,material=wall_material)
        g.write_polygonfile = False
        g.write_files(aux_thickness=aux_thickness)
        print("Running definegeometry2d...")
        subprocess.run(d2_dir+"definegeometry2d dg2d.in",shell=True)

    write_dummy_bg_files_aux(Ntri,Nwall,Ntri+1)
    if runtest:
        sgroup = source.Source(10000,"plate",spec,spec+"+",specify_flux=False,sourcefile="sourcefile.txt")
    else:
        sgroup = source.Source(10000,"plt_e_bins",spec,spec+"+",specify_flux=False,sourcefile="sourcefile.txt",e_bin_num=ebin_num,e_bin_max=ebin_max,e_bin_min=ebin_min,e_bin_log=ebin_log)

    source.write_db_input([sgroup])
    print("Running defineback...")
    subprocess.run(d2_dir+"defineback db.in",shell=True)

    subprocess.run(d2_dir+"tallysetup")

    print("Done.")

def gen_problem_for_xgc(wall_material,ionmass):

    eps = 1.0e-3
    if np.abs(ionmass-1.0) < eps:
        spec = "H"
    elif np.abs(ionmass-2.0) < eps:
        spec = "D"
    else:
        raise Exception("ionmass = %f. Currently can only handle H or D main ions."%ionmass)

    testSps = ["0",spec,spec+"2",spec+"2+"]
    backSps = ["e",spec+"+"]
    reactions = ["hionize5",spec.lower()+"chex_const","h2dis","h2ion","h2dision","h2pdision","h2pdis","h2pdisrec"]
    if wall_material == "mirror":
        testSps = ["0",spec]
        reactions = ["hionize5",spec.lower()+"chex_const"]
        pmis = ["hmirror"]
    elif wall_material == "C":
        pmis = ["hdesorbc","h2desorbc",spec.lower()+"reflc"]
    elif wall_material == "mo":
        pmis = ["hdesorbmo","h2desorbmo",spec.lower()+"reflmo"]
    elif wall_material == "Li":
        if spec == "D":
            pmis = ["hdesorbLi","D_refl_vftrim_Li"]
        elif spec == "H":
            pmis = ["hdesorbLi","H_refl_svftrim_Li"]
    else:
        raise Exception("wall_material = %s. Currently can only handle mirror, C, or mo"%wall_material)

    problem.generateProblemInput(testSps,backSps,reactions,[wall_material],pmis)
    subprocess.run("problemsetup",shell=True)

    return spec
    
def get_tri_mesh(basename):
    f = open(basename+".node","r")
    line = f.readline().split()
    Nnode = int(line[0])
    rz = np.zeros([Nnode,2])
    nodeidx = np.zeros([Nnode],dtype=int)
    wallflag = np.zeros([Nnode],dtype=int)
    for i in range(0,Nnode):
        line = f.readline().split()
        nodeidx[i] = int(line[0])
        rz[i,0] = float(line[1])
        rz[i,1] = float(line[2])
        wallflag[i] = int(line[3])
    f.close()

    wallnodes = np.argwhere(wallflag > 0)
    shift = np.min(nodeidx)

    f = open(basename+".ele","r")
    line = f.readline().split()
    Ntri = int(line[0])
    conn = np.zeros([Ntri,3],dtype=int)
    for i in range(0,Ntri):
        line = f.readline().split()
        conn[i,0] = int(line[1])
        conn[i,1] = int(line[2])
        conn[i,2] = int(line[3])
    f.close()

    conn = conn - shift

    return rz, conn, wallnodes



def get_bp_mesh(filename="xgc.mesh.bp",oldfile=False):
    meshfile = adios2.Stream(filename,"rra")
    if oldfile:
        coords = meshfile.read("/coordinates/values")
        connections = meshfile.read("/cell_set[0]/node_connect_list")
    else:
        coords = meshfile.read("rz")
        connections = meshfile.read("nd_connect_list")
    wallnodes = meshfile.read("grid_wall_nodes")
    # wall_nodes refers to list of nodes by starting count from 1,
    # which is confusingly inconsistent with node_connect_list.
    # Make them consistent here.
    for i in range(0,len(wallnodes)):
        wallnodes[i] -= 1
    meshfile.close()
    return coords, connections, wallnodes

def isclockwise(tri):
    if not len(tri.vertices) == 3:
        sys.exit("Not a triangle. isclockwise logic not valid.")
    s = 0.0
    s += (tri.vertices[1].coords[0] - tri.vertices[0].coords[0])*(tri.vertices[1].coords[1] + tri.vertices[0].coords[1])
    s += (tri.vertices[2].coords[0] - tri.vertices[1].coords[0])*(tri.vertices[2].coords[1] + tri.vertices[1].coords[1])
    s += (tri.vertices[0].coords[0] - tri.vertices[2].coords[0])*(tri.vertices[0].coords[1] + tri.vertices[2].coords[1])
    if s >= 0.0:
        return True
    else:
        return False

def write_geometry_files(material="C",recyc=0.99,walltemp=300,use_xgc_mesh=True,polygonfilename="none",trifile_base=None):

    if trifile_base==None:
        coords,connections,wallnode_ids = get_bp_mesh()
        triang = mtri.Triangulation(coords[:,0],coords[:,1],connections)
    else:
        nodefile=trifile_base+".node"
        f=open(nodefile,"r")
        Nnode = int(f.readline().split()[0])
        coords = np.zeros([Nnode,2])
        wallnode_ids = []
        rwall = []
        zwall = []
        for i in range(0,Nnode):
            line = f.readline().split()
            coords[i,0] = float(line[1])
            coords[i,1] = float(line[2])
            if int(line[3]) == 1:
                wallnode_ids.append(i)
                rwall.append(coords[i,0])
                zwall.append(coords[i,1])
        rwall.append(rwall[0])
        zwall.append(zwall[0])
        rwall = np.array(rwall)
        zwall = np.array(zwall)
        f.close()
        elefile=trifile_base+".ele"
        f=open(elefile,"r")
        Ntri = int(f.readline().split()[0])
        connections = np.zeros([Ntri,3],dtype=int)
        for i in range(0,Ntri):
            line = f.readline().split()
            connections[i,0] = int(line[1])-1
            connections[i,1] = int(line[2])-1
            connections[i,2] = int(line[3])-1
        f.close()       
        triang = mtri.Triangulation(coords[:,0],coords[:,1],connections)

        #plt.triplot(triang)
        #plt.plot(rwall,zwall)
        #plt.savefig("mesh.png")
        #plt.close()

    Nnode = len(coords[:,0])
    Ntri = len(connections[:,0])
    Nwall = len(wallnode_ids)

    # Record nodes in local representation
    # Nodes start counting from 0. Triangles from 1.
    vertices = []
    wallnodes = []
    for inode in range(0,Nnode):
        vertices.append(Vertex(inode,coords[inode,0],coords[inode,1]))
        if np.isin(inode,wallnode_ids):
            vertices[-1].wall = True
            wallnodes.append(vertices[-1])
    if not len(wallnodes) == Nwall:
        sys.exit("Found %d wall nodes, but XGC had %d"%(len(wallnodes),Nwall))

    ####################################################
    # Write wallfile
    wallfile = open("wallfile.txt","w")
    wallfile.write("# Wallfile automatically generated by script\n")
    wallfile.write("1\n")
    wallfile.write("#\n")
    wallfile.write("%d\n"%(Nnode))
    wallfile.write("#\n")
    for i in range(0,Nnode):
        x = coords[i,0]
        z = coords[i,1]
        wallfile.write("%f %f\n"%(x,z))
    wallfile.close()

    Rmin = np.min(coords[:,0])
    Rmax = np.max(coords[:,0])
    Zmin = np.min(coords[:,1])
    Zmax = np.max(coords[:,1])

    Zmin_tot = Zmin - (Zmax-Zmin)
    Zmax_tot = Zmax + (Zmax-Zmin)
    Rmax_tot = 1.5*Rmax
    Rmin_tot = max(0.0001,0.5*Rmin)

    dg2dfile = open("dg2d.in","w")
    dg2d.write_dg2d_header(dg2dfile,"cylindrical",Rmin_tot,Rmax_tot,Zmin_tot,Zmax_tot,wallfile_name="wallfile.txt")

    polys = []
    wall_triangles = []
    for itri in range(0,Ntri):
        polys.append(Polygon(id=itri+1))
        for j in range(0,3):
            idx = connections[itri,j]
            polys[-1].add_vertex(vertices[idx])
        if not isclockwise(polys[-1]):
            polys.pop()
            polys.append(Polygon(increment=False,id=itri+1))
            for j in range(0,3):
                idx = connections[itri,2-j]
                polys[-1].add_vertex(vertices[idx])
        nwallnodes = polys[-1].reorder_wallnodes_first()
        if use_xgc_mesh:
            polys[-1].write_plasma_polygon_dg2d(dg2dfile,stratum=itri+1,wallid=1,commonzone=True)

#        print(polys[-1].id,polys[-1].vertices[0].id,polys[-1].vertices[1].id,polys[-1].vertices[2].id)
        if nwallnodes >= 2:
            wall_triangles.append(polys[-1])
#            print(wall_triangles[-1].id,wall_triangles[-1].vertices[0].id,wall_triangles[-1].vertices[1].id,wall_triangles[-1].vertices[2].id)


    # Now all the plasma elements are defined as their own zones (Polygon type).
    # Each polygon's id is the zone number
    # Nodes are oriented clockwise
    # The first wall node is placed first so that the first segment is 
    # along the wall if the polygon is adjacent to the wall

    # Find the appropriate node to start from for breaking up outer 2 polygons
    # (closest to origin)
    first = True
    for wallnode in wallnodes:
        dist = np.linalg.norm(wallnode.coords)
        if first:
            first = False
            mindist = dist
            firstwallnode_id = wallnode.id
        elif dist < mindist:
            mindist = dist
            firstwallnode_id = wallnode.id

    wallnodes_ordered=[]
    wallnodes_ordered.append(vertices[firstwallnode_id])
    wallpoly = Polygon(increment=False,id=-1)
    wallpoly.add_vertex(wallnodes_ordered[-1])
    idx_ordered = []
    idx_ordered.append(list(wallnode_ids).index(wallnodes_ordered[-1].id))
    for i in range(0,len(wallnodes)-1):
#        print(i,wallnodes_ordered[i].id)
        if i==0:
            prevnode = vertices[firstwallnode_id]
        else:
            prevnode = wallnodes_ordered[i-1]
        next_node = find_next_wall_node(wallnodes_ordered[i],prevnode,wallnodes,polys,i==0)
        idx_ordered.append(list(wallnode_ids).index(next_node.id))
        wallnodes_ordered.append(next_node)
        wallpoly.add_vertex(wallnodes_ordered[-1])

    wallpoly_clockwise=Polygon(increment=False,id=-2)
    wallpoly_clockwise.add_vertex(wallpoly.vertices[0])
    for i in range(0,len(wallnodes)-1):
        node = wallpoly.vertices[len(wallnodes)-1-i]
        wallpoly_clockwise.add_vertex(node)

    if not use_xgc_mesh:
        wallpoly_clockwise.write_plasma_polygon_dg2d(dg2dfile,stratum=1,wallid=1,commonzone=False)
        Polygon.close_in_universal_cell(dg2dfile,wallnodes_ordered,0,2,material,recyc,debug=False,clockwise=False,walltemp=walltemp)
    else:
        Polygon.close_in_universal_cell(dg2dfile,wallnodes_ordered,0,Ntri+1,material,recyc,debug=False,clockwise=False,walltemp=walltemp)

    dg2dfile.write("polygon_nc_file "+polygonfilename+"\n")
    dg2dfile.write("end")
    dg2dfile.close()

    return wallnodes_ordered, wall_triangles, idx_ordered

def write_background_files(dt,tstep,tstep_neut,wallnodes_ordered,wall_triangles,ionmass=3.34e-27,xgc1=False):
    coords,connections,wallnode_ids = get_bp_mesh()
    Ntri_xgc = len(connections[:,0])
    Nnode_xgc = len(coords[:,0])
    Nwall = len(wallnode_ids)

    if not len(wallnodes_ordered) == Nwall:
        sys.exit("Size of wallnodes_ordered passed to write_background_files not consistent with XGC data.")

    # Get plasma node data
    if xgc1:
        f=adios2.open("xgc.f3d.%05d.bp"%tstep,"r")
        ne=np.average(f.read("e_den"),axis=1)
        Te_para=np.average(f.read("e_T_para"),axis=1)
        Te_perp=np.average(f.read("e_T_perp"),axis=1)
        Ti_para=np.average(f.read("i_T_para"),axis=1)
        Ti_perp=np.average(f.read("i_T_perp",axis=1))
        ui_para=np.average(f.read("i_u_para"),axis=1)
    else:
        f=adios2.Stream("xgc.f2d.%05d.bp"%tstep,"rra")
        ne=f.read("e_den")
        Te_para=f.read("e_T_para")
        Te_perp=f.read("e_T_perp")
        Ti_para=f.read("i_T_para")
        Ti_perp=f.read("i_T_perp")
        ui_para=f.read("i_u_para")
    f.close()
    f=adios2.Stream("xgc.bfield.bp","rra")
    Bfield = f.read("bfield")
    f.close()

    f=adios2.Stream("xgc.neutrals.%05d.bp"%tstep_neut,"rra")
    raw_source = f.read("wall_source")

    Te = (2.0*Te_perp + Te_para)/3.0
    Ti = (2.0*Ti_perp + Ti_para)/3.0

    # Get zones from geometry file

    # Convert to triangle-averages
    ne_zone = np.zeros(Ntri)
    Te_zone = np.zeros(Ntri)
    Ti_zone = np.zeros(Ntri)
    ui_zone = np.zeros(Ntri)
    b_zone = np.zeros([Ntri,3])
    for itri in range(0,Ntri):
        nodes = connections[itri,:]
        for j in range(0,3):
            ne_zone[itri] += ne[nodes[j]]/3.0
            Te_zone[itri] += Te[nodes[j]]/3.0
            Ti_zone[itri] += Ti[nodes[j]]/3.0
            ui_zone[itri] += np.abs(ui_para[nodes[j]])/3.0
            Bnorm = np.linalg.norm(Bfield[nodes[j],:])
            b_zone[itri,0] += Bfield[nodes[j],0]/(3.0*Bnorm)
            b_zone[itri,1] += Bfield[nodes[j],2]/(3.0*Bnorm)
            b_zone[itri,1] += -Bfield[nodes[j],1]/(3.0*Bnorm)

    # Source is specified by stratum (Ntri+1) and segment
    # Use the outer polygons and the single stratum for them
    # First segment is the small one on the inboard side.
    # Segments go around COUNTER-clockwise, starting with 
    # index 0 of wallnodes_ordered.
    # Challenge is to find which triangle has this segment 
    # (and there should always be only one).
    # This subroutine gets passed a list of triangles which
    # have at least two wall nodes so that this search is 
    # minimized.
    # Also, DEGAS2 has its own sheath model, so we don't want to use
    # the source directly from XGC, else it would get doubly-accelerated.
    seg_to_tri=np.zeros(Nwall,dtype=np.int8)
    area_norm=np.zeros(Nwall)
    area_tot=np.zeros(Nwall)
    for i in range(0,Nwall):
        node1 = wallnodes_ordered[i]
        node2 = wallnodes_ordered[np.mod(i+1,Nwall)]

        diff = np.array(node2.coords) - np.array(node1.coords)
        seglen = np.linalg.norm(diff)
        centerpt = 0.5*(np.array(node2.coords) + np.array(node1.coords))
        area_tot[i] = seglen*2.0*np.pi*centerpt[0]

        # Unit vector normal to segment
        a_unit = [-diff[1],0.0,diff[0]]
        a_unit = a_unit/np.linalg.norm(a_unit)

        found = False
        for walltri in wall_triangles:
            hasnode1 = (walltri.vertices[0].id == node1.id)
            hasnode1 = hasnode1 or (walltri.vertices[1].id == node1.id)
            hasnode1 = hasnode1 or (walltri.vertices[2].id == node1.id)
            hasnode2 = (walltri.vertices[0].id == node2.id)
            hasnode2 = hasnode2 or (walltri.vertices[1].id == node2.id)
            hasnode2 = hasnode2 or (walltri.vertices[2].id == node2.id)
            if hasnode1 and hasnode2:
                if not found:
                    found = True
                    tri_id = walltri.id
                else:
                    sys.exit("Found more than one compatible triangle for wall segment %d"%i)
        if not found:
            sys.exit("Could not find a suitable triangle adjacent to wall segment %d"%i)
        seg_to_tri[i] = tri_id

        area_norm[i] = area_tot[i]*np.abs(np.dot(a_unit,b_zone[tri_id,:]))

    defineback.write_plasmafile(ne_zone,Te_zone,Ti_zone,ui_zone=ui_zone,b=b_zone,plasmafilename="plasmafile.txt")

    # Enforce Bohm criterion?
    #cs = np.sqrt(1.602e-19*Te_zone/ionmass)
    #ui_zone = np.minimum(ui_zone,cs)

    # If using specify_flux, use this
    #source_strength = wall_source / (dt*area_tot)
    # Using specify_current
    #source_strength = raw_source / dt

    # Source strength calculation. Use directly from adjacent plasma properties.
    # Use specify_current for this
    source_strength = np.zeros(Nwall)
    for iwall in range(0,Nwall):
        itri = seg_to_tri[iwall]
        source_strength[iwall] = area_norm[iwall]*ne_zone[itri]*ui_zone[itri]

    strata = [Nwall+1]*Ntri
    segments = list(range(0,Nwall))
    defineback.generate_sourcefile(strata,segments,source_strength)

# Writes plasma background and neutral source files for DEGAS2 from XGC's output.
# Uses a mesh generated by DEGAS2 and interpolates from XGC mesh.
# (another routine uses XGC's mesh directly)
# Input parameters:
# - dt: the size of a neutral timestep in seconds. sml_dt(s) from units.txt 
#   times neu_background_period. Can also contain a factor which divides the 
#   neutral source, rescaling it by a constant factor of (1/dt).
# - tstep: integer index of the XGC timestep
# - tstep_neut: integer index of the neutral step (index of xgc.neutrals can 
#   differ from xgc.f2d
# - wallnodes_ordered: array of indices for the wall nodes as they go around the 
#   limiter counterclockwise.
# - wallnode_order: When using XGC's source, this provides the mappgng that makes
#   its way around the limiter. Calculated by write_geometry_files.
# - use_xgc_source: whether the neutral source is used directly from 
#   neu_weight_wall_lost (now output to xgc.neutrals.*.bp as "wall_source"
# - ion_mass: used in previous iterations where the Bohm criterion was used.
def write_background_files_for_own_mesh(dt,tstep,tstep_neut,wallnodes_ordered,wallnode_order,use_xgc_source=True,ionmass=3.34e-27,xgc1=False):
    coords,connections,wallnode_ids = get_bp_mesh()
    Ntri = len(connections[:,0])
    Nnode = len(coords[:,0])
    Nwall = len(wallnode_ids)

    if not len(wallnodes_ordered) == Nwall:
        sys.exit("Size of wallnodes_ordered passed to write_background_files not consistent with XGC data.")

    # Get plasma node data
    if xgc1:
        f=adios2.open("xgc.f3d.%05d.bp"%tstep,"r")
        ne=np.average(f.read("e_den"),axis=1)
        Te_para=np.average(f.read("e_T_para"),axis=1)
        Te_perp=np.average(f.read("e_T_perp"),axis=1)
        Ti_para=np.average(f.read("i_T_para"),axis=1)
        Ti_perp=np.average(f.read("i_T_perp"),axis=1)
        ui_para=np.average(f.read("i_u_para"),axis=1)
        f.close()
        f=adios2.open("xgc.bfield.bp","r")
        Bfield = f.read("/node_data[0]/values")
        f.close()
    else:
        f=adios2.Stream("xgc.f2d.%05d.bp"%tstep,"rra")
        ne=f.read("e_den")
        Te_para=f.read("e_T_para")
        Te_perp=f.read("e_T_perp")
        Ti_para=f.read("i_T_para")
        Ti_perp=f.read("i_T_perp")
        ui_para=f.read("i_u_para")
        f.close()
        f=adios2.Stream("xgc.bfield.bp","rra")
        Bfield = f.read("bfield")
        f.close()

#    f=adios2.open("xgc.neutrals.%05d.bp"%tstep_neut,"r")
#    raw_source = f.read("wall_source")
#    f.close()

    # Altnernative raw source:
    f=adios2.Stream("xgc.oneddiag.bp","rra")
    nstep =int(f.available_variables()["step"]['AvailableStepsCount'])
    tstep_map = f.read("step",start=[],count=[],step_selection=[0,nstep])
    f.close()

    f=adios2.Stream("xgc.sheathdiag.bp","r")
    i = 0
    found = False
    for step in f:
        if tstep_map[i] == tstep:
            print("Using sheathdiag source...")
            raw_source = step.read("sheath_ilost")[0,:]/1.602e-19
            print(np.shape(raw_source))
            found = True
        i+=1
    f.close()

    if not found:
        print("Could not find the right timestep")


    Te = (2.0*Te_perp + Te_para)/3.0
    Ti = (2.0*Ti_perp + Ti_para)/3.0

#    Te = Te-np.amin(Te)+1.0
#    Ti = Ti-np.amin(Ti)+10.0
#    ne = ne-np.amin(ne)+1.0e18

    # Get interpolant function in order to 
    r_xgc = coords[:,0]
    z_xgc = coords[:,1]
    ne_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),ne)
    Te_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),Te)
    Ti_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),Ti)
    ui_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),ui_para)
    # XGC's coordinate system is (R,Z,phi_xgc)
    # DEGAS2 is (x,y,z) = (R,-phi,Z) = 
    Bx_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),Bfield[0,:])
    Bz_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),Bfield[1,:])
    By_interp = LinearNDInterpolator(list(zip(r_xgc,z_xgc)),-Bfield[2,:])

    # Get coordinates of zone centers from generated geometry.nc file
    geomfilename="geometry.nc"
    ncdata = nc.Dataset(geomfilename)
    zone_coords_3D = ncdata["zone_center"]
    zone_type = ncdata["zone_type"]
    plasma_sector = ncdata["plasma_sector"]
    sector_zone = ncdata["sector_zone"]
    strata = ncdata["strata"]
    sector_strata_segment = ncdata["sector_strata_segment"]
    sector_points = ncdata["sector_points"]

    # Store the plasma zones in zone_idx and their center locations in zone_coords
    zone_idx = []
    zone_coords = []
    for point in range(0,len(zone_coords_3D)):
        # Only count zones in the plasma volume
        if zone_type[point] == 2:
            zone_idx.append(point)
    Nzone = len(zone_idx)

    zone_coords = np.zeros([Nzone,2])
    for i in range(0,Nzone):
        idx = zone_idx[i]
        zone_coords[i,0] = zone_coords_3D[idx,0]
        zone_coords[i,1] = zone_coords_3D[idx,2]

    r_d2 = np.array(zone_coords[:,0])
    z_d2 = np.array(zone_coords[:,1])
    Ntri = len(r_d2)

    # r_d2 and z_d2 now contain the centers of DEGAS2's plasma zones inthe R-Z plane
    # Interpolate from XGC mesh to get DEGAS2's description of the plasma as
    #   volume averages by zone
    Bfield_zone = np.zeros([Ntri,3])
    b_zone = np.zeros([Ntri,3])
    ne_zone = ne_interp(r_d2,z_d2)
    Te_zone = Te_interp(r_d2,z_d2)
    Ti_zone = Ti_interp(r_d2,z_d2)
    ui_zone = ui_interp(r_d2,z_d2)
    Bfield_zone[:,0] = Bx_interp(r_d2,z_d2)
    Bfield_zone[:,1] = By_interp(r_d2,z_d2)
    Bfield_zone[:,2] = Bz_interp(r_d2,z_d2)

    # b_zone = 3D unit vector of magnetic field in R-Z plane
    b_zone = np.zeros([Ntri,3])
    for itri in range(0,Ntri):
        nodes = connections[itri,:]
        Bnorm = np.linalg.norm(Bfield_zone[itri,:])
        b_zone[itri,:] = Bfield_zone[itri,:]/Bnorm

    # Source is specified by stratum (Ntri+1) and segment
    # Use the outer polygons and the single stratum for them
    # First segment is the small one on the inboard side.
    # Segments go around COUNTER-clockwise, starting with 
    # index 0 of wallnodes_ordered.
    # Challenge is to find which triangle has this segment 
    # (and there should always be only one).
    # This subroutine gets passed a list of triangles which
    # have at least two wall nodes so that this search is 
    # minimized.
    # TODO: DEGAS2 has its own sheath model. Careful so we don't want to use
    # the source directly from XGC, else it would get doubly-accelerated.
    seg_to_tri=np.zeros(Nwall,dtype=int)
    area_norm=np.zeros(Nwall)
    area_tot=np.zeros(Nwall)
    theta=np.zeros(Nwall)
    # Ad-hoc "center point" from which to calculate theta for diagnostic purposes
    Rcenter = np.average(r_d2)


    # Go around the limiter and calculate the source from each wall segment
    for i in range(0,Nwall):
        node1 = wallnodes_ordered[i]
        node2 = wallnodes_ordered[np.mod(i+1,Nwall)]

        # Vector connecting nodes in a segment
        diff = np.array(node2.coords) - np.array(node1.coords)
        seglen = np.linalg.norm(diff)
        centerpt = 0.5*(np.array(node2.coords) + np.array(node1.coords))
        # Area of this segment
        area_tot[i] = seglen*2.0*np.pi*centerpt[0]

        # Unit vector normal to segment
        a_unit = [-diff[1],0.0,diff[0]]
        a_unit = a_unit/np.linalg.norm(a_unit)

        # Just for outputting source vs theta
        vec = centerpt - [Rcenter,0.0]
        theta[i] = np.arctan2(vec[1],vec[0])

        # Find the zone whose center is closest to the centerpoint of this segment
        dists = np.linalg.norm(zone_coords[:]-centerpt,axis=1)
        tri_id = dists.argmin() 

        # Mapping between segment number and mesh element index
        seg_to_tri[i] = tri_id

        # The "parallel" area that the particle stream sees
        area_norm[i] = area_tot[i]*np.abs(np.dot(a_unit,b_zone[tri_id,:]))

    # Write the plasma zone data
    # Here, ui is the scalar parallel flow velocity, so b unit vector is needed to
    # construct flow velocity in DEGAS2 coordinates.
    defineback.write_plasmafile(ne_zone,Te_zone,Ti_zone,ui_zone=ui_zone,b=b_zone)

    ui_zone = np.abs(ui_zone)

    # Enforce Bohm criterion?
    cs = np.sqrt(1.602e-19*Te_zone/ionmass)
    vti = np.sqrt(2.0*1.602e-19*Ti_zone/ionmass)
    ui_cutoff = 0.5*(ui_zone*(1.0+sp.erf(ui_zone/vti)) + vti*np.exp(-(ui_zone/vti)**2)/np.sqrt(np.pi))
#    ui_cutoff = 0.5*vti/np.sqrt(np.pi)

    ui_zone_use = ui_zone
#    ui_zone_use = np.maximum(ui_zone,cs)
#    ui_zone_use = np.maximum(ui_zone,ui_cutoff)

    source_strength = np.zeros(Nwall)
    # Source strength calculation.
    # Make sure specify_current is set in defienback input file
    if use_xgc_source:
        for i in range(0,Nwall):
            # Since neu_weight_wall_lost is accumulated over a neutral background 
            # step, divide by the appropriate dt here. 
            source_strength[i] = raw_source[wallnode_order[i]]/dt
    else:
        source_strength_xgc = np.zeros(Nwall)
        source_flux = np.zeros(Nwall)
        source_totflux = np.zeros(Nwall)
        ne_wall = np.zeros(Nwall)
        upar_wall = np.zeros(Nwall)
        upar_wall_raw = np.zeros(Nwall)
        cs_wall = np.zeros(Nwall)
        wallflux = np.zeros(Nwall)
        for iwall in range(0,Nwall):
            itri = seg_to_tri[iwall]
#            source_strength[iwall] = area_norm[iwall]*ne_zone[itri]*np.abs(cs[itri])
            wallflux[iwall] = ne_zone[itri]*np.abs(ui_zone[itri])*area_norm[iwall]/area_tot[iwall]
            ne_wall[iwall] = ne_zone[itri]
            cs_wall[iwall] = cs[itri]
            upar_wall_raw[iwall] = np.abs(ui_zone[itri])
            upar_wall[iwall] = np.abs(ui_zone_use[itri])
            source_strength[iwall] = area_norm[iwall]*ne_zone[itri]*np.abs(ui_zone_use[itri])
            source_flux[iwall] = (area_norm[iwall]/area_tot[iwall])*ne_zone[itri]*np.abs(ui_zone_use[itri])
#            source_strength_xgc[iwall] = raw_source[wallnode_order[iwall]]/dt

        plt.plot(theta*180.0/np.pi,ne_wall)
        plt.xlabel("Theta (deg.)")
        plt.ylabel("n_e (wall)")
        plt.tight_layout()
        plt.savefig("newall.pdf")
        plt.clf()
    
        plt.plot(theta*180.0/np.pi,upar_wall_raw)
        plt.plot(theta*180.0/np.pi,cs_wall)
        plt.plot(theta*180.0/np.pi,upar_wall)
        plt.legend(["i_u_par",r"$c_s$",r"max(i_u_par,$\int_0^\infty v f_M dv$)"])
        plt.xlabel("Theta (deg.)")
        plt.ylabel(r"$u_{\|,i}$ @ wall)")
        plt.tight_layout()
        plt.savefig("uiwall.pdf")
        plt.clf()
    
        plt.plot(theta*180.0/np.pi,area_tot)
        plt.xlabel("Theta (deg.)")
        plt.ylabel("Area = 2*pi*R*l_i (wall)")
        plt.tight_layout()
        plt.savefig("areawall.pdf")
        plt.clf()
    
        plt.plot(theta*180.0/np.pi,area_norm/area_tot)
        plt.xlabel("Theta (deg.)")
        plt.ylabel("cos(alpha)")
        plt.tight_layout()
        plt.savefig("cosawall.pdf")
        plt.clf()

        plt.plot(theta*180.0/np.pi,source_strength_xgc)
        plt.plot(theta*180.0/np.pi,source_strength)
        plt.legend(["neu_weight_wall_lost","inferred source"])
        plt.xlabel("Theta (deg.)")
        plt.ylabel("cos(alpha)")
        plt.tight_layout()
        plt.savefig("sourcecomp.pdf")
        plt.clf()

#        plt.semilogy(theta*180.0/np.pi,source_strength_xgc)
#        plt.semilogy(theta*180.0/np.pi,source_strength)
#        plt.legend(["neu_weight_wall_lost","inferred source"])
#        plt.xlabel("Theta (deg.)")
#        plt.ylabel("cos(alpha)")
#        plt.tight_layout()
#        plt.savefig("logsourcecomp.pdf")
#        plt.clf()



    # If using specify_flux, use this instead
    #source_strength = wall_source / (dt*area_tot)
    

    plt.plot(theta*180.0/np.pi,source_strength,"o")
    plt.xlabel("Theta (deg.)")
    plt.ylabel("Neutral source (particles / second)")
    plt.tight_layout()
    plt.savefig("sourceVStheta.pdf")
    plt.clf()

#    plt.plot(theta*180.0/np.pi,source_flux)
#    plt.xlabel("Theta (deg.)")
#    plt.ylabel("Neutral source (particles / m^2 / second)")
#    plt.tight_layout()
#    plt.savefig("fluxVStheta.pdf")
#    plt.clf()
    
    # Write the file while specifies the source for DEGAS2
    strata = np.array([2]*Nwall)
    segments = np.array(range(0,Nwall),dtype=int)

    # Select only those segments which have nonzero source.
    idx = source_strength > 0.0
    source_strength=source_strength[idx]
    segments=segments[idx]
    strata=strata[idx]
    defineback.generate_sourcefile(strata,segments,source_strength)
#    defineback.generate_sourcefile(strata,segments,source_flux)

    return source_strength
#    return source_flux

def write_dummy_bg_files(nzone,nwall_input ,nebins=-1):
    if (nebins > 0):
        nwallsegs = nwall_input*nebins
    else:
        nwallsegs = nwall_input
    nperline = 10
    stratum=nzone+2
    
    f=open("plasmafile.txt","w")
    f.write("%10s %10s %10s %10s %10s\n"%("zone","T(1)","N(1)","T(2)","N(2)"))
    for i in range(0,nzone):
        f.write("%i %e %e %e %e\n"%(i+1,100.0,1e19,100.0,1e19))
    f.close()
    
    f=open("sourcefile.txt","w")
    f.write("#\nstratum\n#\n")
    for i in range(0,nwallsegs):
        f.write("%d "%(stratum))
        if ((i+1)%nperline == 0) and i != nwallsegs-1:
            f.write("\n")
    f.write("\n")
    f.write("#\nsegment\n#\n")
    for i in range(0,nwallsegs):
        if nebins > 0:
            f.write("%d "%(int(i/nwall_input)))
        else:
            f.write("%d "%(i))
        if ((i+1)%nperline == 0) and i != nwallsegs-1:
            f.write("\n")
    f.write("\n")
    f.write("#\nF\n#\n")
    for i in range(0,nwallsegs):
        f.write("%e "%(1e20))
        if ((i+1)%nperline == 0) and i != nwallsegs-1:
            f.write("\n")

def write_dummy_bg_files_aux(nzone,nwall_input,stratum_start):
    nwallsegs = nwall_input
    nperline = 10
    stratum=nzone+2
    
    f=open("plasmafile.txt","w")
    f.write("%10s %10s %10s %10s %10s\n"%("zone","T(1)","N(1)","T(2)","N(2)"))
    for i in range(0,nzone):
        f.write("%i %e %e %e %e\n"%(i+1,100.0,1e19,100.0,1e19))
    f.close()
    
    f=open("sourcefile.txt","w")
    f.write("#\nstratum\n#\n")
    for i in range(0,nwallsegs):
        f.write("%d "%(stratum_start+i))
        if ((i+1)%nperline == 0) and i != nwallsegs-1:
            f.write("\n")
    f.write("\n")
    f.write("#\nsegment\n#\n")
    for i in range(0,nwallsegs):
        f.write("%d "%(0))
        if ((i+1)%nperline == 0) and i != nwallsegs-1:
            f.write("\n")
    f.write("\n")
    f.write("#\nF\n#\n")
    for i in range(0,nwallsegs):
        f.write("%e "%(1e20))
        if ((i+1)%nperline == 0) and i != nwallsegs-1:
            f.write("\n")

def get_xgc_energy_change():
    # To define:
    bk_nc = nc.Dataset("background.nc","r")
    ne_zone = np.array(bk_nc["background_n"][:,0])
    Te_zone_ev = np.array(bk_nc["background_temp"][:,0])/1.602e-19
    Nzone = len(ne_zone)

    nn = postprocess.get_output("neutral density")[0:Nzone,1]

    E_iz = 30.0*1.602e-19

    xgc_izrate = 0.8e-8*ne_zone*np.sqrt(Te_zone_ev)*(np.exp(-np.minimum(13.56/Te_zone_ev, 1.0e6))/(1.0 + 0.01*Te_zone_ev))*1.0e-6

    return xgc_izrate*E_iz*nn








