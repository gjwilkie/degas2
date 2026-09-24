import matplotlib.pyplot as plt
import numpy as np
from polygon import *
import sys
from importlib import reload
import geomutils
import scipy.interpolate as interpolate
import netCDF4 as nc
import matplotlib.tri as mtri
import subprocess
#from tqdm import tqdm

class DG2D:
    """ 
    Class representing all input for definegeometry2d.
    Conventions: First N=len(polys) zones are the internal plasma and vacuum zones, bounded by Nwall=len(wallpoly.vertices). 
    Strata are labelled sequentially starting from 1. Next Nwall zones (and strata) are the embedded auxiliary zones.
    Segments around closed wall polygon are indexed corresponding to the first vertex of the segment when traversing clockwise.

    Steps to building the appropriate data:
    1a. Define the limiter polygon. Triangulation will be done automatically. Option to specify minimum wall triangle size.
    OR
    1b. Import triangular mesh. Limiter shape will be inferred.
    2. Specify wall properties: material, temperature, recycling coefficient. Can be specified by segment.
    3. Write files

    Attributes:
        nodes: list of Vertices, each representing an entry in the wall file from which to build the geometry
        polys: list of Polygons, each representing a plasma or vacuum zone, densely covering a single polygon (see wallnodes)
        wallpoly: Polygon representing the boundary wall, vertices ordered and enforced to be clockwise.
        outerpoly: Polygon representing the outermost polygon embedded in the wall.
        materials: list of strings of length equal to wallpoly.vertices corresponding to the material of that wall segment
        walltemps: list of floats of length equal to wallpoly.vertices corresponding to the wall temeprature of that wall segment in Kelvin
        Rcoeffs: list of floats of length equal to wallpoly.vertices corresponding to the recycling coefficient of that wall segment
        exits: list of booleans of length equal to wallpoly.vertices corresponding whether this wall segment is an exit
        symmetry: string for the type of symmetry in problem
        Rbounds: two-element array for minimum and maximum R in the universal cell
        Zbounds: two-element array for minimum and maximum Z in the universal cell
        vertex_list: "global" list of Vertices to be included in wallfile.
        infile_name: string for the name of the general input file. Deafults to dg2d.in and likely to not change.
        wallfile_name: string for the name of the "wallfile". Deafults to dg2d.in and likely to not change.
    """

    def __init__(self):
        """
        Constructor for DG2D object. Generates for minimal defaults.
        """
        self.polys = []
        self.wallpoly = None
        self.outerpoly = None
        self.materials = None
        self.walltemps = None
        self.Rcoeffs = None
        self.write_polygonfile = True
        self.symmetry = "cylindrical"
        self.Rbounds = [0.01,10.0]
        self.Zbounds = [-10.0,10.0]
        self.vertex_list = []
        self.current_stratum = 0
        self.infile_name = "dg2d.in"
        self.wallfile_name = "wallfile.txt"
        self.polygonfile_name = "polygon.nc"

    def finalize(self):
        subprocess.run("definegeometry2d "+self.infile_name,shell=True)

    def set_symmetry(self,sym):
        self.symmetry = sym
   
    def set_wallprops(self,walltemp=300.0,Rcoeff=1.0,material="mirror",exitzone=False,wallidx=None):
        """
        Sets boundary properties for domain. Can be applied uniformly or by segment depending on if wallidx is provided.

        Args:
            walltemp: (optional) wall temperature in Kelvin. Default = 300.
            Rcoeff: (optional) float for recycling coefficient. Default = 1.0.
            material: (optional) string for the wall material name. Default = "mirror"
            exitzone: (optional) boolean default to False. If True, segment is a exit stratum rather than a plasma-facing component.
            wallidx: (optional) integer index of the wall segment to apply properties to. If not provided, will apply the provided properties around entire wall.
        """
        # Apply to all segments as a default
        if wallidx == None:
            self.walltemps= np.array([walltemp]*len(self.wallpoly.vertices))
            self.Rcoeffs= np.array([Rcoeff]*len(self.wallpoly.vertices))
            self.materials= [material]*len(self.wallpoly.vertices)
            self.exits= [exitzone]*len(self.wallpoly.vertices)
        else:
            self.walltemps[wallidx]= walltemp
            self.Rcoeffs[wallidx]= Rcoeff
            self.materials[wallidx]= material
            self.exits[wallidx]= exitzone

    def define_mesh(self,coords,conn,wallnodes=[-1],trust_wallnodes=True,wallpolys=[-1],progress=False):
        """
        Imports a set of coordinates and connectivity to generate a mesh for definegeometry2d.

        Args:
            coords: float array of size (Nnode, 2), where the first and second columns are the R and Z coordinates of the nodes, respectively.
            conn: integer array of size (Npoly, Nseg), where each row are the node indices (in coords) that connect to form a polygon (usually a triangle). Must densely cover a domain.
            wallnodes: (optional) integer array that specifies which indices of the nodes (in coords) form the boundary. Optional for triangular meshes; mandatory for higher-order polygons.
        """
        self.polys=[]
        Nnode = len(coords[:,0])
        Nseg = len(conn[0,:])
        Npoly = len(conn[:,0])

        for i in range(0,Nnode):
            self.vertex_list.append(Vertex(i,coords[i,0],coords[i,1])) 
        
        iterator = range(0,Npoly)
        if progress:
            iterator = tqdm(iterator)
        for i in iterator:
            poly = Polygon()
            for j in range(0,Nseg):
                if conn[i,j] >= 0:
                    poly.add_vertex(self.vertex_list[conn[i,j]])
            self.polys.append(poly)

        print("Done defining mesh. Now finding wall nodes...")

        self.wallpoly = Polygon(increment=False)
        startidx = get_first_idx(self.vertex_list) 
        if (wallnodes[0] == -1 or not trust_wallnodes):
            if (Nseg > 3):
                print("ERROR: for quadrilateral or higher basic polygons in call to define_mesh, wallnodes must be specified")
                
            wallnodes_out=[]
            wallnodes_out.append(self.vertex_list[startidx])
            self.wallpoly.add_vertex(wallnodes_out[-1])
            closed = False
            prevnode=wallnodes_out[-1]
            first = True
            while not closed:
                if progress:
                    print(len(wallnodes_out))
                if wallpolys == [-1]:
                    next_node = find_next_wall_node(wallnodes_out[-1],prevnode,self.vertex_list,self.polys,first,check_clockwise=True)
                else:
                    next_node = find_next_wall_node(wallnodes_out[-1],prevnode,wallnodes,wallpolys,first,check_clockwise=True)
                first = False
                if next_node == wallnodes_out[0]:
                    closed = True
                else:
                    prevnode = wallnodes_out[-1]
                    wallnodes_out.append(next_node)
                    self.wallpoly.add_vertex(wallnodes_out[-1])
            self.wallpoly.is_clockwise(force=True)
        else:
            Nwall = len(wallnodes)
            for i in range(0,Nwall):
                self.wallpoly.add_vertex(self.vertex_list[wallnodes[np.mod(startidx+i,Nwall)]])

    def define_limiter(self,Rlim_in,Zlim_in,maxdist=-1):
        """
        Defines the boundary of a domain to be triangulated by definegeometry2d.

        Args:
            Rlim: array of floats specifying R coordinates of boundary nodes.
            Zlim: array of floats specifying Z coordinates of boundary nodes.
        """

        limpoly = Polygon()
        if len(Rlim_in) != len(Zlim_in):
            print("ERROR: Rlim and Zlim must have the same dimensions in call to define_limiter")
        Nlim = len(Rlim_in)

        mindist = 999999
        minidx = -1
        for i in range(0,len(Rlim_in)):
            dist = np.linalg.norm([Rlim_in[i],Zlim_in[i]])
            if dist < mindist:
                minidx = i
        
        Rlim = Rlim_in[np.mod( np.array(range(minidx,minidx+Nlim),dtype=int),Nlim)]
        Zlim = Zlim_in[np.mod( np.array(range(minidx,minidx+Nlim),dtype=int),Nlim)]
        if maxdist > 0:
            Rlim, Zlim = refine_limiter(Rlim,Zlim,maxdist)

        for i in range(0,Nlim):
            self.vertex_list.append(Vertex(i,Rlim[i],Zlim[i]))
            limpoly.add_vertex(self.vertex_list[-1])
        limpoly.is_clockwise(force=True)
        limpoly.split = True
        self.wallpoly= limpoly
        self.polys.append(limpoly)

    def write_files(self,aux_thickness=0.005):
        """
        Writes definegeometry2d input files from data in DG2D object.
        All data must be specified except for auxilliary polygons.
        """

        def write_internal_polygon(poly,newzone=True):
            f = open(self.infile_name,"a")
    
            if newzone:
                f.write("new_zone ")
            if poly.vacuum:
                f.write("vacuum\n")
            else:
                f.write("plasma\n")
    
            f.write("new_polygon\n")
            self.current_stratum += 1
            f.write("  stratum "+str(self.current_stratum)+"\n")
            if poly.is_clockwise():
                for vertex in poly.vertices:
                    f.write("  wall 1 "+str(vertex.id)+" "+str(vertex.id)+"\n")
            else:
                for vertex in reversed(poly.vertices):
                    f.write("  wall 1 "+str(vertex.id)+" "+str(vertex.id)+"\n")
    
            if poly.split:
                f.write("  triangulate_to_zones\n")
            else:
                f.write("  triangulate_polygon\n")
            f.write("\n")
            f.close()
    
        def write_aux_polygon(poly,material=None,walltemp=300.0,Rcoeff=1.0,exitzone=False):
            f = open(self.infile_name,"a")
    
            if exitzone:
                f.write("new_zone exit\n")
            else:
                f.write("new_zone solid\n")
    
            f.write("new_polygon\n")
            self.current_stratum += 1
            f.write("  stratum "+str(self.current_stratum)+"\n")
            if not exitzone:
                f.write("  material "+material+"\n")
                f.write("  temperature "+str(walltemp)+"\n")
                f.write("  recyc_coef "+str(Rcoeff)+"\n")
            if poly.is_clockwise():
                for vertex in poly.vertices:
                    f.write("  wall 1 "+str(vertex.id)+" "+str(vertex.id)+"\n")
            else:
                for vertex in reversed(poly.vertices):
                    f.write("  wall 1 "+str(vertex.id)+" "+str(vertex.id)+"\n")
    
            f.write("  triangulate_polygon\n")
            f.write("\n")
            f.close()


        # Calculate bounds
        Nvertex = len(self.vertex_list)
        Rmin = 9e30
        Rmax = -9e30
        Zmin = 9e30
        Zmax = -9e30

        for i in range(0,Nvertex):
            if self.vertex_list[i].coords[0] < Rmin:
                Rmin = self.vertex_list[i].coords[0]
            if self.vertex_list[i].coords[0] > Rmax:
                Rmax = self.vertex_list[i].coords[0]
            if self.vertex_list[i].coords[1] < Zmin:
                Zmin = self.vertex_list[i].coords[1]
            if self.vertex_list[i].coords[1] > Zmax:
                Zmax = self.vertex_list[i].coords[1]
    
        dR = Rmax-Rmin
        dZ = Zmax-Zmin
        if self.symmetry == "cylindrical":
            self.Rbounds = [max(0.005,0.5*Rmin),Rmax+0.5*dR]
            self.Zbounds = [Zmin - 0.5*dZ, Zmax + 0.5*dZ]
        else:
            self.Rbounds = [Rmin-0.5*dR,Rmax+0.5*dR]
            self.Zbounds = [Zmin - 0.5*dZ, Zmax + 0.5*dZ]


        # Write header
        f = open(self.infile_name,"w")
        f.write("symmetry "+self.symmetry+"\n")
        f.write("bounds %f %f  %f %f\n"%(self.Rbounds[0],self.Rbounds[1],self.Zbounds[0],self.Zbounds[1]))
        f.write("wallfile "+self.wallfile_name+"\n")
        f.write("end_prep\n")
        f.write("\n")
        f.close()

        for poly in self.polys:
            write_internal_polygon(poly)
        
        if True:
            aux_polys, outpoly, newvertices = self.wallpoly.build_aux_wall_polygons(self.vertex_list[-1].id+1,aux_thickness)
            self.outerpoly = outpoly

            Nnew = len(newvertices)
            for i in range(0,Nnew):
                self.vertex_list.append(newvertices[i])

            for i in range(0,len(aux_polys)):
                write_aux_polygon(aux_polys[i],material=self.materials[i],walltemp=self.walltemps[i],Rcoeff=self.Rcoeffs[i],exitzone=self.exits[i])
        else:
            outpoly = self.wallpoly
            self.outerpoly = self.wallpoly

        f = open(self.infile_name,"a") 
        Polygon.close_in_universal_cell(f,outpoly.vertices,0,self.current_stratum+1,self.materials[-1],self.Rcoeffs[-1],walltemp=self.walltemps[-1],clockwise=outpoly.is_clockwise())

        if self.write_polygonfile:
            f.write("polygon_nc_file "+self.polygonfile_name+"\n")
        else:
            f.write("polygon_nc_file none\n")
        f.write("end\n")
        f.close()

        Nnode = len(self.vertex_list)
        wallfile = open(self.wallfile_name,"w")
        wallfile.write("# Wallfile automatically generated by script\n")
        wallfile.write("1 \n")
        wallfile.write("#\n")
        wallfile.write("%d\n"%(Nnode))
        wallfile.write("#\n")
        for i in range(0,Nnode):
            wallfile.write("%26.20f   %26.20f \n"%(self.vertex_list[i].coords[0],self.vertex_list[i].coords[1]))
        wallfile.close()



def setup(material,Rcoeff,Rlim=None,Zlim=None,gfile=None,bpfile=None,triang=None,walltemp=300.0,bindir="",run_dg2d=True,polygonfile_name="polygon.nc"):
    """
    Consolidated workflow routine to generate geometry for typical cases and run definegeometry2d. Generates in one of several ways depending on the keyword arguments provided.

    Args:
        material: string for the wall material (to be used uniformly around boundary)    
        Rcoeff: float for the wall recycling coefficient (to be used uniformly around boundary)
        walltemp: (optional) float for the wall temperature (to be used uniformly around boundary) in Kelvin. Default 300.
        Rlim, Zlim: (optional) arrays of floats that define the limiter. Can be clockwise or counter-clockwise. If provided, definegeometry2d itself will perform triangulation and no other data is needed.
        gfile: (optional) string for path and name to a EQDSK "g" file. Used here only to extract the limiter shape. Behvaior follows as if Rlim and Zlim were provided manually.
        bpfile: (optional) string for path and name to an ADIOS2 file that contains a mesh which will be mimiced in the definegeometry2d input files with each triangle being a unique zone.
        triang: (optional) a matlotlib.tri.Triangulation object whose triangulation will be imported and each triangle will be its own zone in the same order defined.
        run_dg2d: (optional) boolean for whether definegeometry2d is to be run again. Default True. 
        polygonfile_name: (optional) string for the polygon filename. Defaults to "polygon.nc". Set to "none" for particularly large (>10k node) meshes.
    """

    self = DG2D()

    coords = [0]
    if Rlim != None:
        self.define_limiter(Rlim,Zlim)
    elif gfile != None:
        from geomutils import read_geqdsk
        g = read_geqdsk(gfile)
        Rlim = g.lim[:,0]
        Zlim = g.lim[:,1]
        self.define_limiter(Rlim,Zlim)
    elif bpfile != None:
        coords,conn,wallnodes = get_bp_mesh(bpfile) 
        self.define_mesh(coords,conn,wallnodes=wallnodes)
    elif triang != None:
        r = triang.x
        z = triang.y
        nnode = len(r)
        coords = np.zeros((nnode,2))
        coords[:,0] = r
        coords[:,1] = z
        conn = triang.triangles
        self.define_mesh(coords,conn)
    else:
        print("ERROR: calculate method requires some geometry data")

    if (material == None) or (Rcoeff == None):
        print("ERROR: calculate method requires material and recycling coefficient")

    self.set_wallprops(walltemp=walltemp,material=material,Rcoeff=Rcoeff)
    if polygonfile_name == "polygon.nc" and np.shape(coords)[0] > 10000:
        print("WARNING: definegeometry2d is instructed to write polygon_nc_file even though mesh is very large. May take too long.\n")
    self.polygonfile_name = polygonfile_name

    self.write_files()
    if bindir != "" and bindir[-1] != '/':
        bindir = bindir+"/"

    if run_dg2d:
        subprocess.run(bindir+"definegeometry2d dg2d.in",shell=True)


def get_first_idx(vertex_list):
    # Start with vertex whose corresponding segment is closest to origin
    startidx = -1
    Nvertex = len(vertex_list)
    mindist = 9.0e30
    for i in range(0,Nvertex):
        segdist = np.linalg.norm(0.5*np.array(vertex_list[i].coords + vertex_list[np.mod(i+1,Nvertex)].coords))
        if segdist < mindist:
            mindist = segdist 
            startidx = i
    return startidx


def write_dg2d_header(f,symmetry,Xmin,Xmax,Zmin,Zmax,wallfile_name="wallfile.txt"):
    f.write("symmetry "+symmetry+"\n")
    f.write("bounds %f %f  %f %f\n"%(Xmin,Xmax,Zmin,Zmax))
    f.write("wallfile "+wallfile_name+"\n")
    f.write("end_prep\n")
    f.write("\n")

# Not sure this works. Use poly_to_try executable instead.
def get_triangulation(tris,nodes):
    Nnode = len(nodes)
    r = np.zeros(Nnode)
    z = np.zeros(Nnode)
    for i in range(0,Nnode):
        r[i] = nodes[i].coords[0]
        z[i] = nodes[i].coords[1]
    Ntri = len(tris)
    conn = np.zeros((Ntri,3),dtype=int)
    for i in range(0,Ntri):
        if isclockwise(tris[i]):
            conn[i,0] = tris[i].vertices[2].id
            conn[i,1] = tris[i].vertices[1].id
            conn[i,2] = tris[i].vertices[0].id
        else:
            conn[i,0] = tris[i].vertices[0].id
            conn[i,1] = tris[i].vertices[1].id
            conn[i,2] = tris[i].vertices[2].id
    return mtri.Triangulation(r,z,triangles=conn)

def get_triangulation_from_polygons(polyfile,geomfile="geometry.nc"):
    ncdata = nc.Dataset(polyfile)
    coords = ncdata["g2_polygon_xz"]
    coords = np.array(coords)

    ncdata = nc.Dataset(geomfile)
    zone_type = ncdata["zone_type"]
    Nzone = 0
    # Count number of plasma zones
    for i in range(0,len(zone_type)):
        if zone_type[i] == 2:
            Nzone+=1
    # Assumes all non-plasma zones are the end

    # Searches through array of vertices,looks for matching coordinates,
    # and returns index. If not found, returns -1
    def find_vertex(coords,vertices):
        eps = 1.0e-8
        idx = -1
        for i in range(0,len(vertices[:,0])):
            if (np.abs(coords[0] - vertices[i,0]) < eps) and (np.abs(coords[1] - vertices[i,1]) < eps):
                idx = i
        return idx


    # Initially populate array of vertices
    allvertices=coords[0,0:3,0:2]
    allconn = []
    # Go through polygons and record vertices
    for i in range(0,Nzone):
        vertices=coords[i,0:3,0:2]
        conn = np.zeros(3,dtype=int)
        for j in range(0,3):
            idx = find_vertex(vertices[j,:], allvertices)
            if idx == -1:
                allvertices = np.append(allvertices,[vertices[j,:]],axis=0)
                idx = len(allvertices[:,0]) - 1
            conn[j] = idx
        allconn.append(conn)
    allconn = np.array(allconn,dtype=int)

    return mtri.Triangulation(allvertices[:,0],allvertices[:,1],triangles=allconn), len(allconn[:,0])

def write_cylinder_dg2d_input(Rmax,NR,wallfile_name="wallfile.txt",dg2dfile_name="dg2d.in",recyc=1.0,walltemp=300.0,Rmin=1.0e-3,Zfac=10.0,material="C"):
    Xmin=0.5*Rmin
    Xmax=1.5*Rmax
    Zmin_tot=-Rmax
    Zmax_tot=(Zfac+1)*Rmax
    Zmin = 0.0
    Zmax = Zfac*Rmax

    dr = Rmax/(NR-0.5)
    rgrid=np.linspace(0.5*dr,Rmax,NR)
    rgrid = np.insert(rgrid,0,0.0)
    rgrid_ctr = np.linspace(0.0,Rmax-0.5*dr,NR)

    # Write wallfile
    wallfile = open(wallfile_name,"w")

    wallfile.write("# Wallfile automatically generated by script\n")
    wallfile.write("%d\n"%(NR+1))
    wallfile.write("#\n")
    for ir in range(0,NR+1):
        wallfile.write("2 ")
    wallfile.write("\n#\n")
    
    # First wall is at r=Rmin, not 0
    wallfile.write("%e %e\n"%(Rmin,Zmin))
    wallfile.write("%e %e\n"%(Rmin,Zmax))
    for ir in range(0,NR):
        wallfile.write("%e %e\n"%((ir+0.5)*dr,Zmin))
        wallfile.write("%e %e\n"%((ir+0.5)*dr,Zmax))
    wallfile.close()

    dfile=open(dg2dfile_name,"w")
    write_dg2d_header(dfile,"cylindrical",Xmin,Xmax,Zmin_tot,Zmax_tot,wallfile_name=wallfile_name)

    # Zones 1 through NR are the plasma zones of interest. Solid zones come later.

    for ir in range(0,NR):
        dfile.write("new_zone plasma\n")
        dfile.write("new_polygon\n")
        dfile.write("  stratum %d\n"%(ir+1))
        dfile.write("  wall %d 0 1\n"%(ir+1))
        dfile.write("  wall %d 1 0\n"%(ir+2))
        dfile.write("  triangulate_polygon\n\n")

    dfile.write("new_zone solid\n")
    dfile.write("new_polygon\n")
    dfile.write("  stratum %d\n"%(NR+1))
    dfile.write("  material %s\n"%material)
    dfile.write("  recyc_coef %f\n"%recyc)
    dfile.write("  temperature %f\n"%walltemp)
    dfile.write("  wall %d 0 1\n"%(NR+1))
    dfile.write("  outer 2 3\n")
    dfile.write("  breakup_polygon\n\n")

    dfile.write("new_zone solid\n")
    dfile.write("new_polygon\n")
    dfile.write("  stratum %d\n"%(NR+2))
    dfile.write("  material mirror\n")
    dfile.write("  recyc_coef 1.0\n")
    dfile.write("  outer 0 1 2\n")
    for ir in range(0,NR+1):
        dfile.write("  wall %d 1 1\n"%(NR+1-ir))
    for ir in range(0,NR+1):
        dfile.write("  wall %d 0 0\n"%(ir+1))
    dfile.write("  outer 3 4\n")
    dfile.write("  breakup_polygon\n\n")

    dfile.write("polygon_nc_file polygon.nc\nend\n")
    dfile.close

    return rgrid_ctr, NR+1

    

def write_cylindrical_polygon_dg2d_input(R_tot,NR,material,Ntheta_min=12,Ntheta_max=200,wallfile_name="wallfile.txt",dg2dfile_name="dg2d.in",recyc=1.0,walltemp=300.0):
    Xmin = -1.5*R_tot
    Xmax = 1.5*R_tot
    Zmin = -1.5*R_tot
    Zmax = 1.5*R_tot

    delta_r = 2*R_tot/(2*NR-1.0)
    r_grid = np.linspace(0.0,R_tot-0.5*delta_r,num=NR)

    plasma_polys = []

    # First polygon (0) is simple: a polygon with Ntheta_min sides. 
    #    Corresponds to the center of the plasma  at rgrid[0]=0.0 
    # Next (1:NR) polygons are the "main" plasma ones. Consists of a cresent-shaped polygon,
    #    leaving a small quadrilateral, which gets defined next. This is necessary because
    #    polygons with holes is not supported.
    # Next (NR:2*NR-1) are the "fillers". Statistics from these will not be reported.
    # Finally, (2*NR:2*NR+1) are the wall polygons

    # ith polygon has Ntheta_min*(r_grid[i]+0.5*delta_r)/(0.5*delta_r) nodes sides
    Ntheta = ( Ntheta_min*(r_grid+0.5*delta_r)/(0.5*delta_r)).astype(int)
    Ntheta = np.minimum(Ntheta,Ntheta_max)


    # Gets the wall id (counting from 1) of a node based on its global index (counting from 0)
    def get_node_wall_id(node_id):
        wall_id = -1

        n_node = 0
        for i in range(0,NR):
            if node_id >= n_node and node_id < n_node+Ntheta[i]:
                wall_id = i+1
            n_node += Ntheta[i]

        if wall_id < 1:
            sys.exit("ERROR: index algebra failed in get_node_wall_id; node_id = $d, wall_id = $d"%(node_id,wall_id))
        return wall_id

    def get_node_wall_ids_from_global_index(i):
        wall_id = -1
        node_id = -1

        n_node = 0
        for i in range(0,NR):
            if i >= n_node and i < n_node+Ntheta[i]:
                wall_id = i+1
                node_id = i - n_node
            n_node += Ntheta[i]

        if wall_id < 1 or node_id < 0:
            sys.exit("ERROR: index algebra failed in get_node_wall_ids_from_global_index; i = $d"%(i))
        return wall_id, node_id

    ####################################################
    # Define surfaces (for "wallfile")
    surfs = []
    for i in range(0,NR):
        r = delta_r * (i + 0.5)
        dtheta = 2.0*np.pi/Ntheta[i]
        theta = np.linspace(0.0,2.0*np.pi,num=Ntheta[i],endpoint=False)
        
        # By convention (which is relied upon later in this function), 
        #  start with the leftmost point and go clockwise
        x = -r*np.cos(theta)
        z = r*np.sin(theta)

        # In dg2d, walls are counted from 1, but nodes within them are counted from 0
        surfs.append(Surface(i+1))

        inode = 0
        # Now build the nodes
        for j in range(0,Ntheta[i]):
            node = Vertex(inode,x[j],z[j])
            inode+=1
            surfs[-1].add_vertex(node)
            surfs[-1].vertices[-1].wall_id = i+1

    ####################################################
    # Write wallfile using the defined surfaces
    wallfile = open(wallfile_name,"w")
    wallfile.write("# Wallfile automatically generated by script\n")
    wallfile.write("%d\n"%(NR))
    wallfile.write("#\n")
    for i in range(0,NR):
        wallfile.write("%d "%Ntheta[i])
    wallfile.write("\n")
    wallfile.write("# \n")
    for i in range(0,NR):
        wallfile.write("# coordinates for wall %d - %d points\n"%(i+1,Ntheta[i]))
        for j in range(0,Ntheta[i]):
            x = surfs[i].vertices[j].coords[0]
            z = surfs[i].vertices[j].coords[1]
            wallfile.write("%f %f\n"%(x,z))
    wallfile.close()

    ####################################################
    # Build polygons from wallfile nodes
    polys = []
     
    # First zone is straightforward
    # Next (NR-1) zones consist of two polygons: the "main" cresent-shaped ones, and the small quadrliateral that closes it to an annulus.
    for isurf in range(0,NR):
        polys.append(Polygon())
        for j in range(0,Ntheta[isurf]):
            polys[-1].add_vertex(surfs[isurf].vertices[j])
        # Now go backwards along previous surface
        if isurf > 0:
            for j in range(Ntheta[isurf-1]-1,-1,-1):  # Too many "-1"s!
                polys[-1].add_vertex(surfs[isurf-1].vertices[j])
            # Just as a sanity check, add the first point again. dg2d is fine with this.
            polys[-1].add_vertex(surfs[isurf].vertices[0])

            # Now make the little quadrilaterals that close the "flux surface"
            polys.append(Polygon())
            polys[-1].add_vertex(surfs[isurf].vertices[0])
            polys[-1].add_vertex(surfs[isurf].vertices[-1])
            polys[-1].add_vertex(surfs[isurf-1].vertices[-1])
            polys[-1].add_vertex(surfs[isurf-1].vertices[0])

    ####################################################
    # Write dg2d input file
    dg2dfile = open(dg2dfile_name,"w")

    write_dg2d_header(dg2dfile,"plane",Xmin,Xmax,Zmin,Zmax,wallfile_name=wallfile_name)

    for ipoly in range(0,len(polys)):
        if ipoly < NR:
            wallid = ipoly+1
        else:
            wallid = ipoly - NR + 2

        if (ipoly == 0) or ( (ipoly+1)%2 == 0):
            newzone = True
        else:
            newzone = False

        polys[ipoly].write_plasma_polygon_dg2d(dg2dfile,stratum=ipoly+1,wallid=wallid,commonzone=True,newzone=newzone)

    wallnodes = surfs[-1].vertices

    # Syntax is a little confusing right now. "wallid" in this routine actually assumes 
    # they're counted from zero and corrects to write dg2d input file
    Polygon.close_in_universal_cell(dg2dfile,wallnodes,NR-1,2*NR,material,recyc,clockwise=True,walltemp=walltemp)
    dg2dfile.write("polygon polygon.nc\n")
    dg2dfile.write("end\n")

    dg2dfile.close()

    return r_grid

def write_dg2d_input_from_triangle_file(trifile_base,material,recyc,dg2dfile_name="dg2d.in",polygon_filename="polygons.nc",debug=False,trust_wallflags=True,ionmass=1.67e-27,newformat=False,minarea=0.0):

    def next_noncomment_line(f):
        found=False
        while not found:
            line = f.readline().strip()
            if not (line[0] == '#' or line[0:2] == "//"):
                found = True
        return line

    nodefile = open(trifile_base+".node",'r')
    line1 = next_noncomment_line(nodefile).split(",")
    nnode = int(line1[0])
    nattr = int(line1[2])

    node_coords = np.zeros([nnode,2])
    Te_node = np.zeros([nnode])
    Ti_node = np.zeros([nnode])
    ne_node = np.zeros([nnode])
    mach_node = np.zeros([nnode])
    Br_node = np.zeros([nnode])
    Bt_node = np.zeros([nnode])
    Bz_node = np.zeros([nnode])
    wallflag_node = np.zeros([nnode],dtype=int)
    lowestRnode = 0
    Rmin = -1.0
    Rmax = -1.0
    Zmin = 99999.0
    Zmax = -99999.0
    nnode_read=0
    vertices = []
    wallvertices = []
    while nnode_read < nnode:
        line = next_noncomment_line(nodefile).split(",")
        inode = int(line[0])
        node_coords[inode,0] = float(line[1])
        node_coords[inode,1] = float(line[2])
        Te_node[inode] = 1.602e-19*float(line[3])
        Ti_node[inode] = 1.602e-19*float(line[4])
        ne_node[inode] = float(line[5])
        mach_node[inode] = float(line[6])
        Br_node[inode] = float(line[7])
        Bt_node[inode] = float(line[8])
        Bz_node[inode] = float(line[9])
        if newformat:
#            wallflag_node[inode] = int(line[11])
            wallflag_node[inode] = int(line[13])
        else:
            wallflag_node[inode] = int(line[10])

        vertices.append(Vertex(inode,node_coords[inode,0],node_coords[inode,1]))

        if not trust_wallflags:
            vertices[-1].wall = False  # Using invalid triangles later to add wall nodes
        elif wallflag_node[inode] == 1:
            vertices[-1].wall = True
            wallvertices.append(vertices[-1])
            wallvertices[-1].id = vertices[-1].id
        nnode_read += 1

    nodefile.close()

    elefile = open(trifile_base+".ele",'r')
    line1 = next_noncomment_line(elefile).split(",")
    ntri = int(line1[0])
    nattr = int(line1[2])

    trinodes = np.zeros([ntri,3],dtype=int)
    Te_tri = np.zeros([ntri])
    Ti_tri = np.zeros([ntri])
    ne_tri = np.zeros([ntri])
    mach_tri = np.zeros([ntri])
    Br_tri = np.zeros([ntri])
    Bt_tri = np.zeros([ntri])
    Bz_tri = np.zeros([ntri])
    validflag = np.zeros([ntri],dtype=int)
    lowestRtri = -1
    ntri_read=0
    polys = []
    wallpolys = []
    invalidpolys = []
    ne_zone = []
    Te_zone = []
    Ti_zone = []
    mach_zone = []
    zone_map = []
    iline = 0
    while iline < ntri:
        line = next_noncomment_line(elefile).split(",")
        itri = int(line[0])
        trinodes[itri,0] = int(line[1])
        trinodes[itri,1] = int(line[2])
        trinodes[itri,2] = int(line[3])
        Te_tri[itri] = 1.602e-19*float(line[4])
        Ti_tri[itri] = 1.602e-19*float(line[5])
        ne_tri[itri] = float(line[6])
        mach_tri[itri] = float(line[7])
        Br_tri[itri] = float(line[8])
        Bt_tri[itri] = float(line[9])
        Bz_tri[itri] = float(line[10])
        if newformat:
#            validflag[itri] = int(line[12])
            validflag[itri] = int(line[14])
            area = node_coords[trinodes[itri,0],0]*(node_coords[trinodes[itri,1],1]-node_coords[trinodes[itri,2],1])
            area += node_coords[trinodes[itri,1],0]*(node_coords[trinodes[itri,2],1]-node_coords[trinodes[itri,0],1])
            area += node_coords[trinodes[itri,2],0]*(node_coords[trinodes[itri,0],1]-node_coords[trinodes[itri,1],1])
            area = 0.5*abs(area)
            if area < minarea:
                validflag[itri] = 0

        else:
            validflag[itri] = int(line[11])

        if (validflag[itri] == 0):
            zone_map.append(-1)
            invalidpolys.append(Polygon())
            invalidpolys[-1].id = itri
            invalidpolys[-1].add_vertex(vertices[trinodes[itri,0]])
            invalidpolys[-1].add_vertex(vertices[trinodes[itri,1]])
            invalidpolys[-1].add_vertex(vertices[trinodes[itri,2]])
        else:
            zone_map.append(len(polys))
            polys.append(Polygon())
            polys[-1].id = itri
            polys[-1].add_vertex(vertices[trinodes[itri,0]])
            polys[-1].add_vertex(vertices[trinodes[itri,1]])
            polys[-1].add_vertex(vertices[trinodes[itri,2]])
   
            ne_zone.append(ne_tri[itri])
            Te_zone.append(Te_tri[itri])
            Ti_zone.append(Ti_tri[itri])
            mach_zone.append(mach_tri[itri])
        iline += 1
    elefile.close()
    ntri = len(polys)

    ne_zone = np.array(ne_zone)
    Te_zone = np.array(Te_zone)
    Ti_zone = np.array(Ti_zone)
    mach_zone = np.array(mach_zone)

    if not trust_wallflags:
#        wallvertices = infer_wall_nodes(polys,invalidpolys)
        vertices = purge_invalid_nodes(vertices,polys)
        wallvertices = infer_wall_nodes(polys,invalidpolys)

    # Now that we know the wall vertices for sure, 
    # go back and calculate some needed things:
    rtemp = []
    ztemp = []
    for node in vertices:
        if node in wallvertices:
            node.wall = True

#            print("wallnode id = %d"%node.id)

            if (node.coords[0] < Rmin or Rmin < 0.0) :
                Rmin= node.coords[0]
            if node.coords[0] > Rmax:
                Rmax= node.coords[0]
            if node.coords[1] < Zmin:
                Zmin= node.coords[1]
            if node.coords[1] > Zmax:
                Zmax= node.coords[1]
    
    #plt.plot(rtemp,ztemp,".")
    #plt.plot(lowestRnode.coords[0],lowestRnode.coords[1],"o")
    #plt.savefig("test.png")
    #plt.close()

    Zmin_tot = Zmin - (Zmax-Zmin)
    Zmax_tot = Zmax + (Zmax-Zmin)
    Rmax_tot = 1.5*Rmax
    Rmin_tot = max(0.0001,0.5*Rmin)

    mindist = 1.0e20
    innerSplitGuidePt = np.array([Rmin_tot,0.5*(Zmax_tot+Zmin_tot)])
    for node in vertices:
        for node in wallvertices:
            dist = (innerSplitGuidePt[0]-node.coords[0])**2 + (innerSplitGuidePt[1]-node.coords[1])**2
            if dist < mindist:
                lowestRnode = node
                mindist = dist

    for tri in polys:
        nwallnodes = tri.reorder_wallnodes_first()
        if nwallnodes >= 2:
            tri.alongwall = True
            wallpolys.append(tri)
#            print("wallpoly.id = %d"%wallpolys[-1].id)
            ntest = wallpolys[-1].reorder_wallnodes_first()
            if not nwallnodes == ntest:
                print(nwallnodes)
                print(ntest)
                sys.exit("Wrong number of wall nodes. Something wrong with the logic.")

#    print("Number of wall triangles= %d"%len(wallpolys))

#    # wallpolys contains all polygons that have a segment along the wall
#    # and the wall segment is first
#    wall_strata = []
#    source_strength = []
#    for tri in wallpolys:
#        point1 = np.array(tri.vertices[0].coords)
#        point2 = np.array(tri.vertices[1].coords)
#        diff = point2 - point1
#        area = 2.0*np.pi*0.5*(point1[0]+point2[0])*np.linalg.norm(diff)
#        normal = [-diff[1],0.0,diff[0]]
#        a_unit = normal/np.linalg.norm(normal)
#
#        # TODO: pretty sure this logic is wrong if the nodes are not listed in order
#        Br = 0.5*(Br_node[tri.vertices[0].id] + Br_node[tri.vertices[1].id])
#        Bt = 0.5*(Bt_node[tri.vertices[0].id] + Bt_node[tri.vertices[1].id])
#        Bz = 0.5*(Bz_node[tri.vertices[0].id] + Bz_node[tri.vertices[1].id])
#        b_unit = [Br,Bt,Bz]/np.linalg.norm([Br,Bt,Bz])
#
#        wetarea = area*np.abs(np.dot(b_unit,a_unit))
#
#        cs = 0.5*(np.sqrt(Te_node[[tri.vertices[0].id]]/ionmass) + \
#                np.sqrt(Te_node[[tri.vertices[1].id]]/ionmass))
#
#        ne = 0.5*(ne_node[[tri.vertices[0].id]] + ne_node[[tri.vertices[1].id]])
#
#        source_strength.append(ne*cs*wetarea)
#        wall_strata.append(tri.id+1)

    dg2dfile = open(dg2dfile_name,'w')
    dg2dfile.write("symmetry cylindrical\n")
    dg2dfile.write("bounds "+str(Rmin_tot)+" "+str(Rmax_tot)+" "+\
            str(Zmin_tot)+" "+str(Zmax_tot)+"\n")
    dg2dfile.write("wallfile wallfile.txt\n")
    dg2dfile.write("end_prep\n\n")

    def isclockwise(tri):
        if not len(tri.vertices) == 3:
            sys.exit("Not a triangle isclockwise logic not valid.")
        s = 0.0
        s += (tri.vertices[1].coords[0] - tri.vertices[0].coords[0])*(tri.vertices[1].coords[1] + tri.vertices[0].coords[1])
        s += (tri.vertices[2].coords[0] - tri.vertices[1].coords[0])*(tri.vertices[2].coords[1] + tri.vertices[1].coords[1])
        s += (tri.vertices[0].coords[0] - tri.vertices[2].coords[0])*(tri.vertices[0].coords[1] + tri.vertices[2].coords[1])
        if s >= 0.0:
            return True
        else:
            return False

    for poly in polys:

        if not isclockwise(poly):
            poly.vertices[:] = poly.vertices[-1::-1]

    stratum = 0
    for poly in polys:
        stratum += 1
        poly.write_plasma_polygon_dg2d(dg2dfile,stratum=stratum,wallid=1,debug=debug,commonzone=True)

    wallvertices_ordered = []
    wallvertices_ordered.append(lowestRnode)

    for i in range(0,len(wallvertices)-1):
        if i==0:
            prevnode = lowestRnode
        else:
            prevnode = wallvertices_ordered[i-1]

#        next_node = find_next_wall_node(wallvertices_ordered[i],prevnode,wallvertices,wallpolys,i==0)
#        if i == 85:
        if False:
            rtemp = []
            ztemp = []
            for node in wallvertices_ordered:
                rtemp.append(node.coords[0])
                ztemp.append(node.coords[1])
            plt.plot(rtemp,ztemp,"-")
            plt.plot(rtemp,ztemp,"-")
            rtemp = []
            ztemp = []
            for node in vertices:
                rtemp.append(node.coords[0])
                ztemp.append(node.coords[1])
            plt.plot(rtemp,ztemp,".")
            plt.xlim([1.3,1.5])
            plt.ylim([-1.2,-0.5])
            plt.savefig("test.png")
            plt.close()
        next_node = find_next_wall_node(wallvertices_ordered[i],prevnode,vertices,polys,i==0)
        wallvertices_ordered.append(next_node)

    # wallvertices_ordered is now the list of wall vertices, starting at 
    # the minimum R value, and going around counter-clockwise
#    print("wallvertices:")
#    for node in wallvertices:
#        print(node.coords,node.id)
#    print("wallvertices_ordered:")
#    for node in wallvertices_ordered:
#        print(node.coords,node.id)

    nwall = len(wallvertices_ordered)
    segments=np.array(range(0,nwall),dtype=int)
    rwall = np.zeros(nwall)
    zwall = np.zeros(nwall)
    source_strength = []
    for iseg in range(0,nwall):
        vertex1 = wallvertices_ordered[iseg]
        point1 = np.array( vertex1.coords )
        iseg2 = (iseg+1)%nwall
        vertex2 = wallvertices_ordered[iseg2]
        rwall[iseg] = vertex1.coords[0]
        zwall[iseg] = vertex1.coords[1]
        point2 = np.array( vertex2.coords )
        diff = point2 - point1
        area = 2.0*np.pi*0.5*(point1[0]+point2[0])*np.linalg.norm(diff)
        normal = [-diff[1],0.0,diff[0]]
        a_unit = normal/np.linalg.norm(normal)

        Br = 0.5*(Br_node[vertex1.id] + Br_node[vertex2.id])
        Bt = 0.5*(Bt_node[vertex1.id] + Bt_node[vertex2.id])
        Bz = 0.5*(Bz_node[vertex1.id] + Bz_node[vertex2.id])
        b_unit = [Br,Bt,Bz]/np.linalg.norm([Br,Bt,Bz])

        wetarea_fraction = np.abs(np.dot(b_unit,a_unit))

        vpar = 0.5*(np.sqrt(Te_node[vertex1.id]/ionmass)*mach_node[vertex1.id] + \
                np.sqrt(Te_node[vertex2.id]*mach_node[vertex2.id]/ionmass))

        ne = 0.5*(ne_node[[vertex1.id]] + ne_node[[vertex2.id]])

        source_strength.append(ne*vpar*wetarea_fraction)

    source_strength = np.array(source_strength).flatten()

    plt.plot(rwall,zwall)
    plt.tight_layout()
    plt.savefig("wallpoly.pdf")
    plt.close()

    # Enclose in universal cell

    Polygon.close_in_universal_cell(dg2dfile,wallvertices_ordered,0,stratum+1,material,recyc,clockwise=False,walltemp=300.0)

    dg2dfile.write("polygon_nc_file polygon.nc\n")
    dg2dfile.write("end")
    dg2dfile.close()

    ####################################################
    # Write wallfile
    wallfile = open("wallfile.txt","w")
    wallfile.write("# Wallfile automatically generated by script\n")
    wallfile.write("1\n")
    wallfile.write("#\n")
    wallfile.write("%d\n"%(nnode))
    wallfile.write("#\n")
    for i in range(0,nnode):
        x = node_coords[i,0]
        z = node_coords[i,1]
        wallfile.write("%f %f\n"%(x,z))
    wallfile.close()

    segments=np.array(range(0,len(source_strength)),dtype=int)

    strata = np.array([stratum+2]*len(segments))

    return ne_zone, Te_zone/1.602e-19, Ti_zone/1.602e-19, strata, segments, source_strength, zone_map

def write_dg2d_input_from_single_wall(wallfile_name,material,recyc,walltemp=300.0,minarea=-1.0,dg2dfile_name="dg2d.in",
                                      polygon_filename="polygons.nc",debug=False,def_separatrix=False,clockwise=False,exitnodes=[]):
    """ Function to write the 'dg2d.in' file from 'wallfile.txt' containing a single wall.

    Args:
        wallfile_name: string filename for the wallfile. (usually='wallfile.txt')
        mat: string wall material name.
        recyc: (float) recycling coefficient.
        walltemp: (float, optional) temperature of wall in Kelvin.
        minarea: (float, optional) ???
        dg2dfile_name: (str, optional) 
        polygon_filename: (str, optional)
        debug: (bool, optional)
        def_separatrix: (bool, optional)
        clockwise: (bool, optional)
        exitnodes: (list, optional)
    Returns:
        polys: 
        wall:
        if def_separatrix:
            Rcoords_sep:
            Zcoords_sep:
    """
    # Read the wallfile
    wallfile = open(wallfile_name,'r')
    
    def next_noncomment_line(f):
        """ local helper function to incrementally read lines from a file.
        In DEGAS2, lines beginning with whitespace or '#' are comments.
        """
        found=False
        while not found:
            line = f.readline().strip()
            if not line[0] == '#':
                found = True
        return line

    # As per the definegeometry2d docs,
    #   the first line is the number of walls.
    nwalls = int(next_noncomment_line(wallfile))
    #   the next line(s) indicate the number of points comprising each wall.
    #   e.g. 100 200 100 # for three walls.
    line = next_noncomment_line(wallfile).split()
    # We only consider the first wall here,
    npoints = int(line[0]) 

    # If a second wall is given, it's taken to define the separatrix, 
    if def_separatrix:
        npoints_sep = int(line[1])

    nlines = npoints
    iwall = 0
    ipoint = 0
    Rcoords = []
    Zcoords = []
    wall = Surface(0)
    for ipoint in range(0,npoints):
        line = next_noncomment_line(wallfile)
        line_floats = [float(n) for n in line.split()]
        Rcoord = line_floats[0]
        Zcoord = line_floats[1]
        Rcoords.append(Rcoord)
        Zcoords.append(Zcoord)
        wall.add_vertex(Vertex(ipoint,Rcoord,Zcoord))

    if def_separatrix:
        sep = Surface(1)
        Rcoords_sep = []
        Zcoords_sep = []
        for ipoint in range(0,npoints_sep):
            line = next_noncomment_line(wallfile)
            line_floats = [float(n) for n in line.split()]
            Rcoord = line_floats[0]
            Zcoord = line_floats[1]
            Rcoords_sep.append(Rcoord)
            Zcoords_sep.append(Zcoord)
     
    # Now we have read the wallfile. We know how many points are in each wall and the R/Z coordinates for each point defining these walls. We can close the wallfile now.
    wallfile.close()

    polys=[]

    # There will be only three polygons:
    # 1: the vacuum vessel interior
    # 2: Small polygon representing exit on LFS, connects to right side of universal cell
    # 3: Most of the wall to the edge of the universal cell. Shares zone with #2. 

    # Polygon 1
    poly = Polygon()
    poly.add_whole_surface(wall,backward=True)
    
    # Now we build the wall specification in the dg2d input file
    dg2dfile = open(dg2dfile_name,'w')

    ####################################
    # First, set the header
    dg2dfile.write("symmetry cylindrical\n")

    ####################################
    # Find and write appropriate bounds for the universal cell,
    # Xmin Xmax Zmin Zmax (for 2d cases)
    Rrange = np.amax(Rcoords) - np.amin(Rcoords)
    Zrange = np.amax(Zcoords) - np.amin(Zcoords)
    Zmin = np.amin(Zcoords) - Zrange
    Zmax = np.amax(Zcoords) + Zrange
    Rmin = max(0.0001,0.5*np.amin(Rcoords))
    Rmax = np.amax(Rcoords) + Rrange
    dg2dfile.write("bounds     %f %f    %f %f \n" % (Rmin, Rmax, Zmin, Zmax))

    ####################################
    # Write out name of wallfile
    dg2dfile.write("wallfile "+wallfile_name.strip()+"\n")

    dg2dfile.write("end_prep\n")
    dg2dfile.write("\n")
    # End of dg2d.in preparatory section.
    ####################################
    # Start of dg2d.in construction section.
    
    poly.write_plasma_polygon_dg2d(dg2dfile,minarea=minarea,wallid=1,debug=debug)

    if exitnodes == []:
        #wall.write_solid_polygon_dg2d(2,dg2dfile,material,recyc,walltemp,debug=debug)
        Polygon.close_in_universal_cell(dg2dfile,wall.vertices,0,2,material,recyc,walltemp=walltemp,debug=False,clockwise=clockwise)
    else:
        wall.write_solid_polygon_with_exit(2,dg2dfile,exitnodes,material,recyc,walltemp,debug=debug)

    dg2dfile.write("polygon_nc_file "+polygon_filename+"\n")

    if debug:
        dg2dfile.write("quit\n")

    dg2dfile.write("end\n")

    dg2dfile.close()

    Polygon.clear_numPolygon()

    if def_separatrix:
        return polys,wall, Rcoords_sep, Zcoords_sep
    else:
        return polys,wall

def refine_limiter(r_in,z_in,maxdist):
    """
    Method that refines a polygon to have more points. Useful for increasing resolution at the wall. Currently the only way to refine a triangulation.

    Args:
        r: List of floats for the R coordinate of the polygon to refine.
        z: List of floats for the Z coordinate of the polygon to refine.
        maxdist: Float; maximum allowed distance between vertices. Method halves distance until this condition is specified.
    """

    r = np.copy(r_in)
    z = np.copy(z_in)
    if not (np.abs(r[0]-r[-1]) < 1.0e-4 and np.abs(z[0] - z[-1]) < 1.0e-4):
        r = np.append(r,r[0])
        z = np.append(z,z[0])

    N = len(r)
    rnew = []
    znew = []
    for i in range(0,N):
        if i == 0:
            rnew.append(r[i])
            znew.append(z[i])
        else:
            # Original r[i],z[i] are candidates for next wall node
            # 
            dist = np.sqrt( (r[i]-rnew[-1])**2 + (z[i]-znew[-1])**2)
            if dist <= maxdist:
                rnew.append(r[i])
                znew.append(z[i])
            else:
                ndiv = 0 
                while dist > maxdist:
                    dist = 0.5*dist
                    ndiv += 1
                # Had to divide the segment ndiv times. So each segment has a length dist/2**ndiv
                # and we have to insert (2**ndiv-1) points before moving on
                unit = np.array([r[i]-rnew[-1],z[i]-znew[-1]])
                unit = unit/np.linalg.norm(unit)
                for j in range(0,2**ndiv-1):
                    rnew.append(rnew[-1]+dist*unit[0])
                    znew.append(znew[-1]+dist*unit[1])
    return rnew, znew

def refine_limiter_interp(r, z, maxdist):
    """ Alternative method using interpolation to refine the limiter.
    Args:
        r: List of floats
        z: List of floats
        maxdist: Float, 
    Returns:
        rnew: List of floats
        znew: List of floats
    """
    r = np.array(r)
    z = np.array(z)
    d = np.cumsum(np.sqrt(np.diff(r)**2 + np.diff(z)**2))
    d = np.concatenate(([0],d)) # parameterize limiter curve by distance
    # New target grid based on maxdist
    dl = 0.9*maxdist
    di = np.arange(0, max(d) + dl, dl)
    rnew = np.interp(di, d, r, period=max(d))
    znew = np.interp(di, d, z, period=max(d))

    return rnew.tolist(), znew.tolist()

def generateGeometryFromEFITfile(efitfile,mat,recyc_coef=1.0,Twall=300.0,wfilename="wallfile.txt",
                                 clockwise=False,RZlim_start=None,def_separatrix=False,
                                 dlim_max=0.05,dsep_max=0.05,minarea=-1.0):
    """
    Writes dg2d.in from an EFIT g file using definegeometry2d's built-in triangulation.

    Args:
        efitfile: string filename for the EFIT g file.
        mat: string wall material name.
        recyc_coeff: (optional) float for recycling coefficient.
        Twall: (optional) float for temperature of wall in Kelvin.
        wfilename: (optional) string for the name of the "wallfile". Included for compatibility with legacy scripts.
        clockwise: (optional) boolean flag for whether the limiter points in g file are clockwise or not. Default False.
        RZlim_start: (optional) [R, Z] point at which to start writing the limiter.
                     The closest point on the limiter is used to roll the rlim, zlim arrays after refinement.
        def_separatrix: (bool, optional) pased to write_dg2d_input_from_single_wall
        dlim_max: (float, optional) float for maximum distance between limiter points in meters. Used to refine triangulation near wall.
        dsep_max: (float, optional) float for maximum distance between separatrix points in meters. Used to refine triangulation near separatrix.
        minarea: (float, optional) passed to write_dg2d_input_from_single_wall
    Returns:
        psi_func: function that takes R,Z  as arguments and returns normalized psi coordinate
        nodes: list of Vertices representing the wall. Used to pass to other routines.
    """
    # ----------------
    # read_geqdsk() returns a class with attrs read from the EFIT gEQDSK file. 
    g = geomutils.read_geqdsk(efitfile)
    rlim = g.lim[:,0]
    zlim = g.lim[:,1]
    print(f"gEQDSK limiter has {len(rlim)} points")
    rgrid = g.rgrid
    zgrid = g.zgrid
    # Defines psi_n(R,Z) = (psi - min(psi))/(max(psi) - min(psi)), Normalized poloidal flux.
    psi_rz = (g.psirz-g.ssimag)/(g.ssibry-g.ssimag)
    
    # R = list of R points, Z list of Z point
    def psi_func(R,Z):
        # interp2d is deprecated, use RectBivariateSpline, 
        f = interpolate.RectBivariateSpline(rgrid, zgrid, psi_rz.T)
        return f(R,Z)[0]

    # ----------------
    # refine/reorient the limiter,
    rlim, zlim = refine_limiter_interp(rlim,zlim,dlim_max)
    N_lim = len(rlim)
    print(f"refined limiter has {N_lim} points")
    
    if RZlim_start is not None:
        print(f"* finding point closest to: {RZlim_start}")
        rlim = np.array(rlim)
        zlim = np.array(zlim)
        dist = np.sqrt((rlim - RZlim_start[0])**2 + (zlim - RZlim_start[1])**2)
        i_start = np.argmin(dist) 
        # roll such that i_start = 0
        rlim = np.roll(rlim, -i_start)
        zlim = np.roll(zlim, -i_start)
        rlim = rlim.tolist()
        zlim = zlim.tolist()
        print(f"* limiter starts at: [{rlim[0]:.3f},{zlim[0]:.3f}]")
            
    # ----------------
    # Open and begin writing the wallfile,
    wallfile = open(wfilename,"w")
    wallfile.write("# Wallfile automatically generated by script\n")
    
    N_walls = 1 if not def_separatrix else 2
    wallfile.write(f"{int(N_walls)}\n") # First line= number of walls in the file (=1)
    wallfile.write("#\n")
    # For each limiter point,
    #   if the i'th point is a unique point (not overlapping with prev or next point) OR i==0,
    #       add it to a list of 'nodes' containing Vertex objects.
    wall_nodes = []
    for i in range(0,N_lim):
        if i == 0:
            # Always write the 0th point,
            wall_nodes.append(Vertex(i, rlim[i], zlim[i]))
        else:
            # Store only unique limiter points,
            if ((rlim[i]-rlim[i-1])**2 + (zlim[i]-zlim[i-1])**2 > 1.0e-6) and ((rlim[i]-rlim[0])**2 + (zlim[i]-zlim[0])**2 > 1.0e-6):
                wall_nodes.append(Vertex(i,rlim[i],zlim[i]))
    N_wall_nodes = len(wall_nodes)
    print(f"found {N_wall_nodes} valid wall nodes")
    N_total_points = [N_wall_nodes]

    if def_separatrix:
        # Get the separatrix directly from gEQDSK file,
        rsep = g.sep[:,0]
        zsep = g.sep[:,1]
        N_sep = len(rsep)
        print(f"gEQDSK separatrix has {N_sep} points")
        rsep, zsep = refine_limiter_interp(rsep, zsep, dsep_max)
        N_sep = len(rsep)
        print(f"Refined separatrix has {N_sep} points")
        sep_nodes = []
        for i in range(0,N_sep):
            if i == 0:
                # Always write the 0th point,
                sep_nodes.append(Vertex(i, rsep[i], zsep[i]))
            else:
                # Store only unique separatrix points,
                if ((rsep[i]-rsep[i-1])**2 + (zsep[i]-zsep[i-1])**2 > 1.0e-6) and ((rsep[i]-rsep[0])**2 + (zsep[i]-zsep[0])**2 > 1.0e-6):
                    sep_nodes.append(Vertex(i,rsep[i],zsep[i]))
        # ---
        N_sep_nodes = len(sep_nodes)
        print(f"found {N_sep_nodes} valid separatrix nodes")
        N_total_points += [N_sep_nodes]
        
    # Next line(s)= number of points comprising each wall.
    line = " ".join([f"{int(n)}" for n in N_total_points])
    wallfile.write(f"{line}\n") 
    wallfile.write("#\n")
    count = 0
    for i in range(0,N_wall_nodes):
        wallfile.write("%f %f\n"%(wall_nodes[i].coords[0],wall_nodes[i].coords[1]))
        count += 1
    if def_separatrix:
        for i in range(0,N_sep_nodes):
            wallfile.write("%f %f\n"%(sep_nodes[i].coords[0],sep_nodes[i].coords[1]))
            count += 1
            
    wallfile.close()
    print(f"wrote {count} floats to wallfile")

    polys, wall = write_dg2d_input_from_single_wall(wfilename,mat,recyc_coef,walltemp=Twall,minarea=minarea,
                                                    dg2dfile_name="dg2d.in",polygon_filename="polygons.nc",
                                                    debug=False,exitnodes=[],def_separatrix=def_separatrix,
                                                    clockwise=clockwise)

    return psi_func, wall_nodes

def generateSimpleCylinderGeometry(rgrid,material,Twall,recyc,Lz=None,wallfile_name="wallfile.txt",dg2dfile_name="dg2d.in"):
    if not Lz:
        Lz = 5.0*rgrid[-1]
    Xmin = 0.5*rgrid[0]
    Xmax = 1.5*rgrid[-1]
    Zmin = -Lz
    Zmax = Lz

    ####################################################
    # Write wallfil
    wallnodes = []
    wallfile = open(wallfile_name,"w")
    wallfile.write("# Wallfile automatically generated by script\n")
    wallfile.write("#\n")
    wallfile.write("1\n")
    wallfile.write("#\n")
    wallfile.write("%d\n"%(2*len(rgrid)))
    wallfile.write("# \n")
    wallfile.write("%f %f\n"%(rgrid[0],-0.5*Lz))
    wallnodes.append(Vertex(0,rgrid[0],-0.5*Lz))
    for ir in range(0,len(rgrid)):
        wallfile.write("%f %f\n"%(rgrid[ir],0.5*Lz))
        wallnodes.append(Vertex(ir+1,rgrid[ir],0.5*Lz))
    for ir in range(len(rgrid)-1,0,-1):
        wallfile.write("%f %f\n"%(rgrid[ir],-0.5*Lz))
        wallnodes.append(Vertex(len(rgrid)+(len(rgrid)-ir),rgrid[ir],-0.5*Lz))
    wallfile.close()

    ####################################################
    # Build polygons from wallfile nodes
    polys = []
    for ir in range(0,len(rgrid)-1):
        polys.append(Polygon())

        polys[-1].add_vertex(wallnodes[ir+1])
        polys[-1].add_vertex(wallnodes[ir+2])
        polys[-1].add_vertex(wallnodes[2*len(rgrid)-ir-1])
        polys[-1].add_vertex(wallnodes[(2*len(rgrid)-ir)%(2*len(rgrid))])
     
    ####################################################
    # Write dg2d input file
    dg2dfile = open(dg2dfile_name,"w")

    write_dg2d_header(dg2dfile,"cylindrical",Xmin,Xmax,Zmin,Zmax,wallfile_name=wallfile_name)

    for ipoly in range(0,len(polys)):
        polys[ipoly].write_plasma_polygon_dg2d(dg2dfile,stratum=ipoly+1,wallid=1,commonzone=True)

    dg2dfile.write("new_zone solid\n")
    dg2dfile.write("new_polygon\n")
    dg2dfile.write("  material "+material+"\n")
    dg2dfile.write("  temperature "+str(Twall)+"\n")
    dg2dfile.write("  recyc_coef "+str(recyc)+"\n")
    dg2dfile.write("  stratum "+str(len(rgrid))+"\n")
    dg2dfile.write("  wall 1 %d %d \n"%(len(rgrid)+1,len(rgrid)+1))
    dg2dfile.write("  wall 1 %d %d \n"%(len(rgrid),len(rgrid)))
    dg2dfile.write("  outer 2 3\n")
    dg2dfile.write("  triangulate_polygon\n\n")

    dg2dfile.write("new_zone solid\n")
    dg2dfile.write("new_polygon\n")
    dg2dfile.write("  material mirror\n")
    dg2dfile.write("  recyc_coef 1.0\n")
    dg2dfile.write("  stratum "+str(len(rgrid)+1)+"\n")
    dg2dfile.write("  outer 0 1 2\n")
    dg2dfile.write("  wall 1 %d %d\n"%(len(rgrid),0))
    dg2dfile.write("  wall 1 %d * reverse\n"%(len(rgrid)+1))
    dg2dfile.write("  outer 3 4\n")
#    dg2dfile.write("  triangulate_polygon\n\n")
    dg2dfile.write("  breakup_polygon\n\n")

    dg2dfile.write("polygon_nc_file polygons.nc\n")
    dg2dfile.write("end")
    # Syntax is a little confusing right now. "wallid" in this routine actually assumes 
    # they're counted from zero and corrects to write dg2d input file
    

    dg2dfile.close()

def write_wallfile(nodes,wallfile_name="wallfile.txt"):
    Nnode = len(nodes)
    wallfile = open(wallfile_name,"w")
    wallfile.write("# Wallfile automatically generated by script\n")
    wallfile.write("1 \n")
    wallfile.write("#\n")
    wallfile.write("%d\n"%(Nnode))
    wallfile.write("#\n")
    for i in range(0,Nnode):
        wallfile.write("%f   %f \n"%(nodes[i].coords[0],nodes[i].coords[1]))
    wallfile.close()
 
def get_triangulation(Nzone,polygonfile="polygon.nc"):
    p = nc.Dataset(polygonfile,"r")
    poly_zone = p["g2_polygon_zone"][:]
    poly_rz = p["g2_polygon_xz"][:,0:3,:]
    p.close()

    conn = -np.ones([Nzone,3],dtype=int)
    r = []
    z = []
    k=0
    eps = 1.0e-5

    for i in range(0,Nzone):
        izone = poly_zone[i]-1
        for j in range(0,3):
            idx_r = np.argwhere( np.abs(r[:]-poly_rz[i,j,0])<eps)
            if idx_r.size == 0:
                idx_z = np.array([])
            else:
                idx_z = np.argwhere( np.abs(z[idx_r]-poly_rz[i,j,1])<eps)

            if idx_z.size == 0:
                np.append(r, poly_rz[i,j,0])
                np.append(z, poly_rz[i,j,1])
                conn[izone,j] = len(r)-1
            elif idx_z.size != 1:
                print("ERROR: found more than one point already stored")
            else:
                conn[izone,j] = idx_r[idx_z[0]] 
        
    if np.any(conn == -1):
        print("ERROR: did not fill in conn array in get_triangulation")

    rz = np.zeros([len(r),2])
    rz[:,0] = r
    rz[:,1] = z
        
    return rz, conn
