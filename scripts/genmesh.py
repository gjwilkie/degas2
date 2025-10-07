import matplotlib.tri as mtri
import numpy as np
from tqdm import tqdm
import sys
from multiprocessing import Pool
from tqdm.contrib.concurrent import process_map

# Gets triangulation object from triangle files .node and .ele. Start here.
def get_triangulation(trifile_base):
    f=open(trifile_base+".node","r")
    nnode=int(f.readline().split()[0])
    rz=np.zeros([nnode,2])
    for i in range(0,nnode):
        line=f.readline().split()
        rz[i,0]=float(line[1])
        rz[i,1]=float(line[2])
    f.close()

    f=open(trifile_base+".ele","r")
    ntri=int(f.readline().split()[0])
    conn=np.zeros([ntri,3],dtype=int)
    for i in range(0,ntri):
        line=f.readline().split()
        conn[i,0]=int(line[1])-1
        conn[i,1]=int(line[2])-1
        conn[i,2]=int(line[3])-1
    f.close()

    triang = mtri.Triangulation(rz[:,0],rz[:,1],conn)
    return triang, rz, conn

def extend_mesh_into_wall(rz_in,conn_in,edgemap_in):
    # Create small quadrlitaerals just outside wall for each outer segment. 

    # Find the edges that have only one adjacent triangle

    edgemap_in[:,:]

    if not self.is_clockwise():
       thickness = -thickness 
        # Go around plasma polygon and accumulate new vertices 
        # offset outward by thickness. "Outward" is defined in a vector sense:
        # find the bisecting vector between adjacent segments and negate it.
    for i in range(0,Nvertex):

        # v0, v1, and v2 are in order around the wall polygon.
        # Parity is accounted for by the sign on thickness.
        v1 = self.vertices[i]
        if self.is_clockwise():
            v0 = self.vertices[np.mod(i-1,Nvertex)]
            v2 = self.vertices[np.mod(i+1,Nvertex)]
        else:
            v2 = self.vertices[np.mod(i-1,Nvertex)]
            v0 = self.vertices[np.mod(i+1,Nvertex)]

        # diff1 and diff2 are unit vectors pointing away form v1
        diff1 = np.array(v0.coords) - np.array(v1.coords)
        diff2 = np.array(v2.coords) - np.array(v1.coords)
        diff1 = diff1/np.linalg.norm(diff1)
        diff2 = diff2/np.linalg.norm(diff2)

        # Sum of two unit vectors bisects the angle between them
        bisect = diff1 + diff2

        # Now, determine if the bisection points "inward" or "outward"
        # For clockwise (default), "inward" is "to the right" of diff2.
        convex = np.sign(np.cross(diff1,diff2))
        

        if np.linalg.norm(bisect) < 1.0e-6:
            # For colinear or nearly colinear points,
            # the outward unit vector is the cross product between diff1
            # and phi=(r cross z). 
            outward = [diff1[1],-diff1[0]]
        else:
            # Otherwise, outward is the negation of the bisecting vector.
            outward = -convex*bisect/np.linalg.norm(bisect)

            newvertices.append(Vertex(firstid+i,v1.coords[0]+thickness*outward[0],v1.coords[1]+thickness*outward[1]))
            outpoly.add_vertex(newvertices[-1])

   



# Generate an equation for a cone for every pair of points x1 and x2 (size [Nsurfs,2]).
# Returns [Nsurfs,10]  coefficients
def gen_cones(x0,x1):
   Nsurfs = np.shape(x1)[0] 
   coeffs = np.zeros([Nsurfs,10])

   eps = 1.0e-8
   eps_angle = np.sqrt(2.0*1.0e-10)

   cylcond = (np.abs(x1[:,0]-x0[:,0]) < eps_angle*np.abs(x1[:,1]-x0[:,1])) 
   cylcond = np.logical_or(cylcond, (np.abs(x1[:,0]-x0[:,0]) < eps))

   planecond = (np.abs(x1[:,1]-x0[:,1]) < eps_angle*np.abs(x1[:,0]-x0[:,0])) 
   planecond = np.logical_or(planecond, (np.abs(x1[:,1]-x0[:,1]) < eps))

   m = np.divide((x1[:,1]-x0[:,1]),(x1[:,0]-x0[:,0]),out=np.ones_like(x1[:,0]),where=np.logical_not(cylcond))
   # m=1 for cylinder, slope for a cone, and 0 for a plane
   m2 = np.where(planecond, np.zeros_like(m), m*m)
   m2 = np.where(cylcond, np.ones_like(m), m2)
   # b = 0 for a cylinder, 1/2 for a plane, and intercept for a cone
   b = np.where(cylcond, np.zeros_like(m), x1[:,1] - m*x1[:,0])
   b = np.where(planecond, 0.5*np.ones_like(m), b)

   # c0 is -b^2 for a cone, -R^2 for a cylinder, -Z0 for a plane
   c0 = np.where( cylcond , -x0[:,0]**2, -b**2)
   c0 = np.where( planecond, -x0[:,1], c0)

   coeffs[:,0] = c0
   coeffs[:,3] = 2.0*b
   coeffs[:,4] = m2
   coeffs[:,5] = m2
   coeffs[:,6] = np.where(np.logical_or(planecond,cylcond),np.zeros(Nsurfs),-np.ones(Nsurfs))

   return coeffs

def orient_surfaces(coeffs, points):
    eps = 1.0e-16

    # Planes have no coefficients on x^2 or y^2
    planecond = np.abs(coeffs[:,4]) < eps

    # Only cones have a coefficient on z^2, which will be -1
    conecond = np.abs(coeffs[:,6]+1.0) < eps

    # Cylinders are neither planes nor cones
    cylcond = np.logical_not(np.logical_or(planecond,conecond))

    b = 0.5*coeffs[:,3]

    # How to determine sign of surface relevant for each point:
    #  - Find a point on the surface (cone_r) at the same Z as the point
    #  - If R of the input point is > the surface point's R, it's outside the surface.
    #  - If surface is a plane, just compare Z(point) to Z(surface)

    # If surface is a cone, use R at common Z as the coordinate to compare
    #  otherwise (eventually, plane only), use Z0
    m = np.where(conecond, np.sqrt(coeffs[:,4]), np.ones_like(b) )
    surfpoint = np.where(conecond,  np.abs(points[:,1] - b)/m,-coeffs[:,0])

    # If surface is a cylinder, use R of surface to compare
    c0 = np.where(cylcond,-coeffs[:,0],np.ones_like(b))
    surfpoint = np.where(cylcond,np.sqrt(c0), surfpoint)

    # What is this for a cylinder: b=0, m=1, so cone_r is Z(point)
    # Should instead use c0 to get either b^2 or R^2

    if np.any( np.logical_or(np.isnan(surfpoint), np.isinf(surfpoint))):
        print("ERROR: surface points to be compared is populated with NaN, which should have been ruled out by plane condition.")

    sign_cond = np.where(planecond,points[:,1]>surfpoint,points[:,0] > surfpoint)

    return np.where(sign_cond, 1, -1)

# First try: loop over triangles. Not as performant
#    def populate_edge_map(itri):
#        for iedge in range(0,3):
#            idx1 = iedge
#            idx2 = (iedge+1)%3
#            # Every pair in edges array has indices in decreasing order
#            edgemap[itri,iedge] = edgelist.index([max(conn[itri,idx1],conn[itri,idx2]),min(conn[itri,idx1],conn[itri,idx2])])
#
#    pool = Pool()
#    pool.imap(populate_edge_map,tqdm(range(0,Ntri)))
#    for itri in range(0,Ntri):
#        populate_edge_map(itri)
#    for itri in tqdm(range(0,Ntri)):
#        populate_edge_map(itri)
#    for itri in tqdm(pool.imap(populate_edge_map,range(0,Ntri))):
#        populate_edge_map(itri)

# This way (looping over edges and assigning them to triangles rather than vice versa)
# is more efficient even though there are more edges. We can make better use of numpy this way.

def populate_edge_map(iedge):
    global edgemap
    candidates = np.argwhere(conn == edges[iedge,0])[:,0]
    new_candidates = np.argwhere(conn[candidates,:] == edges[iedge,1])[:,0]
#     print(iedge, candidates, new_candidates)
    tri_idxs = candidates[new_candidates]
    if len(tri_idxs) > 2:
        print("ERROR: more than 2 triangles adjacent to an edge; doesn't make sense")
    if len(tri_idxs) == 0:
        print("ERROR: could not find a triangle with this edge")
    for itri in tri_idxs:
        if (edgemap[itri,0] < 0):
            edgemap[itri,0] = iedge
        elif (edgemap[itri,1] < 0):
            edgemap[itri,1] = iedge
        elif (edgemap[itri,2] < 0):
            edgemap[itri,2] = iedge

if __name__=="__main__":
    # Generate a geometry.nc file for a large mesh
    tribase = str(sys.argv[1])

    triang,rz,conn=get_triangulation(tribase)

    # For every pair of connected points, there's only one edge. A good proxy for surfaces.
    edges = triang.edges
    edgelist = edges.tolist()
    Nedge = len(edges[:,0])
    Ntri = len(conn[:,0])
    edgemap = -1*np.ones([Ntri,3],dtype=int)

    rzedges1 = rz[edges[:,0],:]
    rzedges2 = rz[edges[:,1],:]
    coeffs = gen_cones(rzedges1,rzedges2)
    # We now have surface_coeffs array!

    # Now, map to triangles
    # Brute force: search for point pair in edges array
#    edgemap = np.zeros([Ntri,3],dtype=int)


# Uncomment is not using saved edgemap:
#    pool = Pool()
#    pool.imap(populate_edge_map,range(0,Nedge))
#    pool.imap(populate_edge_map,tqdm(range(0,Nedge)))

    def pairint(a,b):
        return (a+b)*(a+b+1)/2 + a

    print(Ntri,Nedge)
    print("Populating edgemap...")
    print("  sorting vertices...")
    # We can count on edges being sorted in decending node index order
    conn_sort = np.sort(conn,axis=1)
    edge1 = conn_sort[:,[1,0]]
    edge2 = conn_sort[:,[2,0]]
    edge3 = conn_sort[:,[2,1]]
    print("  matching pairs...")
    edges_c = pairint(edges[:,0],edges[:,1])
    edge1_c = pairint(edge1[:,0],edge1[:,1])
    edge2_c = pairint(edge2[:,0],edge2[:,1])
    edge3_c = pairint(edge3[:,0],edge3[:,1])
    for i in tqdm(range(0,Ntri)):
        edgemap[i,0] = np.argwhere(np.equal(edge1_c[i], edges_c))[0][0]
        edgemap[i,1] = np.argwhere(np.equal(edge2_c[i], edges_c))[0][0]
        edgemap[i,2] = np.argwhere(np.equal(edge3_c[i], edges_c))[0][0]
#    print(np.shape(edge1_c), np.shape(edges_c))
#    idx1 = np.isin(edge1_c,edges_c)
#    idx2 = np.isin(edge2_c,edges_c)
#    idx3 = np.isin(edge3_c,edges_c)
#    print(np.shape(idx1), np.shape(idx2),np.shape(idx3))
#    edgemap[:,0] = np.where(idx1)[0]
#    edgemap[:,1] = np.where(idx2)[0]
#    edgemap[:,2] = np.where(idx3)[0]


#    for iedge in tqdm(range(0,Nedge)):
#    for iedge in pool.imap(populate_edge_map,range(0,Nedge)):
#    for iedge in tqdm(range(0,Nedge)):
#        populate_edge_map(iedge)
#    process_map(populate_edge_map,range(0,Nedge),chunksize=100)
    np.save("edgemap.npy",edgemap)
#    edgemap = np.load("edgemap.npy")
#    edgemap = np.array(edgemap,dtype=int)

    coeffs_all = np.zeros([Nedge+4+2*Ntri,10])
    coeffs_all[4:Nedge+4,:] = coeffs

    # Create universal cell coefficients
    # Original order: rmin, zmin, rmax, zmax
    rmin = 0.5*np.min(rz[:,0])
    dr = (np.max(rz[:,0]) - np.min(rz[:,0]))
    rmax = 2*rmin+dr
    zmax = np.max(rz[:,1])+dr
    zmin = np.min(rz[:,1])-dr
    # Surface directions don't make sense for universal cell, but I'm not sure they need to
    coeffs_all[0,:] = [-rmin**2,0,0,0,1,1,0,0,0,0]
    coeffs_all[1,:] = [1,0,0,-(1.0/zmax),0,0,0,0,0,0]
    coeffs_all[2,:] = [1,0,0,0,-(1.0/rmax**2),-(1.0/rmax**2),0,0,0,0]
    coeffs_all[3,:] = [1,0,0,-(1.0/zmin),0,0,0,0,0,0]
    
    # Universal cell is first 4 surfaces
    # Next Nedge surfaces are the triangle edges
    # Next 2*NTri surfaces are cut surfaces

    # Every triangular cell gets two cut surfaces: planes at minimum and maximum z
    # Doesn't work yet, but is fast enough.
    boundaries = np.zeros(4+5*Ntri,dtype=int)
    cells = np.zeros([Ntri+1,4],dtype=int)
    cells[0,0:4] = [1,4,4,0]

    cells[1:Ntri+1,0] = 1+4*np.array(range(1,Ntri+1),dtype=int)
    cells[1:Ntri+1,1] = 3
    cells[1:Ntri+1,2] = 5
    cells[1:Ntri+1,3] = -1

    boundaries[0:4] = range(0,4)

    print("Building surface linking...")
    for i in tqdm(range(1,Ntri+1)):
        bdy_start = 4+(i-1)*5
        cut_start = 4+Nedge+2*(i-1)
        boundaries[bdy_start:bdy_start+3] = edgemap[i-1,:] + 1
        boundaries[bdy_start+3:bdy_start+5] = [cut_start+1,cut_start+2]

        zmin=np.min(rz[edges[edgemap[i-1,:]],1])
        zmax=np.max(rz[edges[edgemap[i-1,:]],1])

        coeffs_all[cut_start,:] = [zmin,0,0,1.0,0,0,0,0,0,0]
        coeffs_all[cut_start+1,:] = [zmax,0,0,-1.0,0,0,0,0,0,0]

# Nontrivial nc variables: 
# surfaces: four data for each surface counting faces and directions
# surface_sectors: four data for each surface counting sectors and directions
# sector_points: R,Z points on each side of sector segment. Needed for sources.
# sector_surface: surface for each sectoa (signed)
# sector_zone: zone for each sector
# sector_type_pointer:  points to one of vacuum, plasma, target, wall, exit, diagnostic (indexed accordingly). Manages different types of sectors as subclasses
# zone_index
# strata
# sector_strata_segment
# diagnostic_sector_tab
# diagnostic_delta
# diagnostic_num_sectors

# Example sector definition from the code:
#@m define_sector_plasma(sc)
#      sc_plasma_num=sc_plasma_num+1;
#      var_realloca(plasma_sector);
#      plasma_sector[sc_plasma_num]=sc;
#      sector_type_pointer[sc][sc_plasma]=sc_plasma_num;

# sector:strata::wallsegment:wall



#dimensions:
# 	vector = 3 ;
# 	string = 300 ;
# 	cell_info_ind = 4 ;
# 	cell_ind = 24977 ;          = Ntri+Nwall+1
# 	surface_ind = 74068 ;       = Nedge_tot
# 	boundary_ind = 111680 ;     = Nuniv(4?) + (Ntri+Nwall)*5
# 	neighbor_ind = 25569 ;      = (Ntri+Nwall)*3
# 	neg_pos = 2 ;
# 	surface_info_ind = 2 ;
# 	surface_tx_ind = 2 ;
# 	tx_ind_1 = 3 ;
# 	tx_ind_2 = 4 ;
# 	transform_ind = 1 ;
# 	coeff_ind = 10 ;
# 	zone_type_ind = 4 ;
# 	zone_index_ind = 4 ;
# 	zone_ind = 24131 ;
# 	sector_ind = 551 ;
# 	sector_neg_pos_ind = 2 ;
# 	sector_type_ind = 17 ;
# 	vacuum_ind = 1 ;
# 	plasma_ind = 276 ;
# 	target_ind = 276 ;
# 	wall_ind = 1 ;
# 	exit_ind = 1 ;
# 	sc_diag_name_string = 40 ;
# 	diag_grp_ind = 4 ;
# 	sc_diag_ind = 825 ;
# 	de_symbol_string = 24 ;
# 	de_name_string = 100 ;
# 	de_grp_ind = 1 ;
# 	de_zone_frags_ind = 10000 ;
# 	de_tot_view_ind = 1 ;
# 	de_start_end_ind = 2 ;
# 	de_view_ind = 1 ;
    




