import numpy as np
import copy
import sys
import matplotlib.pyplot as plt

# A polygon is an ordered set of vertices, connected by segments which close in on itself.
class Polygon:
    """
    A helper class to build polygons for definegeometry2d input. Not the most concise class it could be.

    Attributes:
        num_polygons: Global integer to keep track of the total number of polygons defined
        vertices: A list of points representing a closed polygon. Accessed through add_vertex method.
        id: Integer identification for the polygon. Usually corresponds to "stratum.
        alongwall: Boolean for whether this is an external (wall) point. Incomplete feature.
    """

    numPolygons = 0 

    def __init__(self,increment=True,id=None):
        """
        Constructor for an instance of Polygon
    
        Args:
            increment: If true (default), increments numPolygons.
            id: If specified, will set this to be the id of the polygon. Otherwise, will use numPolygons to increment. Inconsistent implementation.
        """
        self.vertices = []
        if id:
            self.id = id
        else:
            self.id = Polygon.numPolygons
        if increment:
            Polygon.numPolygons += 1
        self.alongwall = False

    # Polygons are equal if their vertices are equal
    def __eq__(self,other):
        return self.vertices == other.vertices

    def __ne__(self,other):
        return self.vertices != other.vertices

    def add_vertex(self,vertex):
        self.vertices.append(vertex)

    def clear_numPolygon():
        Polygon.numPolygons = 0

    def get_n_wallnodes(self):
        """
        Method to count the number of Vertices along Polygon where wall=True
        Returns:
            nwallnodes: number of wall nodes in this polygon
        """
        nwallnodes = 0
        for vertex in self.vertices:
            if vertex.wall:
                nwallnodes += 1
        return nwallnodes

    def get_first_wallnode(self):
        """
        Method to get the first vertex along the Polygon with the wall flag True.
        Returns:
            first: index of the first wall vertex
        """

        # Expects one contiguous set of nodes that are designated as wall nodes
        # Also expects at least one node that is not a wall node.
        first = -1
        count = 0
        nvertex = len(self.vertices)
        wallarray = [0]*nvertex 
        idx = 0
        for idx in range(0,nvertex):
            if self.vertices[idx].wall:
                count += 1
                if not self.vertices[idx-1].wall:
                    first = idx
        if first == -1 and (count < len(self.vertices)):
            sys.exit("Error: get_first_wallnode could not find the first wall node.")
        return first

    def reorder_wallnodes_first(self):
        """
        Method that reorders the nodes in the polygon so that those adjacent to the wall are listed first. Usually used for triangles.

        Returns:
            nwallnodes: integer for the number of wall nodes in this polygon
        """

        nwallnodes = self.get_n_wallnodes()
        nvertex = len(self.vertices)
        if nwallnodes >= 1:
            first = self.get_first_wallnode()
            temp = copy.deepcopy(self.vertices)
            for idx in range(0,nvertex):
                temp[idx] = self.vertices[(idx+first)%nvertex]
            self.vertices = copy.deepcopy(temp)

        # use return value if expecting more than one wall segment
        return nwallnodes 

    def get_next_vertex_extrapolated_to_wall(self,wall):
        """
        Used for niche cases.
        Linearly extrapolate the last two vertices defined for this polygon
        and find the point on the specified surface which is closest to said line.
        Returns index of said point
        Method that reorders the nodes in the polygon so that those adjacent to the wall are listed first. Usually used for triangles.

        Returns:
            closest_idx: index of the wall vertex 
        """


        # These last two points defined thus far for the polygon define
        # the line to which we are seeking the closest point on wall
        point1 = self.vertices[-2].coords
        point2 = self.vertices[-1].coords

        dir = [np.sign(point2[0]-point1[0]),np.sign(point2[1]-point1[1])]

        first = True
        idx = 0
        closest_idx = 0
        for vertex in wall.vertices:
            # Only check points that are on the correct half-plane according to 
            # the directionality of the last two points
            if ( (dir[0]*(vertex.coords[0] - point1[0])>0) or (dir[0] == 0)) and \
               ( (dir[1]*(vertex.coords[1] - point1[1])>0) or (dir[1] == 0)):
                current_distance_sq = ( (point2[1]-point1[1])*(point1[0]-vertex.coords[0] )
                    - (point2[0]-point1[0])*(point1[1]-vertex.coords[1] ) )**2 \
                    / ( (point2[1]-point1[1])**2 + (point2[0]-point1[0])**2 )
                    
                if first:
                    closest_vertex=vertex
                    closest_idx = idx
                    closest_distance_sq = current_distance_sq
                    first = False
                elif current_distance_sq < closest_distance_sq:
                        closest_vertex = vertex
                        closest_idx = idx
                        closest_distance_sq = current_distance_sq
            idx +=1

        return closest_idx

    # Adds surfaces to construct a polygon. Entire is an array to determine which surfaces
    # get added in their entirity. 
    # If entire[i] = True, the entire surface[i] gets added (e.g. an open or closed flux surface)
    # If entire[i] = False, only a subset that gets added. Which points to add are determined by 
    #                extrapolating from the previous and next entire surface
#    def construct_from_surfaces(self,surfaces,entire):
#        if (len(surfaces) == 2) and surfaces[0].closed and surfaces[1].closed:
#            # To deal with two closed surfaces, connect the first points of each to form polygon
#        else:
#            reverse = False
#            for isurf in range(0,len(surfaces)):
#                if entire[isurf]:
#                    for vertex in surfaces[isurf].vertices:
#                        self.add_vertex(vertex)
#                    # Reverse direction for next entire surface
#                    reverse = not reverse 
#                else:
#                    if isurf == 0:
#                        print("Error. Need to start with an entire surface in call to construct_from_surfaces. Don't know how to begin defining polygon.")
#                        sys.exit(1)
#                    idx = self.add_next_vertex_extrapolated_to_wall(surfaces[isurf])

    # wallnodes is a collection of integer identifiers that make up the outer wall
    # they must go *counter*-clockwise, start with the innermost, share a common wall (wallid)

    def close_in_universal_cell(f,wallnodes,wallid,stratum,material,recyc,debug=False,clockwise=False,walltemp=300.0):
        """
        A widely used method to envelop the plasma/vacuum domain within the universal cell; otherwise closed by a solid surface with specified properties

        Args:
            f: opened dg2d.in file ready to write the solid zones.
            wallnodes: list of Vertices that enclose the plasma zones.
            stratum: integer or string; first unused stratum to write
            material: string for the wall material to use
            recyc: recycling coefficient for the wall
            debug: (optional) If true, intends for definegeometry2d to write out the polygon for debugging instead of generating the geometry.
            clockwise: (optional) If true, wallnodes are specified in clockwise order.
            walltemp: (optional) float for the wall temperature in Kelvin
        """

        f.write("new_zone solid\n")
        f.write("new_polygon\n")
        f.write("  material "+material+"\n")
        f.write("  recyc_coef "+str(recyc)+"\n")
        f.write("  temperature "+str(walltemp)+"\n")
        f.write("  stratum "+str(stratum)+"\n")
        if clockwise:
            f.write("  wall "+str(wallid+1)+" "+\
                str(wallnodes[1].id)+" "+\
                str(wallnodes[1].id)+"\n")
            f.write("  wall "+str(wallid+1)+" "+\
                str(wallnodes[0].id)+" "+\
                str(wallnodes[0].id)+"\n")
        else:
            f.write("  wall "+str(wallid+1)+" "+\
                str(wallnodes[0].id)+" "+\
                str(wallnodes[0].id)+"\n")
            f.write("  wall "+str(wallid+1)+" "+\
                str(wallnodes[1].id)+" "+\
                str(wallnodes[1].id)+"\n")
        f.write("  outer 0 1\n")
        if debug:
            f.write("  print_polygon poly.out1.dat\n")
            f.write("  clear_polygon\n")
        else:
            f.write("  triangulate_polygon\n")
        f.write("\n")
        #f.write("new_zone solid\n")
        f.write("new_polygon\n")
        f.write("  material "+material+"\n")
        f.write("  recyc_coef "+str(recyc)+"\n")
        f.write("  temperature "+str(walltemp)+"\n")
        f.write("  stratum "+str(stratum+1)+"\n")
        written_ids = []
        if clockwise:
            f.write("  wall "+str(wallid+1)+" "+
                "0 0 "+"\n")
            written_ids.append(0)

            for node in wallnodes[-1:0:-1]:
                if not node.id in written_ids:
                    f.write("  wall "+str(wallid+1)+" "+\
                        str(node.id)+" "+\
                        str(node.id)+"\n")
                    written_ids.append(node.id)
        else:
            for node in wallnodes[1:]:
                if not node.id in written_ids:
                    f.write("  wall "+str(wallid+1)+" "+\
                        str(node.id)+" "+\
                        str(node.id)+"\n")
                    written_ids.append(node.id)
            f.write("  wall "+str(wallid+1)+" "+\
                str(wallnodes[0].id)+" "+\
                str(wallnodes[0].id)+"\n")
        f.write("  outer 1 2 3 4\n")
        if debug:
            f.write("  print_polygon poly.out2.dat\n")
            f.write("  clear_polygon\n")
        else:
            f.write("  triangulate_polygon\n")
        f.write("\n")


    # Writes the polygon to file f
    # If debug, include lines that output polygons to poly.X.dat files
    def write_plasma_polygon_dg2d(self,f,stratum=None,wallid=None,commonzone=False,minarea=-1.0,debug=False,newzone=True):
        """
        Writes a plasma Polygon to the dg2d.in file

        Args:
            f: opened dg2d.in file ready to write the plasma polygon
            stratum: (optional) Integer stratum to label this polygon. If not specified, uses the polygon's id plus one.
            wallid: (optional) Integer id of the wall label in wallfile. Usually 1 or the polygon id. If not specified, uses the polygon id.
            commonzone: (optional) Boolean flag. If true, entire polygon is a common zone. Otherwise, instructs definegeometry2d to break up polygons into triangles, with each a unique zone. Default False.
            minarea: (optional) Float for the minimum triangle area to pass to definegeometry2d. Does not appear to work. 
            newzone: (optional) Boolean flag. If false, dose not create a new zone but appends polygon to the existing zone. Default True.
            debug: (optional) Boolean flag for whether to instruct definegeometry2d to print out polygon instead of generating geometry for this zone. Default False.
        """
        if not wallid:
            wallnum = self.id+1 
        if not stratum:
            stratum = self.id+1

        if newzone:
            f.write("new_zone plasma\n")
        f.write("new_polygon\n")
        f.write("  stratum "+str(stratum)+"\n")
        # For each point write:
        #   wall {i_wall}, {j_start}, {j_end}
        # Takes points from the i'th wall in the wallfile and adds them to this polygon.
        # With j_start=j_end, a single point is added.
        for vertex in self.vertices:
            if not vertex.wall_id:
                f.write("  wall "+str(wallid)+" "+str(vertex.id)+" "+str(vertex.id)+"\n")
            else:
                f.write("  wall "+str(vertex.wall_id)+" "+str(vertex.id)+" "+str(vertex.id)+"\n")

        if debug:
            f.write("  print_polygon poly."+str(self.id)+".dat\n")
            f.write("  clear_polygon\n")
        elif commonzone:
            # triangulate_polygon is only recommended for solid zones.
            f.write("  triangulate_polygon\n")
        else: 
            f.write("  triangulate_to_zones\n")
        if minarea > 0.0:
            print(f"Setting minarea={minarea:.1e}")
            f.write("  triangle_area "+str(minarea)+"\n")

        f.write("\n")

    def plotpolygon(self,hold=False):
        import matplotlib.pyplot as plt
        R = []
        Z = []
        for vertex in self.vertices:
            R.append(vertex.coords[0])
            Z.append(vertex.coords[1])

        #plt.clf()
        plt.plot(R,Z,"+-")
        #plt.show()

    def plot_polygons(polys):
        import matplotlib.pyplot as plt

        plt.clf()
        for poly in polys:
            R = []
            Z = []
            for vertex in poly.vertices:
                R.append(vertex.coords[0])
                Z.append(vertex.coords[1])
            plt.plot(R,Z,"+-")
        plt.show()

    def add_wall_segment(self,wall,begin_idx,end_idx,backward=False):
        if backward:
            begin=begin_idx
            end=end_idx-1
            step=-1
        else:
            begin=begin_idx
            end=end_idx+1
            step=1
        for vertex in wall.vertices[begin:end:step]:
            self.add_vertex(vertex)
 
    def add_whole_surface(self,surf,backward=False):
        if backward:
            begin=-1
            end=None
            step=-1
        else:
            begin=0
            end=None
            step=1
        for vertex in surf.vertices[begin:end:step]:
            self.add_vertex(vertex)

    def is_clockwise(self,force=False):
        Nvertex = len(self.vertices)
        sumarea = 0.0
        result = True
        for i in range(0,Nvertex):
            r1 = self.vertices[i].coords[0]
            z1 = self.vertices[i].coords[1]
            r2 = self.vertices[np.mod(i+1,Nvertex)].coords[0]
            z2 = self.vertices[np.mod(i+1,Nvertex)].coords[1]
            sumarea += (r2-r1)*(z1+z2)

        if sumarea < 0.0:
            result = False
            if force:
                newpoly = Polygon(increment=False,id=self.id)
                for i in range(0,Nvertex):
                    newpoly.add_vertex(self.vertices[Nvertex-i-1])
                self = newpoly

        return result

            
    def build_aux_wall_polygons(self,firstid,thickness=0.005):
        """
        Uses a wall-defining polygon to build many individual zones adjacent to the segments that make up the limiter. 
        Allows more flexibility in defining the behavior of each individual segment. 
        Use with caution when limiter features are sub-10mm scale.
        
        Args:
            firstid: integer for the first id start building sequentially new vertices.
            thickness: (optional) float for the thickness of each auxiliary polygon in meters. Defaults to 5mm, but not universally consistent for small features.
        Returns:
            aux_polys: list of auxiliary Polygons representing independent wall behavior.
            outpoly: Polygon representing the outer points of the auxiliary polygons
            newvertices: Vertices added to create auxiliary wall polygons. To be appended to dg2d wallfile.
        """
    
        aux_polys = []
        newvertices = []
        outpoly = Polygon()

        self.is_clockwise(force=True)
        Nvertex = len(self.vertices)

        # Start with vertex closest to origin
        startidx = -1
        dist = 9.0e30
        for i in range(0,Nvertex):
            if np.norm(self.vertices[i].coords) < dist:
                dist = np.norm(self.vertices[i].coords)
                startidx = i

        v1 = self.vertices[startidx]
        v2 = self.vertices[np.mod(startidx+1,Nvertex)]
        normal = np.zeros(2)
        normal[0] = v1.coords[1]-v2.coords[1]
        normal[1] = v2.coords[0]-v1.coords[0]
        normal = normal / np.norm(normal)
        newvertices.append(Vertex(firstid,v1[0]+thickness*normal[0],v1[1]+thickness*normal[1]))
        outpoly.add_vertex(newvertices[-1])
            
        # Go around plasma polygon and accumulate new vertices offset outward by thickness
        for i in range(1,Nvertex):
            v1 = self.vertices[i]
            v2 = self.vertices[np.mod(i+1,Nvertex)]
            normal = np.zeros(2)
            normal[0] = v1.coords[1]-v2.coords[1]
            normal[1] = v2.coords[0]-v1.coords[0]
            normal = normal / np.norm(normal)

            newvertices.append(Vertex(firstid+i,v1[0]+thickness*normal[0],v1[1]+thickness*normal[1]))
            outpoly.add_vertex(newvertices[-1])

        # Build new polygons
        for i in range(0,Nvertex):
            newpoly = Polygon()
            newpoly.add_vertex(self.vertices[np.mod(i+1,Nvertex)])
            newpoly.add_vertex(self.vertices[i])
            newpoly.add_vertex(newvertices[i])
            newpoly.add_vertex(newvertices[np.mod(i+1,Nvertex)])
            aux_polys.append(newpoly)

        return aux_polys, outpoly, newvertices

# A surface is also an ordered set of vertices, but can be open or closed
# The vertices that make up various surfaces are combined to make a polygon
class Surface:
    numSurfaces = 0
    def __init__(self,id_in):
        self.vertices = []
        Surface.numSurfaces +=1
        self.id = Surface.numSurfaces
        self.closed = True

    def make_open(self):
        self.closed = False
    def make_closed(self):
        self.closed = True

    # Find the N vertices on surface that are closest to given vertex
    def get_N_closest_vertices(self,N,vertex_in):
        surf_use = copy.deepcopy(self)
        closest_idxs = []
        for i in range(0,N):
            first = True
            for iv in range(0,len(self.vertices)):
                distance_sq = (surf_use.vertices[iv].coords[0]-vertex_in.coords[0])**2 + \
                        (surf_use.vertices[iv].coords[1]-vertex_in.coords[1])**2 
                if first:
                    first=False
                    closest_distance_sq= distance_sq
                    closest_vertex_idx = iv
                    closest_vertex = self.vertices[iv]
                elif distance_sq < closest_distance_sq:
                    closest_distance_sq= distance_sq
                    closest_vertex_idx = iv
                    closest_vertex = self.vertices[iv]
            closest_idxs.append(closest_vertex_idx)
            surf_use.vertices[closest_vertex_idx].coords = [np.nan]*2
        return closest_idxs


    def add_vertex(self,vertex):
        self.vertices.append(vertex)

    # Take the given index, and add it the surface at the appropriate point
    # (as determined by the two closest vertices on the surface)
    def insert_vertex_into_surface(self,vertex_in):
        closest_idxs = self.get_N_closest_vertices(2,vertex_in)

        if not abs(closest_idxs[1] - closest_idxs[0]) == 1:
            print("Error: could not find two adjacent indices in which to insert new point.")
            sys.exit(1)

        self.vertices.insert(max(closest_idxs),vertex_in)

        return max(closest_idxs)

    def clear_numSurface():
        Surface.numSurfaces=0

    # Return two polygons defined by two closed surfaces.
    def genpolys_from_two_closed_surfs(inner,outer):
        poly1 = Polygon(increment=False)
        poly1.add_vertex(outer.vertices[0])
        poly1.add_vertex(outer.vertices[1])
        poly1.add_vertex(inner.vertices[1])
        poly1.add_vertex(inner.vertices[0])
        poly1.add_vertex(outer.vertices[0])

        poly2 = Polygon(increment=False)
        for vertex in inner.vertices[1:]:
            poly2.add_vertex(vertex)
        poly2.add_vertex(inner.vertices[0])
        poly2.add_vertex(outer.vertices[0])
        for vertex in outer.vertices[-1:0:-1]:
            poly2.add_vertex(vertex)
        poly2.add_vertex(inner.vertices[1])

        return poly1, poly2
            
    # Here we wish to include the whole surface and the wall between where it intersects
    def genpoly_from_rightwall_and_surface(wall,surface):

        poly = Polygon(increment=False)
        poly.add_whole_surface(surface,backward=True)

        wallidx_intersect_upper = poly.get_next_vertex_extrapolated_to_wall(wall)
        dummy = Polygon(increment=False)
        dummy.add_whole_surface(surface)
        wallidx_intersect_lower = dummy.get_next_vertex_extrapolated_to_wall(wall)

        poly.add_wall_segment(wall,wallidx_intersect_upper,wallidx_intersect_lower,backward=True)
        
        poly.add_vertex(surface.vertices[-1])

        return poly
    
    def genpoly_from_two_open_surfs(inner,outer,wall):

        poly = Polygon(increment=False)
        poly.add_whole_surface(outer)

        wallidx_intersect_br = poly.get_next_vertex_extrapolated_to_wall(wall)
        dummy = Polygon(increment=False)
        dummy.add_whole_surface(inner)
        wallidx_intersect_bl = dummy.get_next_vertex_extrapolated_to_wall(wall)

        poly.add_wall_segment(wall,wallidx_intersect_br,wallidx_intersect_bl,backward=True)

        poly.add_whole_surface(inner,backward=True)
        wallidx_intersect_ul = poly.get_next_vertex_extrapolated_to_wall(wall)
        dummy = Polygon(increment=False)
        dummy.add_whole_surface(outer,backward=True)
        wallidx_intersect_ur = dummy.get_next_vertex_extrapolated_to_wall(wall)

        poly.add_wall_segment(wall,wallidx_intersect_ul,wallidx_intersect_ur,backward=True)
       
        poly.add_vertex(outer.vertices[0])

        return poly



    # Creates the outer polygon from the solid wall Surface, writes it to file f
    # If debug, include lines that output polygons to poly.X.dat files
    def write_solid_polygon_dg2d(self,stratum,f,material,recyc,walltemp=300.0,debug=False,clockwise=False):
        f.write("new_zone solid\n")
        f.write("new_polygon\n")
        f.write("  stratum "+str(stratum)+"\n")
        f.write("  material "+material+"\n")
        f.write("  recyc_coef "+str(recyc)+"\n")
        f.write("  temperature "+str(walltemp)+"\n")
        f.write("  wall "+str(self.id)+" 0 1 \n ")
        f.write("  outer 1 0 \n")
        #f.write("  outer 0 1 \n")
        f.write("  wall "+str(self.id)+" 0 0 \n ")
        if debug:
            f.write("  print_polygon poly."+str(stratum)+".dat\n")
            f.write("  clear_polygon\n")
        else:
            f.write("  triangulate_polygon\n")
        f.write("\n")
 
        f.write("new_zone solid\n")
        f.write("new_polygon\n")
        f.write("  stratum "+str(stratum+1)+"\n")
        f.write("  material "+material+"\n")
        f.write("  recyc_coef "+str(recyc)+"\n")
        f.write("  temperature "+str(walltemp)+"\n")
        f.write("  wall "+str(self.id)+" 1 * \n ")
        f.write("  wall "+str(self.id)+" 0 0 \n ")
        f.write("  outer 1 2 3 4\n")
        f.write("  wall "+str(self.id)+" 1 1 \n ")
        if debug:
            f.write("  print_polygon poly."+str(stratum+1)+".dat\n")
            f.write("  clear_polygon\n")
        else:
            f.write("  triangulate_polygon\n")
        f.write("\n")

    def write_solid_polygon_with_exit(self,stratum,f,exitnodes,material,recyc,walltemp=300.0,debug=False):
        f.write("new_zone exit\n")
        f.write("new_polygon\n")
        f.write("  stratum "+str(stratum)+"\n")
        f.write("  wall "+str(self.id)+" "+str(exitnodes[0])+" "+str(exitnodes[1])+"\n")
        f.write("  outer 2 3 \n")
        if debug:
            f.write("  print_polygon poly."+str(stratum)+".dat\n")
            f.write("  clear_polygon\n")
        else:
            f.write("  triangulate_polygon\n")
        f.write("\n")
 
        f.write("new_zone solid\n")
        f.write("new_polygon\n")
        f.write("  stratum "+str(stratum+1)+"\n")
        f.write("  material "+material+"\n")
        f.write("  recyc_coef "+str(recyc)+"\n")
        f.write("  temperature "+str(walltemp)+"\n")
        f.write("  wall "+str(self.id)+" "+str(exitnodes[1])+" * \n")
        f.write("  wall "+str(self.id)+" 0 "+str(exitnodes[0])+"\n ")
        f.write("  outer 3\n")
        f.write("  outer 0 1 2\n")
        if debug:
            f.write("  print_polygon poly."+str(stratum+1)+".dat\n")
            f.write("  clear_polygon\n")
        else:
            f.write("  triangulate_polygon\n")
        f.write("\n")
        
    def plot_surfaces(surfs):
        for surf in surfs:
            R = []
            Z = []
            for vertex in surf.vertices:
                R.append(vertex.coords[0])
                Z.append(vertex.coords[1])
            plt.plot(R,Z,"+-")
            plt.plot(R[0],Z[0],"o")
        plt.show()

class Vertex:
    """
    A helper class to build polygons for definegeometry2d.

    Attributes:
        coords: 2-element array representing R,Z coordinates for the vertex.
        id: Integer identification number for the vertex. Must be passed to constructor.
        wall: Boolean flag for whether this vertex is "outer"; along the wall.
        wall_id: Id for the vertex in the context of the "wallfile". Inconsistent implementation?
    """
    def __init__(self,id_in,R,Z):
        """
        Constructor for Vertex class instance.

        Args:
            id_in: Mandatory integer identifier.
            R: Float R coordinate
            Z: Float Z coordinate
        """

        self.coords = [R,Z]
        self.id = id_in 
        self.wall = False
        self.wall_id = None

    # Nodes are the same the coordinates are equal
    def __eq__(self,other):
        if not isinstance(other,Vertex):
            return NotImplemented
        return (self.coords[0] == other.coords[0]) and (self.coords[1] == other.coords[1])

    def __ne__(self,other):
        if not isinstance(other,Vertex):
            return NotImplemented
        return not ( (self.coords[0] == other.coords[0]) and (self.coords[1] == other.coords[1]) )

def infer_wall_nodes(internal_triangles, external_triangles):
    wallnodes = []

    for tri_ext in external_triangles:
        for ext_node in tri_ext.vertices:
                for tri_int in internal_triangles:
                    if ext_node in tri_int.vertices and not ext_node in wallnodes:
                        wallnodes.append(ext_node)
                        wallnodes[-1].id = ext_node.id
                        ext_node.wall = True

#    # Remove duplicates
#    new_wallnodes = []
#    new_wallnodes.append(wallnodes[0])
#    for node in wallnodes:
#        duplicate = False
#        for newnode in new_wallnodes:
#            if node.id != newnode.id:
#                duplicate = True
#        if not duplicate:
#            new_wallnodes.append(node)
#            new_wallnodes[-1].id = node.id
#
#    return new_wallnodes
    return wallnodes

# x_corners and z_corners are 2D arrays of the same shape
def get_polys_from_corners(x_corners,z_corners):
    Nx=len(x_corners[:,0])
    Nz=len(x_corners[0,:])
    nodes = []
    for ix in range(0,Nx):
        for iz in range(0,Nz):
            nodes.append(Vertex(ix*Nz+iz,x_corners[ix,iz],z_corners[ix,iz]))
            if (ix == 0) or (ix == Nx) or (iz == 0) or (iz == Nz):
                nodes[-1].wall = True

    polys=[]
    centers = np.zeros(((Nx-1)*(Nz-1),2))
    for ix in range(0,Nx-1):
        for iz in range(0,Nz-1):
            poly=Polygon()
            if nodes[ix*Nz+iz].id != ix*Nz+iz:
                exit("Node index not what was expected")
            poly.add_vertex(nodes[ix*Nz+iz])
            poly.add_vertex(nodes[ix*Nz+iz+1])
            poly.add_vertex(nodes[(ix+1)*Nz+iz+1])
            poly.add_vertex(nodes[(ix+1)*Nz+iz])
            polys.append(poly)
            centers[ix*(Nz-1)+iz,0] = 0.25*(x_corners[ix,iz]+x_corners[ix+1,iz]+x_corners[ix+1,iz+1]+x_corners[ix,iz+1])
            centers[ix*(Nz-1)+iz,1] = 0.25*(z_corners[ix,iz]+z_corners[ix+1,iz]+z_corners[ix+1,iz+1]+z_corners[ix,iz+1])
    
    return nodes,polys, centers


def purge_invalid_nodes(allnodes,valid_tris):
    valid_nodes = []
    for node in allnodes:
        valid = False
        for tri in valid_tris:
            if node in tri.vertices:
                valid = True
        if valid:
            valid_nodes.append(node)
    return valid_nodes

def find_next_wall_node(current_node,prev_node,wallnodes,walltriangles,first):

    # Collect triangles that share current_node
    adjacent_wall_triangles = []
    for tri in walltriangles:
        if current_node in tri.vertices:
            adjacent_wall_triangles.append(tri)

    n_adjacent_triangles = len(adjacent_wall_triangles)

#    print("current_node = %d"%current_node.id)

    next_node = -1
    # Loop through the adjacent wall triangles that share the current node
    for tri in adjacent_wall_triangles:
#        print(tri.id)
        # Create list of adjacent triangles that excludes the one under consideration
        other_adjacent_triangles = copy.deepcopy(adjacent_wall_triangles)
        other_adjacent_triangles.remove(tri)
        
#        other_adjacent_triangle_vertices = []
#        other_adjacent_triangle_vertices = []
#        for other_triangle in other_adjacent_triangles:
#            other_adjacent_triangle_vertices.append(other_triangle.vertices[0])
#            other_adjacent_triangle_vertices.append(other_triangle.vertices[1])
#            other_adjacent_triangle_vertices.append(other_triangle.vertices[2])


        # Which vertex of this triangle is the current_node?
        vertex_idx = tri.vertices.index(current_node)

        # Line segments of this triangle that share this vertex
        segments = []
        segments.append( [current_node,tri.vertices[(vertex_idx+1)%3]] )
        segments.append( [current_node,tri.vertices[(vertex_idx+2)%3]] )

        # Loop through segments that share the current node
        for segment in segments:
            # If there are not any other wall triangles that share this segment, we have found a candidate

            segment_shared_with_another_triangle = False
            for othertri in other_adjacent_triangles:
                if segment[1] in othertri.vertices:
                    segment_shared_with_another_triangle = True
                    if not segment[0] in othertri.vertices[:]:
                        sys.exit("Something's very wrong in find_next_wall_node.")
            if (segment[1] in wallnodes) and (not segment_shared_with_another_triangle):
                if first:
                    # Ensure we start by going counterclockwise from the low field side
                    if segment[1].coords[1] < current_node.coords[1]:
                        next_node=segment[1]
                elif (segment[1] != prev_node):
                    if (next_node != -1):
                        sys.exit("Found multiple candidates for next_node. Logic of code fails.")
                    next_node = segment[1]

    if next_node == -1:
        sys.exit("Could not find a candidate next_node.")

    return next_node

# WallVertex is a special case of Vertex, where the vertex ID is an array
# that specifies the Surface ID and vertex ID within that surface. 
class WallVertex(Vertex):
    def __init__(self,wallid,vertexid,R,Z):
        self.id = [wallid,vertexid]
        super().__init__(self.id,R,Z)




