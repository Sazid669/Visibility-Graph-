from matplotlib import pyplot as plt
from matplotlib import path
import numpy as np
import sys
import csv
import math

class Point:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y
    
    def dist(self, p):
        # Distance between point self and point p
        return math.sqrt((self.x - p.x)**2 + (self.y - p.y)**2)
    
    def numpy(self):
        # return the point (x, y) as a numpy array
        return np.array([self.x, self.y])
        
    def dist_line(self, l):
        # return the distance between point self an line l of type Segment.
        return np.linalg.norm(np.cross(l.p2.numpy() - l.p1.numpy(), l.p1.numpy() - self.numpy())) / np.linalg.norm(l.p2.numpy() - l.p1.numpy())

    def __str__(self):
        # returns point self as a string
        return "({}, {})".format(np.round(self.x, 2), np.round(self.y, 2))

    def dot(self, p):
        # Dot product
        return self.x * p.x + self.y*p.y

    def length(self):
        # returns modulus of point self
        return math.sqrt(self.x**2 + self.y**2)

    def vector(self, p):
        # creates a vector of type Point between point self and point p
        return Point(p.x - self.x, p.y - self.y)

    def unit(self):
        # makes the point self unitary if possible
        mag = self.length()
        if mag > 0:
            return Point(self.x/mag, self.y/mag)
        else:
            return Point(0, 0)

    def scale(self, sc):
        # multiplies point self by scalar sc
        return Point(self.x * sc, self.y * sc)

    def __add__(self, p):
        # add point self and point p component by component
        return Point(self.x + p.x, self.y + p.y)

    def __sub__(self, p):
        # substracts point self and point p component by component
        return Point(self.x - p.x, self.y - p.y)

    def __truediv__(self, s):
        # divides point self by scalar s
        return Point(self.x / s, self.y / s)
    
    def __floordiv__(self, s):
        # integer division of point self by scalar s
        return Point(int(self.x / s), int(self.y / s))
    
    def __mul__(self, s):
        return Point(self.x * s, self.y * s)
    
    def __rmul__(self, s):
        return self.__mul__(s)
    
    def __eq__(self, __o: object) -> bool:
        if abs(self.x - __o.x) < 0.0001 and abs(self.y - __o.y) < 0.0001:
            return True
        return False

def ccw(A, B, C):
    return (C.y - A.y) * (B.x - A.x) >= (B.y - A.y) * (C.x - A.x)

def det(a, b):
    return a[0] * b[1] - a[1] * b[0]

class Segment:
    def __init__(self, p1=Point(), p2=Point()):
        # A segment is defined by two Point objects
        self.p1 = p1
        self.p2 = p2

    @classmethod
    def point_angle_length(cls, p1=Point(), angle=0, length=1):
        # A segment can be initialized with a Point object, an angle, and a segment length.
        x2 = p1.x + math.cos(angle) * length
        y2 = p1.y + math.sin(angle) * length
        return cls(p1, Point(x2, y2))
        
    def intersect(self, s):
        # Return true if Segment self and Segment s intersect
        if ccw(self.p1, s.p1, s.p2) != ccw(self.p2, s.p1, s.p2) and ccw(self.p1, self.p2, s.p1) != ccw(self.p1, self.p2, s.p2):
            p = self.intersection_point(s)
            if p == self.p1 or p == self.p2:
                return [False, None]
            else:
                return [True, p]
        else:
            return [False, None]

    def intersection_point(self, line):
        # Returns the point in which line Segment self and line Segment s intersect
        xdiff = (self.p1.x - self.p2.x, line.p1.x - line.p2.x)
        ydiff = (self.p1.y - self.p2.y, line.p1.y - line.p2.y)

        div = det(xdiff, ydiff)
        if div == 0:
            print("Something went wrong!")
            return None

        d = (det((self.p1.x, self.p1.y), (self.p2.x, self.p2.y)), det((line.p1.x, line.p1.y), (line.p2.x, line.p2.y)))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return Point(x, y)
    
    def __str__(self):
        return "[{}, {}]".format(self.p1, self.p2)

    def magnitude_line(self):
        # Calculate the magnitude of the line.  
        return math.sqrt((self.p1.x - self.p2.x)**2 + (self.p1.y - self.p2.y)**2)

    def dot_segment(self, l):
        # Dot products of the lines
        A1= Point(self.p1.x,self.p1.y) 
        A2= Point(self.p2.x,self.p2.y) 
        vect_A=  A1.vector(A2)
        B1= Point(l.p1.x, l.p1.y) 
        B2= Point(l.p2.x,l.p2.y)
        vect_B=  B1.vector(B2) 
        return vect_A.x * vect_B.x + vect_A.y*vect_B.y

class RPS:

    def __init__(self, filename):
        self.vertices = []  # list of vertices .
        self.edges = []     # List of tuples containing the vertex index 
        self.start = []                    
        self.goal = []                 
        self.visble_segments = []
        self.vertex_list=[]
        self.edge_active_list=[]
        self.visibility_graph_=[]
        

        # Load the data from CSV and compute edges      
        with open(filename, newline='\n') as csvfile:
          reader = csv.reader(csvfile)
          header= next(reader)
          self.vertices = [[eval(value) for value in vertex] for vertex in reader]


        
        self.start = self.vertices[0]
        self.goal = self.vertices[-1]

        self.closed_vertices=[] # To store the closed-loop edges of obstacles.
        
        with open(filename, newline='\n') as csvfile:
          reader = csv.reader(csvfile)
          header = next(reader)
          
          self.closed_vertices = []
          intial_check = True
          polygon_id = None
          temp = None

          for vertex in reader:
              values = [eval(value) for value in vertex]
              self.closed_vertices.append(values)

              if intial_check:
                  polygon_id = values[0]
                  temp = values
                  intial_check = False

              if values[0] != polygon_id:
                  self.closed_vertices.insert(-1, temp)
                  temp = values
                  polygon_id = values[0]

          self.closed_vertices.append(self.closed_vertices[-1][0])
        
    def is_Edge(self, v_begin_id, v_end_id):# Verifying whether the segment forms an edge of a polygon.
      v_begin = self.vertices[v_begin_id]
      v_end = self.vertices[v_end_id]

      # Find the index of v_begin in closed_vertices
      index = self.closed_vertices.index(v_begin)

      # Check if v_end is adjacent to v_begin in closed_vertices and their x-coordinates are the same
      next_vertex = self.closed_vertices[(index + 1) % len(self.closed_vertices)]
      prev_vertex = self.closed_vertices[(index - 1) % len(self.closed_vertices)]

      if (next_vertex == v_end or prev_vertex == v_end) and v_begin[0] == v_end[0]:
          return True
      else:
          return False

    def check_visibility(self, v_begin_id, v_end_id):
        #Assess the visibility of a specific edge; return True if it's visible, or False if it's not.

        is_visible = False
        # Is (v_begin, v_end) an edge of an obstacle? If so, visibility is considered true.
        if self.is_Edge(v_begin_id, v_end_id):
          is_visible = True
          return is_visible
        # Part of the same polygon but not an actual edge.
        if self.vertices[v_begin_id][0]== self.vertices[v_end_id][0] and (not self.is_Edge(v_begin_id, v_end_id)):
          ang=self.get_angles(v_begin_id)
          angl=ang[v_end_id] #Determine the angle of the edge in relation to the horizontal line.
          p1= Point(self.vertices[v_begin_id][1],self.vertices[v_begin_id][2])
          p2= Point(self.vertices[v_end_id][1],self.vertices[v_end_id][2])
          dist_ = p1.dist(p2) + 10 # Extending the length of the edge to check for intersections with segments of the polygon.
          edge_1= Segment.point_angle_length(Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]),math.radians(angl + 10),dist_) # increasing the angle of the segment by 10 degrees to check if it will intersect polygon edge
          edge_2= Segment.point_angle_length(Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]),math.radians(angl - 10),dist_) # decreasing the angle of the segment by 10 degrees to check if it will intersect polygon edge
          ploygon_edge=[]
          for i in range(len(self.edges)): # Checking for intersections between edge_1 and edge_2 with all polygon edges, and storing these polygon edges in a list called polygon_edge.
            ploygon_edge.append(Segment(Point(self.edges[i][0][1],self.edges[i][0][2]) , Point(self.edges[i][1][1],self.edges[i][1][2])))
          count=0 # to count the number of intersections
          # To establish a directional intersection check process, meaning that if an intersection is detected to the right once, there is no need to check for intersections in the right direction again, and the same principle applies to the left direction
          right= True 
          left= True 
          for i in range(len(ploygon_edge)):
            if edge_1.intersect(ploygon_edge[i])[0]:
              if right and edge_1.intersect(ploygon_edge[i])[1] != Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2]) and edge_1.intersect(ploygon_edge[i])[1] != Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]):
                count=count+1
                right=False
            if left and edge_2.intersect(ploygon_edge[i])[0]:
              if edge_2.intersect(ploygon_edge[i])[1] != Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2]) and edge_2.intersect(ploygon_edge[i])[1] != Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]):
                count=count+1
                left = False
          if count < 2: # If there are either 0 or 1 recorded intersections, the visibility is considered True; otherwise, it is considered False.
            return True
          else:
            return False
      
        # When 'v_begin' and 'v_end' do not belong to the same polygon, verify if the edge (v_begin, v_end) intersects with any polygon edge. If it does, visibility is considered False; otherwise, it is considered True.
        
        edge_= Segment(Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]) , Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2])) # Segment connecting two nodes and we will check if this segment cuts anyother segment in the map
        for j in range(len(self.edges)):
          ploygon_edge= Segment(Point(self.edges[j][0][1],self.edges[j][0][2]) , Point(self.edges[j][1][1],self.edges[j][1][2]))
          if edge_.intersect(ploygon_edge)[0]:
            if edge_.intersect(ploygon_edge)[1] == Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2]):
              continue
            else:
              return False
        return True
  
    def get_angles(self, current_vertex_id):
    #Compute angles between a horizontal half line and line segments connecting all vertices in the map with the given current_vertex_id.
      a = Segment.point_angle_length(Point(self.vertices[current_vertex_id][1], self.vertices[current_vertex_id][2]), math.radians(0), 1)

      angles = [0] * len(self.vertices)
      for i in range(len(self.vertices)):
          vertex_segment = Segment(Point(a.p1.x, a.p1.y), Point(self.vertices[i][1], self.vertices[i][2]))
          delta_y = -(vertex_segment.p1.y - vertex_segment.p2.y)
          delta_x = -(vertex_segment.p1.x - vertex_segment.p2.x)
          atan2_result = math.degrees(math.atan2(delta_y, delta_x) - math.atan2(a.p2.y - a.p1.y, a.p2.x - a.p1.x))
          angles[i] = atan2_result + 360 if atan2_result < 0 else atan2_result
      return angles
    
    def prune_graph(self, vis_graph):
      #Remove edges that are already in the visibility graph but inverted e.g., (0, 1) and (1, 0)
      return [edge for i, edge in enumerate(vis_graph) if edge not in vis_graph[:i]]
    
    def active_edge_list(self, hor_segment):
    #Get the list of active edges intersecting the given horizontal segment.
      self.edge_active_list = [
          Segment(Point(edge[0][1], edge[0][2]), Point(edge[1][1], edge[1][2]))
          for edge in self.edges
          if hor_segment.intersect(Segment(Point(edge[0][1], edge[0][2]), Point(edge[1][1], edge[1][2])))[0]
      ]
      return self.edge_active_list

    def plot(self, visibility=[]):
        # Given the set of `visibility` edges, plots the visibility graph. Already done.
        x=[]
        y=[]
        for i in visibility:
          x.append(self.vertices[i[0]][1])
          y.append(self.vertices[i[0]][2])
        for v in self.vertices:
            plt.plot(x, y, 'r+')
        for e in self.edges:
            plt.plot([e[0][1],e[1][1]], 
                     [e[0][2],e[1][2]], 'k')

        for e in visibility:
            plt.plot([self.vertices[e[0]] [1], self.vertices[e[1]][1]], 
                     [self.vertices[e[0]][2], self.vertices[e[1]][2]], 'g--')
    
        for i, v in enumerate(self.vertices):
            plt.text(v[1] + 0.2, v[2], str(i))
        plt.title("Visbility Graph Computed By RPS")
        plt.axis('equal') 
    
    def compute_RPS(self):
    # Extract polygon edges (GIVEN in the topological map)
      self.edges = []
      for i in range(2, len(self.closed_vertices) - 2):
          if self.closed_vertices[i][0] == self.closed_vertices[i + 1][0]:
              self.edges.append([self.closed_vertices[i], self.closed_vertices[i + 1]])

      # Iterate over vertices and compute visibility graph
      self.visibility_graph_ = []
      for i in range(len(self.vertices)):
          angle = self.get_angles(i)  # Calculate the angle to the horizontal line
          sorted_indices = sorted(range(len(angle)), key=lambda x: angle[x])

          a = Segment.point_angle_length(Point(self.vertices[i][1], self.vertices[i][2]), math.radians(0), 100)

          active_edges = []
          for index_ in sorted_indices:
              if index_ != i:
                  if self.check_visibility(i, index_):
                      if i < index_:
                          self.visibility_graph_.append([i, index_])
                      else:
                          self.visibility_graph_.append([index_, i])

                  a = Segment.point_angle_length(Point(self.vertices[i][1], self.vertices[i][2]), math.radians(angle[index_] + 10), 100)
                  active_edges = self.active_edge_list(a)

      return self.prune_graph(self.visibility_graph_)

    

    

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:\n ${} file_name.csv".format(sys.argv[0]))
    else:
        vis = RPS(sys.argv[1])
        v = vis.compute_RPS()
        vis.plot(v)
        plt.show()
        print(v)