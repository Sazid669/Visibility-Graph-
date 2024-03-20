from matplotlib import pyplot as plt
from matplotlib.path import Path
import numpy as np
import sys
import csv
import math

class Point:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def distance(self, other):
        return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

    def as_numpy(self):
        return np.array([self.x, self.y])

    def line_distance(self, segment):
        return np.linalg.norm(np.cross(segment.p2.as_numpy() - segment.p1.as_numpy(), segment.p1.as_numpy() - self.as_numpy())) / np.linalg.norm(segment.p2.as_numpy() - segment.p1.as_numpy())

    def __str__(self):
        return f"({self.x:.2f}, {self.y:.2f})"

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def length(self):
        return math.sqrt(self.x**2 + self.y**2)

    def vector(self, other):
        return Point(other.x - self.x, other.y - self.y)

    def unit(self):
        mag = self.length()
        return Point(self.x / mag, self.y / mag) if mag > 0 else Point(0, 0)

    def scale(self, scalar):
        return Point(self.x * scalar, self.y * scalar)

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)

    def __truediv__(self, scalar):
        return Point(self.x / scalar, self.y / scalar)

    def __floordiv__(self, scalar):
        return Point(int(self.x / scalar), int(self.y / scalar))

    def __mul__(self, scalar):
        return Point(self.x * scalar, self.y * scalar)

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __eq__(self, other):
        return abs(self.x - other.x) < 0.0001 and abs(self.y - other.y) < 0.0001

class Segment:
    def __init__(self, p1=Point(), p2=Point()):
        self.p1 = p1
        self.p2 = p2

    @classmethod
    def from_angle_length(cls, p1=Point(), angle=0, length=1):
        x2 = p1.x + math.cos(angle) * length
        y2 = p1.y + math.sin(angle) * length
        return cls(p1, Point(x2, y2))

    def intersect(self, other):
        if ccw(self.p1, other.p1, other.p2) != ccw(self.p2, other.p1, other.p2) and ccw(self.p1, self.p2, other.p1) != ccw(self.p1, self.p2, other.p2):
            p = self.intersection_point(other)
            if p == self.p1 or p == self.p2:
                return [False, None]
            else:
                return [True, p]
        else:
            return [False, None]

    def intersection_point(self, other):
        xdiff = (self.p1.x - self.p2.x, other.p1.x - other.p2.x)
        ydiff = (self.p1.y - self.p2.y, other.p1.y - other.p2.y)
        div = xdiff[0] * ydiff[1] - xdiff[1] * ydiff[0]

        if div == 0:
            print("Something went wrong!")
            return None

        d = (self.p1.x * self.p2.y - self.p1.y * self.p2.x, other.p1.x * other.p2.y - other.p1.y * other.p2.x)
        x = (d[0] * xdiff[0] - d[1] * xdiff[1]) / div
        y = (d[0] * ydiff[0] - d[1] * ydiff[1]) / div
        return Point(x, y)

    def magnitude(self):
        return math.sqrt((self.p1.x - self.p2.x)**2 + (self.p1.y - self.p2.y)**2)

    def dot(self, other):
        vect_a = self.p1.vector(self.p2)
        vect_b = other.p1.vector(other.p2)
        return vect_a.x * vect_b.x + vect_a.y * vect_b.y

def ccw(A, B, C):
    return (C.y - A.y) * (B.x - A.x) >= (B.y - A.y) * (C.x - A.x)

class RPS:
    def __init__(self, filename):
        self.vertices = []
        self.edges = []
        self.start = []
        self.goal = []
        self.visible_segments = []
        self.vertex_list = []
        self.edge_active_list = []
        self.visibility_graph = []

        with open(filename, newline='\n') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            for vertex in reader:
                self.vertices.append([eval(value) for value in vertex])

        self.start = self.vertices[0]
        self.goal = self.vertices[-1]

        self.closed_vertices = []
        initial_check = True
        with open(filename, newline='\n') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)
            for vertex in reader:
                self.closed_vertices.append([eval(value) for value in vertex])
                if initial_check:
                    polygon_id = self.closed_vertices[-1][0]
                    temp = [eval(value) for value in vertex]
                    initial_check = False
                if self.closed_vertices[-1][0] != polygon_id:
                    self.closed_vertices.insert(-1, temp)
                    temp = [eval(value) for value in vertex]
                    polygon_id = self.closed_vertices[-1][0]
            csvfile.close()
        self.closed_vertices.append(self.closed_vertices[-1][0])

    def is_edge(self, v_begin_id, v_end_id):
        index = self.closed_vertices.index(self.vertices[v_begin_id])
        if (self.closed_vertices[index + 1] == self.vertices[v_end_id] or
                self.closed_vertices[index - 1] == self.vertices[v_end_id]) and (self.vertices[v_begin_id][0] == self.vertices[v_end_id][0]):
            return True
        else:
            return False

    def check_visibility(self, v_begin_id, v_end_id):
        is_visible = False
        if self.is_edge(v_begin_id, v_end_id):
            is_visible = True
            return is_visible

        if (self.vertices[v_begin_id][0] == self.vertices[v_end_id][0]) and (not self.is_edge(v_begin_id, v_end_id)):
            angles = self.get_angles(v_begin_id)
            angle = angles[v_end_id]
            p1 = Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2])
            p2 = Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2])
            dist = p1.distance(p2) + 10

            edge_1 = Segment.from_angle_length(Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]),
                                               math.radians(angle + 10), dist)
            edge_2 = Segment.from_angle_length(Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]),
                                               math.radians(angle - 10), dist)

            polygon_edges = [Segment(Point(edge[0][1], edge[0][2]), Point(edge[1][1], edge[1][2])) for edge in self.edges]
            count = 0
            right = True
            left = True

            for i in range(len(polygon_edges)):
                if edge_1.intersect(polygon_edges[i])[0]:
                    if right and edge_1.intersect(polygon_edges[i])[1] != Point(self.vertices[v_end_id][1],
                                                                                  self.vertices[v_end_id][2]) and edge_1.intersect(
                            polygon_edges[i])[1] != Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]):
                        count += 1
                        right = False
                if left and edge_2.intersect(polygon_edges[i])[0]:
                    if edge_2.intersect(polygon_edges[i])[1] != Point(self.vertices[v_end_id][1],
                                                                     self.vertices[v_end_id][2]) and edge_2.intersect(
                            polygon_edges[i])[1] != Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]):
                        count += 1
                        left = False
            if count < 2:
                return True
            else:
                return False

        edge_ = Segment(Point(self.vertices[v_begin_id][1], self.vertices[v_begin_id][2]),
                        Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2]))
        for j in range(len(self.edges)):
            polygon_edge = Segment(Point(self.edges[j][0][1], self.edges[j][0][2]),
                                   Point(self.edges[j][1][1], self.edges[j][1][2]))
            if edge_.intersect(polygon_edge)[0]:
                if edge_.intersect(polygon_edge)[1] == Point(self.vertices[v_end_id][1], self.vertices[v_end_id][2]):
                    continue
                else:
                    return False
        return True

    def prune_graph(self, vis_graph):
        vis_graph_new = []
        for i in range(len(vis_graph)):
            if vis_graph_new.count(vis_graph[i]) == 0:
                vis_graph_new.append(vis_graph[i])
        return vis_graph_new

    def get_angles(self, current_vertex_id):
        angles = []
        a = Segment.from_angle_length(Point(self.vertices[current_vertex_id][1], self.vertices[current_vertex_id][2]), math.radians(0), 1)
        for i in range(len(self.vertices)):
            vertex_segment = Segment(Point(a.p1.x, a.p1.y), Point(self.vertices[i][1], self.vertices[i][2]))
            angle = math.degrees(
                math.atan2(-(vertex_segment.p1.y - vertex_segment.p2.y), -(vertex_segment.p1.x - vertex_segment.p2.x)) -
                math.atan2((a.p2.y - a.p1.y), (a.p2.x - a.p1.x))
            )
            if angle < 0:
                angle = angle + 360
            angles.append(angle)
        return angles

    def get_active_edge_list(self, horizontal_segment):
        for j in range(len(self.edges)):
            edge = Segment(Point(self.edges[j][0][1], self.edges[j][0][2]),
                          Point(self.edges[j][1][1], self.edges[j][1][2]))
            if horizontal_segment.intersect(edge)[0]:
                self.edge_active_list.append(edge)
        return self.edge_active_list

    def plot(self, visibility=[]):
        x = []
        y = []
        for i in visibility:
            x.append(self.vertices[i[0]][1])
            y.append(self.vertices[i[0]][2])
        for v in self.vertices:
            plt.plot(x, y, 'r+')
        for e in self.edges:
            plt.plot([e[0][1], e[1][1]], [e[0][2], e[1][2]], 'k')
        for e in visibility:
            plt.plot([self.vertices[e[0]][1], self.vertices[e[1]][1]],
                     [self.vertices[e[0]][2], self.vertices[e[1]][2]], 'g--')
        for i, v in enumerate(self.vertices):
            plt.text(v[1] + 0.2, v[2], str(i))
        plt.title("Visibility Graph Computed By RPS")
        plt.axis('equal')

    def compute_rps(self):
        self.edges = []
        for i in range(2, len(self.closed_vertices) - 2):
            if self.closed_vertices[i][0] == self.closed_vertices[i + 1][0]:
                self.edges.append([self.closed_vertices[i], self.closed_vertices[i + 1]])

        for i in range(len(self.vertices)):
            self.edge_active_list = []
            self.vertex_list = []
            angle = self.get_angles(i)
            s = np.array(angle)
            sorted_index = np.argsort(s)
            for index_ in sorted_index:
                self.vertex_list.append(self.vertices[index_])
            a = Segment.from_angle_length(Point(self.vertices[i][1], self.vertices[i][2]), math.radians(0), 100)
            active_edges = self.get_active_edge_list(a)

            for index_ in sorted_index:
                if index_ != i:
                    if self.check_visibility(i, index_):
                        if i < index_:
                            self.visibility_graph.append([i, index_])
                        else:
                            self.visibility_graph.append([index_, i])
                        a = Segment.from_angle_length(Point(self.vertices[i][1], self.vertices[i][2]), math.radians((angle[index_] + 10)), 100)
                        active_edges = self.get_active_edge_list(a)

        final_visibility_graph = self.prune_graph(self.visibility_graph)
        self.plot(final_visibility_graph)
        print("The Visible Nodes Computed By RPS are:")
        return final_visibility_graph

if __name__ == "__main":
    if len(sys.argv) != 2:
        print(f"Usage:\n {sys.argv[0]} file_name.csv")
    else:
        vis = RPS(sys.argv[1])
        v = vis.compute_rps()
        vis.plot(v)
        plt.show()
        print(v)
