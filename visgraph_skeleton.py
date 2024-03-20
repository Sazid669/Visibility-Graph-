#!/usr/bin/python3

from matplotlib import pyplot as plt
from matplotlib import path
import numpy as np
import sys
import csv
import math

def visibility(csv_file):
    points = pd.read_csv(csv_file)
    x = points.iloc[:, 1]
    y = points.iloc[:, 2]
    points = []
    for i in range(len(x)):
        point = (x[i], y[i])
        points.append(point)
    v_start = points[0]
    v_goal = points[-1]
    print(v_start)
    print(v_goal)

    vis = points[1: (len(x) - 1)]
    print (vis)
    return vis

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
                return False, None
            else:
                return True, p
        else:
            return False, None

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





if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage:\n ${} file_name.csv".format(sys.argv[0]))
    else:
        # execute visibility graph algorithm
        print(":)")
