
from typing import List, Union, Optional

from scheme_display import setup_display, add_frame, plot_lines

from functools import reduce

import itertools

import math

# CLASS DEFINITIONS ------------------------------------------------------------

# Set a cutoff to avoid missmatches in natural offsets when operating with float values
# e.g. a resolution of 1000 means that different points must be in coordinates separated by 0.001
resolution = 1000

# An x,y coordinate
class Point:

    def __init__(self, x : float, y : float):
        # Apply here the resolution cutoff
        self.x = round( x * resolution ) / resolution
        self.y = round( y * resolution ) / resolution

    def __str__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __repr__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        return False

    def __hash__(self):
        return hash((self.x, self.y))
    
    # Point + Vector -> Point, Point + Point -> Vector
    def __add__(self, other):
        if isinstance(other, self.__class__):
            return Vector(self.x + other.x, self.y + other.y)
        if isinstance(other, Vector):
            return Point(self.x + other.x, self.y + other.y)

    # Point - Vector -> Point, Point - Point -> Vector
    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return Vector(self.x - other.x, self.y - other.y)
        if isinstance(other, Vector):
            return Point(self.x - other.x, self.y - other.y)

    # get the distance from this point to other specified point
    def get_distance (self, other):
        x_distance = self.x - other.x
        y_distance = self.y - other.y
        return math.sqrt( x_distance**2 + y_distance**2 )

class Vector:

    def __init__(self, x : float, y : float):
        self.x = x
        self.y = y

    def __str__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __repr__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        return False

    def __hash__(self):
        return hash((self.x, self.y))

    def __neg__(self):
        return Vector(-self.x, -self.y)
    
    # Vector + Vector -> Vector, Vector + Point -> Point
    def __add__(self, other):
        if isinstance(other, self.__class__):
            return Vector(self.x + other.x, self.y + other.y)
        if isinstance(other, Point):
            return Point(self.x + other.x, self.y + other.y)

    # Vector - Vector -> Vector, Vector - Point -> Point
    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return Vector(self.x - other.x, self.y - other.y)
        if isinstance(other, Point):
            return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, num):
        if isinstance(num, int) or isinstance(num, float):
            x_mul = self.x * num
            y_mul = self.y * num
            return Vector(x_mul, y_mul)
        raise NameError('Vector multiplication is only supported with numeric values')

    __rmul__ = __mul__

    def __div__(self, num):
        if isinstance(num, int) or isinstance(num, float):
            x_div = self.x / num
            y_div = self.y / num
            return Vector(x_div, y_div)
        raise NameError('Vector division is only supported with numeric values')

    __rdiv__ = __div__

    def __truediv__(self, num):
        if isinstance(num, int) or isinstance(num, float):
            x_div = self.x / num
            y_div = self.y / num
            return Vector(x_div, y_div)
        raise NameError('Vector division is only supported with numeric values')

    __rtruediv__ = __truediv__

    def get_magnitude(self) -> float:
        return math.sqrt( self.x**2 + self.y**2 )

    def normalized(self):
        magnitude = self.get_magnitude()
        return self / magnitude

    def get_direction(self) -> float:
        return math.sqrt( self.x**2 + self.y**2 )
        
# A segment defined by 2 coordinates (Points): 'a' and 'b'
class Line:

    def __init__(self, a : Point, b : Point, color : str = 'black'):
        if a == b:
            raise NameError('ERROR: Inexistent line. Points "a" and "b" must be different')
        self.a = a
        self.b = b
        self.color = color
        self.vector = b - a
        self.length = a.get_distance(b)

    def __str__(self):
        return 'A: ' + str(self.a) + ', B: ' + str(self.b)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.a == other.a and self.b == other.b) or (self.a == other.b and self.b == other.a)
        return False

    def __hash__(self):
        return hash(self.get_hash())

    def __contains__(self, other):
        if other.x and other.y:
            distance1 = self.a.get_distance(other)
            distance2 = self.b.get_distance(other)
            return distance1 + distance2 == self.length
        return False

    # Find out if the line is diagonal
    def is_diagonal(self):
        if self.a.x == self.b.x:
            return False
        if self.a.y == self.b.y:
            return False
        return True

    # Get a value which is always the same no matter the order of a and b
    def get_hash(self):
        def sort_by_x(point):
            return point.x
        def sort_by_y(point):
            return point.y
        points = [self.a, self.b]
        sorted_points = sorted( sorted(points, key=sort_by_x), key=sort_by_y )
        return tuple(sorted_points)

    # Check if another line has the same direction than this line
    def same_direction_as (self, line) -> bool:
        nvector1 = self.vector.normalized()
        nvector2 = line.vector.normalized()
        return nvector1 == nvector2 or nvector1 == -nvector2

    # Check if two lines form a corner
    # i.e. only one of their points is the same and both lines have different direction
    def makes_corner_with(self, line):
        if self == line:
            return False
        if self.same_direction_as(line):
            return False
        return line.a == self.a or line.b == self.a or line.a == self.b or line.b == self.b

    def split_at_points(self, points : List[Point]) -> list:
        # Get only thouse points which are cutting the line
        cutting_points = [ point for point in points if point in self and point != self.a and point != self.b ]
        # If no points are cutting the line then return only the intact line
        if len(cutting_points) == 0:
            yield self
            return
        # Sort points by distance with the line 'a' point
        def by_distance(point):
            return self.a.get_distance(point)
        sorted_points = sorted(cutting_points, key=by_distance)
        # Nex create all possible lines with these points
        for a, b in pairwise([self.a, *sorted_points, self.b]):
            yield Line(a, b)

    # Get the intersection between two lines
    # https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines
    # The 'in_extremis' sets if the 'a' and 'b' points which define lines are considered
    # When false, if the intersection point is one of these points then it will be ignored
    def get_intersection_point (self, line, in_extremis : bool = True) -> Optional[Point]:
        xdiff = Vector(self.a.x - self.b.x, line.a.x - line.b.x)
        ydiff = Vector(self.a.y - self.b.y, line.a.y - line.b.y)

        def det(a, b):
            return a.x * b.y - a.y * b.x

        div = det(xdiff, ydiff)
        # Lines are paralel
        if div == 0:
            return None

        d = Vector(det(self.a, self.b), det(line.a, line.b))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        intersection_point = Point(x, y)

        # WARNING: Till this point, we have taken lines as infinite lines, not segments
        # WARNING: Two lines may not intersect, but this function will return the hipotetic intersection point if both lines where infinite
        # Now we must verify that the intersection point is in both lines
        if not intersection_point in self or not intersection_point in line:
            return None

        # In case the 'in_extremis' argument is false check the intersection point to not be one of the 'a' or 'b' points from any line
        if in_extremis == False and intersection_point in [ self.a, self.b, line.a, line.b ]:
            return None

        #print(str(self) + ' / ' + str(line) + ' -> ' + str(Point(x, y)))
        return Point(x, y)

# A rectangular area defined by 2 coordinates (Points): 'max' and 'min'
class Rect:

    def __init__(self, pmin : Point, pmax : Point):
        self.pmin = pmin
        self.pmax = pmax
        self.lines = self.get_lines()
        self._area = None

    # Set the rect from lines
    # The new rect will contain all rects
    @classmethod
    def from_lines(cls, lines : List[Line]):
        # Check that there are at least 2 lines
        if len(lines) < 2:
            raise NameError('It is required at least 2 lines')
        # Get all line points and find the minimum and maximum x and y values of all those points
        points = [ point for line in lines for point in (line.a, line.b) ]
        x_coords = [ point.x for point in points ]
        x_min = min(x_coords)
        x_max = max(x_coords)
        y_coords = [ point.y for point in points ]
        y_min = min(y_coords)
        y_max = max(y_coords)
        # If any minimum and maximum values match then the rectangle has no area
        if x_min == x_max or y_min == y_max:
            raise NameError('The rectangle has no area')
        # Set pmin and pmax and build the rectange
        pmin = Point(x_min, y_min)
        pmax = Point(x_max, y_max)
        return cls(pmin, pmax)

    def __str__(self):
        return 'Min: ' + str(self.pmin) + ', Max: ' + str(self.pmax)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.pmin == other.pmin and self.pmax == other.pmax
        return False

    def __hash__(self):
        return hash((self.pmin, self.pmax))

    def __contains__(self, other):
        if isinstance(other, Point):
            in_x = other.x >= self.pmin.x and other.x <= self.pmax.x
            in_y = other.y >= self.pmin.y and other.y <= self.pmax.y
            return in_x and in_y
        if isinstance(other, Line):
            in_a = other.a in self
            in_b = other.b in self
            return in_a and in_b
        if isinstance(other, Rect):
            in_pmin = other.pmin in self
            in_pmax = other.pmax in self
            return in_pmin and in_pmax
        return False

    # Get rectangle corners appart from pmin and pmax
    def get_upper_left_corner (self):
        return Point(self.pmin.x, self.pmax.y)
    def get_bottom_right_corner (self):
        return Point(self.pmax.x, self.pmin.y)

    # Return all rectangle points in a 'perimeter-friendly' order
    def get_points(self):
        point1 = self.pmin
        point2 = self.get_upper_left_corner()
        point3 = self.pmax
        point4 = self.get_bottom_right_corner()
        return [ point1, point2, point3, point4 ]

    # Return all rectangle lines in a 'perimeter-friendly' order
    def get_lines(self):
        points = self.get_points()
        lines = [ Line(a,b) for a, b in pairwise(points, retro=True) ]
        return lines

    # Return all rectangle lines in a 'perimeter-friendly' order
    def get_corssing_line(self):
        return Line(self.pmin, self.pmax)

    # Calculate the rectangle area
    def get_area(self):
        if self._area:
            return self._area
        x_size = self.pmax.x - self.pmin.x
        y_size = self.pmax.y - self.pmin.y
        self._area = x_size * y_size
        return self._area

    # The area is treated appart since it may be an expensive calculation
    area = property(get_area, None, None, "The rectangle area")

    # Split a rect in as many rects as specified by the 'x' and 'y' cuts
    # DANI: No lo he provado desde que lo moví de abajo
    def split (self, x_splits : list = [], y_splits : list = []):
        # Get the rectangle minimum and maximum values
        pmin = self.pmin
        pmax = self.pmax
        # Sort the split values and discard those values out of range
        def x_in_range(x):
            return x > pmin.x and x < pmax.x
        formatted_x_splits = list(filter(x_in_range, sorted(x_splits)))
        def y_in_range(y):
            return y > pmin.y and y < pmax.y
        formatted_y_splits = list(filter(y_in_range, sorted(y_splits)))
        # Set the steps to form rectangles
        # WARNING: The doble 'for' loop with generators seems to work wrong
        x_steps = list(pairwise([ pmin.x, *formatted_x_splits, pmax.x ]))
        y_steps = list(pairwise([ pmin.y, *formatted_y_splits, pmax.y ]))
        # Create as many rectangles as required
        for xmin, xmax in x_steps:
            for ymin, ymax in y_steps:
                pmin = Point(xmin, ymin)
                pmax = Point(xmax, ymax)
                yield Rect(pmin, pmax)

    # Given 2 rectangles, it returns the overlapping region, if exists, as a new rectangle
    # DANI: No lo he provado desde que lo moví de abajo
    def overlap_rect (self, rect):
        # Find the overlap in the 'x' dimension
        # Get the maximum of the minimums
        x_minimum = max(self.pmin.x, rect.pmin.x)
        # Get the minimum of the maximums
        x_maximum = min(self.pmax.x, rect.pmax.x)
        # Check that the overlap range exists
        x_overlap = ( x_minimum, x_maximum )
        if x_minimum >= x_maximum:
            return None
        # Find the overlap in the 'y' dimension
        # Get the maximum of the minimums
        y_minimum = max(self.pmin.y, rect.pmin.y)
        # Get the minimum of the maximums
        y_maximum = min(self.pmax.y, rect.pmax.y)
        # Check that the overlap range exists
        y_overlap = ( y_minimum, y_maximum )
        if y_minimum >= y_maximum:
            return None
        # Build a rectangle with both dimensional overlaps
        pmin_overlap = Point(x_overlap[0], y_overlap[0])
        pmax_overlap = Point(x_overlap[1], y_overlap[1])
        overlap = Rect(pmin_overlap, pmax_overlap)
        return overlap

    # Join 2 rectangles by returning a list of all rectangles which define the resulting polygon
    # The overlapped region, if exists is transformed to a single rectangle
    # DANI: No lo he provado desde que lo moví de abajo
    def join_rect (self, rect) -> Optional[list]:
        # Find the overlap between these two rectangles
        overlap = self.overlap_rect(rect)
        # If there is no overlap then just return both input rectangles
        if not overlap:
            return [self, rect]
        # Split both input rectangles using the overlap maximum and minimum points as split points
        x_splits = [ overlap.pmin.x, overlap.pmax.x ]
        y_splits = [ overlap.pmin.y, overlap.pmax.y ]
        split1 = self.split(x_splits, y_splits)
        split2 = rect.split(x_splits, y_splits)
        rects = list(split1) + list(split2)
        # Return only unique rectangles
        return set(rects)

    # Substract the second rectangle form the first rectangle and return the rectangles which define the resulting polygon
    # DANI: No lo he provado desde que lo moví de abajo
    def subtract_rect (self, rect) -> Optional[list]:
        # Find the overlap between these two rectangles
        overlap = selfoverlap_rect(rect)
        # If there is no overlap then just return the first rectangle intact
        if not overlap:
            return [self]
        # Split the first input rectangle using the overlap maximum and minimum points as split points
        x_splits = [ overlap.pmin.x, overlap.pmax.x ]
        y_splits = [ overlap.pmin.y, overlap.pmax.y ]
        split = self.split(x_splits, y_splits)
        # Return all splits but the overlapped rectangle
        rects = [ rect for rect in list(split) if rect != overlap ]
        return rects

# A perimeter defined by several lines
# The perimeter parameters must never be modified
# Instead, a new perimeter must be created every time a perimeter needs to be modified
class Perimeter:

    def __init__(self, lines : list):
        self.lines = lines
        # Check the perimeter is closed and lines are not diagonal
        self.check() 
        self._corners = None
        self._rects = None
        self._area = None

    # Set the perimeter from its corners
    # The new rect will contain all rects
    @classmethod
    def from_corners(cls, corners : List[Point]):
        # Check that there are at least 4 points
        if len(corners) < 4:
            raise NameError('It is required at least 4 points')
        # Set the line between each pair of points
        lines = []
        for a,b in pairwise(corners, retro=True):
            new_line = Line(a,b)
            lines.append(new_line)
        # Build the rectangle
        return cls(lines)

    def __str__(self):
        return ', '.join([str(line) for line in self.lines])

    def __contains__(self, other):
        if isinstance(other, Point):
            for rect in self.rects:
                if other in rect:
                    return True
            return False
        if isinstance(other, Line):
            in_a = other.a in self
            cross_any_line = any([ other.get_intersection_point(line, in_extremis=False) for line in self.lines ])
            return in_a and not cross_any_line
        if isinstance(other, Rect):
            return all([ line in self for line in other.get_lines() ])
        if isinstance(other, self.__class__):
            return all([ rect in self for rect in other.rects ])
        return False

    # Get the perimeter corners
    def get_corners(self):
        # If rects are previously calculated then return them
        if self._corners:
            return self._corners
        # Get all rectangles which form the perimeter
        corners = self.set_corners()
        self._corners = corners
        return self._corners

    # The area is treated appart since it may be an expensive calculation
    corners = property(get_corners, None, None, "The perimeter corners")

    # Get all rectangles which form the perimeter
    def get_rects(self):
        # If rects are previously calculated then return them
        if self._rects:
            return self._rects
        # Get all rectangles which form the perimeter
        rects = self.split_in_rectangles()
        self._rects = rects
        return self._rects

    # The area is treated appart since it may be an expensive calculation
    rects = property(get_rects, None, None, "The rectangles which form the perimeter")

    # Get the area inside the perimeter
    def get_area(self):
        # If the area is previously calculated then return it
        if self._area:
            return self._area
        # Otherwise, calculate the area
        # Get all rectangles which form the perimeter
        rects = self.rects
        # Get the area of each rectangle
        areas = [ rect.get_area() for rect in rects ]
        self._area = sum(areas)
        return self._area

    # The area is treated appart since it may be an expensive calculation
    area = property(get_area, None, None, "The area inside the perimeter")

    # Check if the current lines create a closed perimeter or if it is open
    def is_closed(self):
        # Check that each line ends in the same point that the next line starts
        for current, nextone in pairwise(self.lines, retro=True):
            if current.b != nextone.a:
                return False
        return True

    # Check the perimeter to be closed
    # Otherwise return an error, since not closed perimeters are not supported
    # DANI: En principio los perímetros no cerrados no tendrán soporte nunca porque no tienen mucho sentido o no les veo la utilidad
    # DANI: Esto es para el desarrollo. Una vez esté comprobado que los perímteros son estables quitaré esto porque retrasa el cálculo
    # Check each perimeter line to not be diagonal
    # Otherwise return an error, since perimeters with diagonal lines are not yet supported
    # DANI: En principio algun día se podría hacer esto
    def check(self):
        if not self.is_closed():
            raise NameError('The perimeter is not closed')
        for line in self.lines:
            if line.is_diagonal():
                raise NameError('The perimeter has diagonal lines')

    # Set the perimeter corners as points with additional stored values
    # The 'lines' variable defines the two lines of the perimeter which form the corner itself
    # The 'inside' variable defines if the corner is "pointing" to the inside of the perimeter
    # (e.g. rectangles have no inside corners while an 'L' has one inside corner)
    def set_corners(self):
        # For each line, save the end point and if the direction of the next line is left
        # Count how many corners are there in one direction (left in this case)
        precorners = []
        left_count = 0
        for current, nextone in pairwise(self.lines, retro=True):
            point = current.b
            # Given two continue lines which inevitably form a corner, set if the second line goes "left" (true) or "right" (false) relative to the first line
            # https://stackoverflow.com/questions/13221873/determining-if-one-2d-vector-is-to-the-right-or-left-of-another
            a = current.vector
            b = nextone.vector
            dot = a.x * -b.y + a.y * b.x
            lefted_corner = dot < 0
            left_count += int(lefted_corner)
            precorners.append((point, [current, nextone], lefted_corner))

        # Get the numeric difference between corners in one direction and the other
        difference = 2*left_count - len(precorners)

        # There should be always 4 more corners in one direction than in the other (may be left or right)
        # DANI: Esto es para el desarrollo. Una vez esté comprobado que los perímteros son estables quitaré esto porque retrasa el cálculo
        if abs(difference) != 4:
            raise NameError('There is something wrong with the perimeter')

        # Check if are more corners in the counted direction (true) or the other (false)
        lefted_perimeter = difference > 0

        # Now set the real corners and save them in the global variable 'corners'
        # Corners with the minoritarian direction will be set as 'inside = True'
        corners = []
        for precorner in precorners:
            point = precorner[0]
            lines = precorner[1]
            lefted_corner = precorner[2]
            inside = lefted_corner != lefted_perimeter
            point.lines = lines
            point.inside = inside
            corners.append(point)
        return corners

    # Get all perimeter points
    # DANI: Obsoleto
    def get_points(self):
        return [ line.a for line in self.lines ]

    # Get a rectangle which contains the whole perimeter
    def get_box(self):
        points = self.get_points()
        x_coords = [ point.x for point in points ]
        y_coords = [ point.y for point in points ]
        x_min = min(x_coords)
        x_max = max(x_coords)
        y_min = min(x_coords)
        y_max = max(x_coords)
        pmin = Point(x_min, y_min)
        pmax = Point(x_max, y_max)
        return Rect(pmin, pmax)

    # Return all points in the perimeter lines which intersect with a given line
    # Points are sorted according to their distance with the 'a' point of the line (from less to more distance)
    # WARNING: Intersection of paralel (overlapping) lines is not detected
    # WARNING: If the origen point is in a corner/line it may return or not the origen as intersection point
    # DANI: No se ha probado a fondo
    # DANI: Actualmente NO está en uso
    def get_line_intersection_points(self, line : Line) -> Optional[list]:
        # Get the intersection point of the specfied line with each perimeter limit
        intersection_points = []
        for limit in self.lines:
            point = limit.get_intersection_point(line)
            if point:
                intersection_points.append(point)
        # Find out also if the lines intersects any corner
        for corner in self.corners:
            if corner in line:
                intersection_points.append(corner)
        # If no points are found return None
        if len(intersection_points) == 0:
            return None
        # Sort the points
        def by_distance(point):
            return line.a.get_distance(point)
        sorted_points = sorted(intersection_points, key=by_distance)
        return sorted_points

    # Return all lines inside the perimeter which intersect with a given line
    # Lines are sorted according to their distance with the 'a' point of the line (from less to more distance)
    # WARNING: Intersection of paralel (overlapping) lines is not detected
    # DANI: No se ha probado
    # DANI: No está acabado. No se contempla la posibilidad de que la linea que intersecta empieze dentro del perímetro
    # DANI: Actualmente NO está en uso
    def get_line_intersection_lines(self, line : Line) -> Optional[list]:
        # Get the intersection points
        sorted_points = self.get_line_intersection_points(line)
        # If no points are found return None
        points_count = len(sorted_points)
        if points_count == 0:
            return None
        # Create new lines with the intersecting points
        intersecting_lines = []
        for a, b in pairwise(sorted_points):
            intersecting_lines.append(Line(a,b))
        # If the last intersecting point has no pair then use the line 'b' point a the end of the last intersecting line
        if points_count % 2 == 1:
            last_intersection_point = sorted_points[-1]
            last_line_point = line.b
            last_intersection_line = Line(last_intersection_point, last_line_point)
            intersecting_lines.append(last_intersection_line)
        return intersecting_lines

    # Given a perimeter corner which is pointing inside,
    # Create two lines opposed to the corner lines but as long as required to cut the perimeter at onther line or corner
    # This two lines will always be inside the perimeter
    def get_corner_insider_lines(self, corner : Point) -> list:
        # Check that the corner exists
        if corner not in self.corners:
            raise NameError('This is not a perimeter corner')
        # Check that the corner is inside
        if not corner.inside:
            raise NameError('This is only supported por inside corners')
        # Get the length of the most large possible line in the perimeter
        max_length = self.get_box().get_corssing_line().length
        # Get the oppoiste lines to the corner lines but as long as the max_length
        # NEVER FORGET: The corner point is the entry line 'b' point and the exit line 'a' point
        entry_line, exit_line = corner.lines

        tracing1 = Line(corner, corner + entry_line.vector.normalized() * max_length, color='green')
        tracing2 = Line(corner, corner + (-exit_line.vector.normalized()) * max_length, color='green')

        insider_lines = []
        for line in [tracing1, tracing2]:

            # Get the intersection point of the specfied line with each perimeter limit
            intersection_points = []
            for limit in self.lines:
                # Skip the lines of the main corner
                if limit in corner.lines:
                    continue
                point = limit.get_intersection_point(line)
                if point:
                    intersection_points.append(point)
            # Find out also if the line intersects any corner
            for corn in self.corners:
                # Skip the main corner
                if corn == corner:
                    continue
                if corn in line:
                    intersection_points.append(corn)

            # Sort the points by distance
            def by_distance(point):
                return line.a.get_distance(point)
            sorted_points = sorted(intersection_points, key=by_distance)

            # The closest point will be the first point
            cut_point = sorted_points[0]

            insider_line = Line(line.a, cut_point, color='red')
            insider_lines.append(insider_line)

        #return [tracing1, tracing2]
        return insider_lines

    # Split the current perimeter in a list of rectangles
    # If the multisplit option is true it will return as many splitted rectangles as possible
    # Else, it will return the minimum splits to cover all the perimeter
    # (e.g. an 'L' will return 2 rectangles if multisplit is false but 3 rectangles if it is true)
    def split_in_rectangles(self) -> list:

        # Get all inside corners
        inside_corners = [ corner for corner in self.corners if corner.inside ]

        # If there are no inside corners it means our perimeter is a single rectangle
        if len(inside_corners) == 0:
            # Return only the box
            return [ self.get_box() ]

        # Get all inside separator lines
        inside_separators = []
        for corner in inside_corners:
            for line in self.get_corner_insider_lines(corner):
                inside_separators.append(line)

        #add_frame([ *self.lines, *inside_separators ])

        # Remove duplicates
        inside_separators = list(set(inside_separators))

        # Find all points where the inside separator lines intersect each other
        inside_intersections = []
        for line1, line2 in itertools.combinations(inside_separators, 2):
            intersection_point = line1.get_intersection_point(line2)
            if not intersection_point:
                continue
            # All inside corners will be found as intersection points, since their 2 lines intersect
            # We skip these points
            if intersection_point in inside_corners:
                continue
            inside_intersections.append(intersection_point)

        # Remove duplicates
        inside_intersections = list(set(inside_intersections))        

        #print('Intersections: ' + str(len(inside_intersections)))

        # Split the inside separator lines at the intersection points
        inside_lines = []
        for line in inside_separators:
            splitted_lines = line.split_at_points(inside_intersections)
            inside_lines += list(splitted_lines)

        #print('Inside lines: ' + str(len(inside_lines)))

        # Join all inside lines with the perimeter limits in a isngle list
        # WARNING: inside lines MUST be before limit lines
        all_lines = [ *inside_lines, *self.lines ]

        # Finally, for each inside line, try to find 2 rectangles
        # In theory:
        # - All inside lines will be connected to exactly 2 rectangles
        # - All inside rectangles will be defined by at least 1 inside line
        minimum_rectangles = []
        for line in inside_lines:
            this_line_rects = [] # Max 2
            count = 0
            for other_line in all_lines:
                if line.makes_corner_with(other_line):
                    new_rect = Rect.from_lines([ line, other_line ])
                    if new_rect not in this_line_rects:
                        #print(str(line) + ' / ' + str(other_line) + ' -> ' + str(count) + ': ' + str(new_rect))
                        count += 1
                        this_line_rects.append(new_rect)
                    if count == 2:
                        break
            minimum_rectangles += this_line_rects

        #print('Total rectangles: ' + str(len(list(set(minimum_rectangles)))))

        #frame_lines = []
        #for rect in list(set(minimum_rectangles)):
        #    frame_lines += rect.get_lines()
        #    frame_lines.append(rect.get_corssing_line())
        #add_frame([ *self.lines, *frame_lines ])

        # At this point we have the "minimum" rectangles
        # Now it is time to find the "maximum" rectangles

        # First, some functions are defined to find colliding rects
        def get_left_rect (rect : Rect):
            pmax = rect.get_upper_left_corner()
            for r in minimum_rectangles:
                if r.pmax == pmax:
                    return r
            return None
        def get_right_rect (rect : Rect):
            pmin = rect.get_bottom_right_corner()
            for r in minimum_rectangles:
                if r.pmin == pmin:
                    return r
            return None
        def get_upper_rect (rect : Rect):
            pmin = rect.get_upper_left_corner()
            for r in minimum_rectangles:
                if r.pmin == pmin:
                    return r
            return None
        def get_bottom_rect (rect : Rect):
            pmax = rect.get_bottom_right_corner()
            for r in minimum_rectangles:
                if r.pmax == pmax:
                    return r
            return None

        # One by one for each *available rectangle, where available rectangles are the minimum rectangles
        # Get as many rectanges as possible which are connected horizontally to the current rectangle
        # Get as many rows of rectangles as possible which are connected vertically to all previous rectangles
        # Consider all previous rectangles as a single rectange
        # Repeat in the inverse order (first vertically, then horizontally)
        # Remove all previous rectangles from the *available rectangles list
        maximum_rectangles = []
        available_rectangles = [ rect for rect in minimum_rectangles ]
        for available_rect in available_rectangles:
            # Set the first row
            row = [ available_rect ]
            # Append all rectangles at right from current rectangle to the row
            rightest = available_rect
            while (True):
                rightest = get_right_rect(rightest)
                if not rightest:
                    break
                row.append(rightest)
            # Append all rectangles at left from current rectangle to the row
            leftest = available_rect
            while (True):
                leftest = get_left_rect(leftest)
                if not leftest:
                    break
                row.append(leftest)
            # Set the group of rectangles to be joined
            group = [ rect for rect in row ]
            # If all rectangles in the row have a botton rectangle then add all those new rects to a new row
            # This new row is then added to the whole group of rectanges and used to find the next row
            while (True):
                new_row = [ get_bottom_rect(rect) for rect in row ]
                if not all(new_row):
                    break
                group += new_row
                row = [ *new_row ]
            # Create a new rect which contains all group rects
            # Find the most maximum pmax and the most minimum pmin
            def sort_by_x(point):
                return point.x
            def sort_by_y(point):
                return point.y
            pmax_points = [ rect.pmax for rect in group ]
            sorted_pmax_points = sorted( sorted(pmax_points, key=sort_by_x), key=sort_by_y )
            maximum_pmax = sorted_pmax_points[0]
            pmin_points = [ rect.pmin for rect in group ]
            sorted_pmin_points = sorted( sorted( pmin_points, key=sort_by_x, inverse=True ), key=sort_by_y, inverse=True )
            minimum_pmin = sorted_pmin_points[0]
            maximum_rect = Rect(minimum_pmin, maximum_pmax)
            # Add the new maximum rectnagle to the list and update the available rects list
            maximum_rectangles.append(maximum_rect)
            available_rectangles = [ rect not in group for rect in available_rectangles ]

        # DANI: Falta añadir el 'maximum_rectangles' aqui
        return minimum_rectangles

# Auxiliar functions ---------------------------------------------------------------

# Set a special iteration system
# Return each value of the array and the next value
# By default, the final array value is skipped, since it has no next value
# However, if the 'retro' argument is True, the final array value is returned with the first array value as the next value
# By default, values are returned as follows: A with B, B with C, C with D ...
# However, if the 'loyals' argument is True, values are returned as follows: A with B, C with D, E with F ...
def pairwise (values : list, retro : bool = False, loyals = False):
    last = len(values) - 1
    step = 2 if loyals else 1
    for i in range(0, last, step):
        yield values[i], values[i+1]
    if retro:
        yield values[last], values[0]