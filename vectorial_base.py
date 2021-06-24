from typing import List, Union, Optional

from scheme_display import add_frame

from functools import reduce

import itertools

from math import sqrt

import random
colors = [
    'black',
    #'white',
    'blue',
    'green',
    'yellow',
    'red',
    'pink',
    'orange',
    'brown',
    'cyan',
    'purple',
]

# CLASS DEFINITIONS ------------------------------------------------------------

# Float numbers are not allowed internally
number = Union[int, float]
def is_number(var):
    return isinstance(var, int) or isinstance(var, float)

# Set a resoltion limit
# Resolution is the number of decimals to take in count when comparing different coordinates
# When working with squared values (e.g. areas) the resolution should be decreased in 2 orders
base_resolution = 4
def resolute(num, resolution = base_resolution):
    res = 10**resolution
    return round(num * res) / res

# An x,y coordinate
class Point:

    def __init__(self, x : number, y : number):
        # Save the coordinates in Decimal format applying the precision limit
        self.x = resolute(x)
        self.y = resolute(y)
        # Save both coordinates as a tuple
        self.coords = self.x, self.y

    def __str__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __repr__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __eq__(self, other):
        # We also check if the instance belongs to parent classes since this function is inherited by Corner
        if isinstance(other, (self.__class__, *type(self).__bases__)):
            return self.x == other.x and self.y == other.y
        return False

    def __hash__(self):
        return hash((self.x, self.y))
    
    # Point + Vector -> Point, Point + Point -> Vector
    def __add__(self, other):
        if isinstance(other, self.__class__):
            return Vector(other.x - self.x, other.y - self.y)
        if isinstance(other, Vector):
            return Point(self.x + other.x, self.y + other.y)

    # Point - Vector -> Point, Point - Point -> Vector
    def __sub__(self, other):
        if isinstance(other, self.__class__):
            return Vector(self.x - other.x, self.y - other.y)
        if isinstance(other, Vector):
            return Point(self.x - other.x, self.y - other.y)

    # Get the distance from this point to other specified point
    def get_distance_to (self, other : 'Point') -> number:
        x_distance = self.x - other.x
        y_distance = self.y - other.y
        return sqrt( x_distance**2 + y_distance**2 )

class Vector:

    def __init__(self, x : number, y : number):
        # Save the coordinates in Decimal format applying the precision limit
        self.x = resolute(x)
        self.y = resolute(y)

    # Set the vector from a slope
    # If slope is None a vertical vector is returned
    @classmethod
    def from_slope(cls, slope : number):
        if slope == None:
            return cls(1,0)
        y, x = slope.as_integer_ratio()
        return cls(x,y)

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
        if is_number(num):
            x_mul = self.x * num
            y_mul = self.y * num
            return Vector(x_mul, y_mul)
        raise RuntimeError('Vector multiplication is only supported with numeric values')

    __rmul__ = __mul__

    def __div__(self, num):
        if is_number(num):
            x_div = self.x / num
            y_div = self.y / num
            return Vector(x_div, y_div)
        raise RuntimeError('Vector division is only supported with numeric values')

    __rdiv__ = __div__

    def __truediv__(self, num):
        if is_number(num):
            x_div = self.x / num
            y_div = self.y / num
            return Vector(x_div, y_div)
        raise RuntimeError('Vector division is only supported with numeric values')

    __rtruediv__ = __truediv__

    def get_magnitude(self) -> number:
        return sqrt( self.x**2 + self.y**2 )

    def get_slope(self) -> Optional[number]:
        if self.x == 0:
            return None
        return self.y / self.x

    # Return a new vector with identical direction and sense but magnitude = 1
    def normalized(self) -> 'Vector':
        magnitude = self.get_magnitude()
        return self / magnitude

    # Find out if the vector is totally vertical
    def is_vertical(self) -> bool:
        return self.x == 0
    # Find out if the segment is totally horizontal
    def is_horizontal(self) -> bool:
        return self.y == 0

    # Find out if the segment is diagonal
    def is_diagonal(self) -> bool:
        return not self.is_vertical() and not self.is_horizontal()

    # Find out if two vectors are equivalent
    # i.e. they have the same direction and magnitude, no matter the sense
    def is_equivalent_to (self, other : 'Vector') -> bool:
        same_slope = self.get_slope() == other.get_slope()
        same_magnitude = self.get_magnitude() == other.get_magnitude()
        return same_slope and same_magnitude

# A line defined by a point and a directional vector
class Line:

    def __init__(self, point : Point, vector : Vector):
        self.point = point
        self.vector = vector

    def __str__(self):
        return 'x,y = ' + str(self.point) + ' + t' + str(self.vector)

    def __repr__(self):
        return 'x,y = ' + str(self.point) + ' + t' + str(self.vector)

    # In order to compare lines we must use the slope and the intercept
    # This is because lines with different points and vectors may be the same line
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.same_line_as(other)
        return False

    def __contains__(self, other):
        if isinstance(other, Segment):
            return self.same_line_as(other)
        return False

    # Get the vector from the line
    # Return None in case the line is vertical
    def get_slope(self) -> Optional[number]:
        return self.vector.get_slope()

    # Get the line y intercept
    # i.e. the y coordinate from the point in the vertical line x = 0 where this line cuts
    # Return None in case the line is vertical
    def get_y_intercept (self) -> int:
        slope = self.get_slope()
        if slope == None:
            return None
        point = self.point
        y_intercept = point.y - slope * point.x
        return y_intercept

    # Get the line y intercept point
    # i.e. the point in the vertical line x = 0 where this line cuts
    # Return None in case the line is vertical
    def get_y_intercept_point (self) -> Point:
        y_intercept = self.get_y_intercept()
        if y_intercept == None:
            return None
        return Point(0, self.get_y_intercept())

    # Get the line x intercept
    # i.e. the x coordinate from the point in the horizontal line y = 0 where this line cuts
    # Return None in case the line is horizontal
    def get_x_intercept(self) -> int:
        slope = self.get_slope()
        if slope == 0:
            return None
        point = self.point
        if slope == None:
            return point.x
        y_intercept = self.get_y_intercept()
        x_intercept = - y_intercept / slope
        return x_intercept

    # Get the line x intercept point
    # i.e. the point in the horizontal line y = 0 where this line cuts
    # Return None in case the line is horizontal
    def get_x_intercept_point (self) -> Point:
        x_intercept = self.get_x_intercept()
        if x_intercept == None:
            return None
        return Point(self.get_x_intercept(), 0)

    # Get one intercept point
    # Get the y intercept by default
    # Get the x intercept in case it is a vertical line
    def get_intercept_point (self) -> Point:
        return self.get_y_intercept_point() if not self.is_vertical() else self.get_x_intercept_point()

    # Find out if the segment is totally vertical
    def is_vertical(self) -> bool:
        return self.vector.is_vertical()

    # Find out if the segment is totally horizontal
    def is_horizontal(self) -> bool:
        return self.vector.is_horizontal()

    # Find out if the segment is diagonal
    def is_diagonal(self) -> bool:
        return self.vector.is_diagonal()

    # Check if 2 lines are paralel
    def is_paralel_to (self, other : 'Line') -> bool:
        return self.get_slope() == other.get_slope()

    # Check if 2 lines are matematically identical
    # WARNING: This function is independent from __eq__ since it must be inherited by the Segment class
    def same_line_as (self, other : 'Line') -> bool:
        are_paralel = self.is_paralel_to(other)
        same_intercept = self.get_intercept_point() == other.get_intercept_point()
        return are_paralel and same_intercept
        
# A segment defined by 2 coordinates (Points): 'a' and 'b'
class Segment(Line):

    def __init__(self, a : Point, b : Point, color : str = 'black'):
        if a == b:
            raise RuntimeError('Inexistent segment. Points "a" and "b" must be different: ' + str(a))
        self.a = a
        self.b = b
        self.color = color
        super().__init__(a, b - a)
        self.length = a.get_distance_to(b)
        # Save both points as a tuple
        self.points = a, b

    def __str__(self):
        return 'A: ' + str(self.a) + ' -> B: ' + str(self.b)

    def __repr__(self):
        return 'A: ' + str(self.a) + ' -> B: ' + str(self.b)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.a == other.a and self.b == other.b) or (self.a == other.b and self.b == other.a)
        return False

    def __hash__(self):
        return hash(self.get_hash())

    def __contains__(self, other):
        if isinstance(other, Point):
            distance1 = self.a.get_distance_to(other)
            distance2 = self.b.get_distance_to(other)
            return resolute(distance1 + distance2) == resolute(self.length)
        if isinstance(other, self.__class__):
            return other.a in self and other.b in self
        return False

    # Get a value which is always the same no matter the order of a and b
    def get_hash(self) -> tuple:
        def sort_by_x(point):
            return point.x
        def sort_by_y(point):
            return point.y
        points = [self.a, self.b]
        sorted_points = sorted( sorted(points, key=sort_by_x), key=sort_by_y )
        return tuple(sorted_points)

    # Get the inverted segment (i.e. invert a and b points)
    def inverted (self):
        return Segment(self.b, self.a, self.color)

    # Create a new segment identical to this segment but with a specified color
    def get_colored_segment(self, color : str):
        return Segment(self.a, self.b, color)

    # Check if two segments have a point in common
    # WARNING: They may be overlapped
    def is_connected_with(self, other : 'Segment') -> bool:
        return other.a == self.a or other.b == self.a or other.a == self.b or other.b == self.b

    # Check if two segments form a corner
    # i.e. they have a point in common and both segments are not paralel
    def makes_corner_with(self, other : 'Segment') -> bool:
        if self.is_paralel_to(other):
            return False
        return self.is_connected_with(other)

    # Check if a line is fully covered by a list of segments
    # i.e. substract all segments from self and if there is still a remaining segment it is not covered
    def is_covered_by(self, segments : List['Segment']) -> bool:
        remaining_segments = self.substract_segments(segments)
        return len(remaining_segments) == 0

    # Split the segment at multiple points and return the splitted segments
    def split_at_points(self, points : List[Point]) -> list:
        # Get only thouse points which are cutting the segment
        cutting_points = [ point for point in points if point in self and point != self.a and point != self.b ]
        # If no points are cutting the segment then return only the intact segment
        if len(cutting_points) == 0:
            yield self
            return
        # Sort points by distance with the segment 'a' point
        def by_distance(point):
            return self.a.get_distance_to(point)
        sorted_points = sorted(cutting_points, key=by_distance)
        # Nex create all possible segments with these points
        for a, b in pairwise([self.a, *sorted_points, self.b]):
            yield Segment(a, b)

    # Get the intersection point between two segments
    # https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-segments
    # The 'in_extremis' sets if the 'a' and 'b' points which define segments are considered
    # in_extremis = 0 -> Intersections which are the 'extrem' point of any segment are ignored
    # in_extremis = 1 -> Intersections which are the 'extrem' point of only one of the segments are also considered
    # in_extremis = 2 -> All interactions are considered
    def get_intersection_point (self, segment, in_extremis : int = 2) -> Optional[Point]:
        xdiff = Vector(self.a.x - self.b.x, segment.a.x - segment.b.x)
        ydiff = Vector(self.a.y - self.b.y, segment.a.y - segment.b.y)

        def det(a, b):
            return a.x * b.y - a.y * b.x

        div = det(xdiff, ydiff)
        # segments are paralel
        if div == 0:
            return None

        d = Vector(det(self.a, self.b), det(segment.a, segment.b))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        intersection_point = Point(x, y)

        # WARNING: Till this point, we have taken segments as infinite segments, not segments
        # WARNING: Two segments may not intersect, but this function will return the hipotetic intersection point if both segments where infinite
        # Now we must verify that the intersection point is in both segments
        if not intersection_point in self or not intersection_point in segment:
            return None

        # In case the 'in_extremis' argument is 0 check the intersection point to not be one of the 'a' or 'b' points from any segment
        if in_extremis == 0 and intersection_point in [ self.a, self.b, segment.a, segment.b ]:
            return None

         # In case the 'in_extremis' argument is 1 check the intersection point to not be one of the 'a' or 'b' points from both segments
        if in_extremis == 1 and intersection_point in [ self.a, self.b ] and intersection_point in [ segment.a, segment.b ]:
            return None

        #print(str(self) + ' / ' + str(segment) + ' -> ' + str(Point(x, y)))
        return Point(x, y)

    # Get the point in the middle of the segment
    def get_middle_point (self) -> Point:
        return self.a + self.vector / 2

    # Get the overlap segment between two segments
    # Return None if segments are not in the same line
    def get_overlap_segment (self, other : 'Segment') -> Optional['Segment']:
        # If both segments are not in the same line then there is no overlap
        if not self.same_line_as(other):
            return None
        # Otherwise, order both segment points and check how they alternate
        self_points = [self.a, self.b]
        other_points = [other.a, other.b]
        points = list(set(self_points + other_points))
        def sort_by_x(point):
            return point.x
        def sort_by_y(point):
            return point.y
        sorted_points = sorted( sorted( points, key=sort_by_x), key=sort_by_y)
        # Track from which point both lines are in, if any
        # If both lines are in from 1 point then the overlap segment goes from this point to the next
        in_self = False
        in_other = False
        first_point = None
        for point in sorted_points:
            if point in self_points:
                in_self = not in_self
            if point in other_points:
                in_other = not in_other
            if first_point:
                return Segment(first_point, point)
            if in_self and in_other:
                first_point = point
        # If segments do not overlap at any moment then return None
        return None

    # Get current segment after substracting other segments
    # Return an empty list if self segment is fully substracted
    def substract_segments (self, others) -> list:
        # Filter others to be in the same line that self segment
        inline_segments = [ segment for segment in others if self.same_line_as(segment) ]
        # Order segment points and check how they alternate
        self_points = [self.a, self.b]
        other_points = []
        for segment in inline_segments:
            other_points += [segment.a, segment.b]
        points = list(set(self_points + other_points))
        def sort_by_x(point):
            return point.x
        def sort_by_y(point):
            return point.y
        sorted_points = sorted( sorted( points, key=sort_by_x), key=sort_by_y)
        # Create new segments between each pair of sorted points
        splitted_segments = [ Segment(a,b) for a, b in pairwise(sorted_points) ]
        # Get only segments which are in self and are not in others
        result_segments = []
        for segment in splitted_segments:
            if segment not in self:
                continue
            if any(segment in other for other in others):
                continue
            result_segments.append(segment)
        return result_segments

# A corner is a point where 2 non-paralel segments are connected
# i.e. both segments have this point ('a' or 'b') in common
# IMPORTANT: It is a standard that the first segment enters the corner while the second exits the corner
# i.e. first segment is 'a' -> (corner) and second segment is (corner) -> 'b'
# The inside argument is refered to if this corner points towards the 'inside' of the perimeter which belongs
class Corner(Point):

    def __init__(self, x : number, y : number, in_segment : Segment, out_segment : Segment, inside : bool = None):
        super().__init__(x, y)
        if in_segment.b != self or out_segment.a != self:
            raise RuntimeError('The corner does not follow the standard')
        self.in_segment = in_segment
        self.out_segment = out_segment
        # Save also both segments as a tuple
        self.segments = self.in_segment, self.out_segment
        # Save in addition the extreme points of segments
        # i.e. the in segment 'a' point and the out segment 'b' point
        self.in_point = in_segment.a
        self.out_point = out_segment.b
        # Save also both points as a tuple
        self.points = self.in_point, self.out_point
        self.inside = inside

    # Find out if this corner is continuous with other corner
    # Is so, return the merge of both continuous segments
    def get_continuous_segment_with(self, other):
        # If the corner point is not the same we are over
        if other != self:
            return None
        # Calculate the normalized a (from in) and b (from out) points
        # i.e. the a and b points if both segments length was 1
        # Segments which overlap will return identical points
        self_normalized_a = self - self.in_segment.vector.normalized()
        self_normalized_b = self + self.out_segment.vector.normalized()
        self_points = [ self_normalized_a, self_normalized_b ]
        other_normalized_a = other - other.in_segment.vector.normalized()
        other_normalized_b = other + other.out_segment.vector.normalized()
        other_points = [ other_normalized_a, other_normalized_b ]
        # Find self points which do not match any other point
        self_different = [ p for p, point in enumerate(self_points) if point not in other_points ]
        # If there are no identical points or both are identical we are overordered_segmentsnt not in self_points ]
        # And finally we check that both different points belong to aligned segments
        # For this we just check if the corner point is between them
        self_different_point = self.points[self_different[0]]
        other_different_point = other.points[other_different[0]]
        continuous_segment = Segment(self_different_point, other_different_point)
        if self not in continuous_segment:
            return None
        return continuous_segment

# A rectangular area defined by 2 coordinates (Points): 'max' and 'min'
class Rect:

    def __init__(self, pmin : Point, pmax : Point, segments_color : str = 'black', fill_color : str = 'white'):
        self.pmin = pmin
        self.pmax = pmax
        self.segments_color = segments_color
        self.fill_color = fill_color
        self.segments = self.get_segments()
        self._area = None

    # Set the rect from segments
    # The new rect will contain all rects
    @classmethod
    def from_segments(cls, segments : List[Segment], segments_color : str = 'black', fill_color : str = 'white'):
        # Check that there are at least 2 segments
        if len(segments) < 2:
            raise RuntimeError('It is required at least 2 segments')
        # Get all segment points and find the minimum and maximum x and y values of all those points
        points = [ point for segment in segments for point in (segment.a, segment.b) ]
        x_coords = [ point.x for point in points ]
        x_min = min(x_coords)
        x_max = max(x_coords)
        y_coords = [ point.y for point in points ]
        y_min = min(y_coords)
        y_max = max(y_coords)
        # If any minimum and maximum values match then the rectangle has no area
        if x_min == x_max or y_min == y_max:
            print(segments)
            raise RuntimeError('The rectangle has no area')
        # Set pmin and pmax and build the rectange
        pmin = Point(x_min, y_min)
        pmax = Point(x_max, y_max)
        return cls(pmin, pmax, segments_color, fill_color)

    # Set the rect from a corner (i.e. a point with 2 segments)
    # Optionally you can ask for specific x and y sizes
    # If no size is passed then the size of the original corner segment is used
    @classmethod
    def from_corner(cls, corner : Corner, x_size = None, y_size = None, segments_color : str = 'black', fill_color : str = 'white'):
        # Get the two segments from the corner
        segments = corner.segments
        directions = [ segment.vector.normalized() for segment in segments ]
        # NEVER FORGET: The first direction is the 'entring' segment so its vector points from away to the corner
        # NEVER FORGET: The second direction is the 'exiting' segment so its vector points from the corner to away
        # We want both directions 'exiting', so we change the direction of the first vector
        directions[0] = -directions[0]
        # Find out which is the horizontal and which is the vertical direction
        # If some of these steps fail it means there is no horizontal or/and vertical segments
        hdir = next( d for d, direction in enumerate(directions) if direction.is_horizontal() )
        vdir = next( d for d, direction in enumerate(directions) if direction.is_vertical() )
        # Finally set the rectangle segments
        # If the size of any direction is forced then use the vector to build a new segments
        # Otherwise use the original segments
        hsegment = Segment(corner, corner + directions[hdir] * x_size) if x_size else segments[hdir]
        vsegment = Segment(corner, corner + directions[vdir] * y_size) if y_size else segments[vdir]
        return cls.from_segments([hsegment, vsegment], segments_color, fill_color)

    def __str__(self):
        return 'Min: ' + str(self.pmin) + ', Max: ' + str(self.pmax)

    def __repr__(self):
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
        if isinstance(other, Segment):
            in_a = other.a in self
            in_b = other.b in self
            return in_a and in_b
        if isinstance(other, Rect):
            in_pmin = other.pmin in self
            in_pmax = other.pmax in self
            return in_pmin and in_pmax
        return False

    # Get rectangle points appart from pmin and pmax
    def get_upper_left_point (self) -> Point:
        return Point(self.pmin.x, self.pmax.y)
    def get_bottom_right_point (self) -> Point:
        return Point(self.pmax.x, self.pmin.y)

    # Return all rectangle points in a 'perimeter-friendly' order
    def get_points(self) -> list:
        point1 = self.pmin
        point2 = self.get_upper_left_point()
        point3 = self.pmax
        point4 = self.get_bottom_right_point()
        return [ point1, point2, point3, point4 ]

    # Return all rectangle segments in a 'perimeter-friendly' order
    def get_segments(self) -> list:
        points = self.get_points()
        segments = [ Segment(a, b, self.segments_color) for a, b in pairwise(points, retro=True) ]
        return segments

    # Return all rectangle points in a 'perimeter-friendly' order
    # Each point contains its two adjacent segments
    def get_corners(self) -> list:
        points = self.get_points()
        segments = [ Segment(a,b) for a, b in pairwise(points, retro=True) ]
        segment_pairs = list(pairwise(segments, retro=True))
        # Place the last element as the first
        segment_pairs = [ segment_pairs[-1] ] + segment_pairs[0:-1]
        corners = []
        for p, point in enumerate(points):
            corner = Corner(*point.coords, *segment_pairs[p])
            corners.append(corner)
            #point.segments = segment_pairs[p]
        return corners

    # Get horizontal size
    def get_x_size(self) -> number:
        return self.pmax.x - self.pmin.x

    # Get vertical size
    def get_y_size(self) -> number:
        return self.pmax.y - self.pmin.y

    # Get both dimension sizes
    def get_size(self) -> tuple:
        x_size = self.get_x_size()
        y_size = self.get_y_size()
        return x_size, y_size

    # Calculate the rectangle center
    def get_center(self) -> tuple:
        x_size, y_size = self.get_size()
        x_position = self.pmin.x + x_size / 2
        y_position = self.pmin.y + y_size / 2
        return x_position, y_position

    # Calculate the rectangle area
    def get_area(self) -> number:
        if self._area:
            return self._area
        x_size, y_size = self.get_size()
        self._area = resolute(x_size * y_size, base_resolution - 2)
        return self._area

    # The area is treated appart since it may be an expensive calculation
    area = property(get_area, None, None, "The rectangle area")

    # Return a segment which crosses the rectangle in diagonal from min to max point
    def get_crossing_segment(self) -> Segment:
        return Segment(self.pmin, self.pmax, self.segments_color)

    # Get the point in the middle of the rect
    def get_middle_point (self) -> Point:
        # Get the point in the middle of the crossing line
        crossing_segment = self.get_crossing_segment()
        return crossing_segment.get_middle_point()

    # Create a new rectangle identical to this rectangle but with a specified color
    def get_colored_rect(self, segments_color : str = 'black', fill_color : str = 'white') -> 'Rect':
        return Rect(self.pmin, self.pmax, segments_color, fill_color)

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
    def overlap_rect (self, rect) -> Optional['Rect']:
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
    def join_rect (self, rect) -> list:
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

# A perimeter defined by several segments
# The perimeter parameters must never be modified
# Instead, a new perimeter must be created every time a perimeter needs to be modified
class Perimeter:

    def __init__(self, segments : list, segments_color = 'black', fill_color = 'white'):
        self.segments = segments
        self.segments_color = segments_color
        self.fill_color = fill_color
        # Color at this moment all segments
        for segment in segments:
            segment.color = segments_color
        # Check the perimeter is closed and segments are not diagonal
        self.check()
        print(self.segments)
        self._corners = None
        self._rects = None
        self._mrects = None
        self._area = None

    # Set the perimeter from its corners
    # The new rect will contain all rects
    @classmethod
    def from_corners(cls, corners : List[Corner]):
        # Check that there are at least 4 points
        if len(corners) < 4:
            raise RuntimeError('It is required at least 4 points')
        # Set the segment between each pair of points
        segments = []
        for a,b in pairwise(corners, retro=True):
            new_segment = Segment(a,b)
            segments.append(new_segment)
        # Build the rectangle
        return cls(segments)

    def __str__(self):
        return ', '.join([str(segment) for segment in self.segments])

    def __contains__(self, other):
        if isinstance(other, Point):
            for rect in self.rects:
                if other in rect:
                    return True
            return False
        if isinstance(other, Segment):
            in_a = other.a in self
            cross_any_segment = any([ other.get_intersection_point(segment, in_extremis=0) for segment in self.segments ])
            return in_a and not cross_any_segment
        if isinstance(other, Rect):
            return all([ segment in self for segment in other.get_segments() ])
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
        # Apply the perimter color
        for rect in rects:
            rect.color = self.fill_color
        self._rects = rects
        return rects

    # The area is treated appart since it may be an expensive calculation
    rects = property(get_rects, None, None, "The rectangles which form the perimeter")

    # Get all rectangles which form the perimeter
    def get_mrects(self):
        # If rects are previously calculated then return them
        if self._mrects:
            return self._mrects
        # Get all rectangles which form the perimeter
        mrects = get_maximum_rectangles(self.rects)
        self._mrects = mrects
        return mrects

    # The area is treated appart since it may be an expensive calculation
    mrects = property(get_mrects, None, None, "The maximum rectangles posible inside the perimeter")

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
        self._area = resolute(sum(areas), base_resolution - 2)
        return self._area

    # The area is treated appart since it may be an expensive calculation
    area = property(get_area, None, None, "The area inside the perimeter")

    # Check if the current segments create a closed perimeter or if it is open
    def is_closed(self):
        # Check that each segment ends in the same point that the next segment starts
        for current, nextone in pairwise(self.segments, retro=True):
            if current.b != nextone.a:
                return False
        return True

    # Format segments in a way that segments and their points are ordered
    # i.e. each segment 'b' point is the 'a' point of the next segment
    # The first segment 'a' point must be the final segment 'b' point
    # If the list of segments cannot be formatted like this then return an error
    open_perimeter_error = RuntimeError('The perimeter is not closed')
    def close_segments(self):
        # Check that each segment ends in the same point that the next segment starts
        ordered_segments = [ self.segments[0] ]
        available_segments = self.segments[1:]
        while ( len(ordered_segments) < len(self.segments) ):
            # If there are no more available segments return error
            # This may happen in case there was a duplicated segment, which means the perimeter is wrong
            if len(available_segments) == 0:
                raise self.open_perimeter_error
            # Get the last ordered segment to find which is the next connected segment
            last_ordered_segment = ordered_segments[-1]
            last_point = last_ordered_segment.b
            # Get the segment which is connected with the previous segment
            connected_segment = next((segment for segment in available_segments if last_point in segment.points), None)
            if not connected_segment:
                raise self.open_perimeter_error
            # Remove the connected segment from the available 
            available_segments = [ segment for segment in available_segments if segment != connected_segment ]
            # The connected segment must be connected by the 'a' point
            # If it is connected by the 'b' point then get the inverted segment
            if connected_segment.b == last_point:
                connected_segment = connected_segment.inverted()
            ordered_segments.append(connected_segment)
        # Finally, check that the first segment and the last segment are also connected
        if ordered_segments[0].a != ordered_segments[-1].b:
            raise self.open_perimeter_error
        self.segments = ordered_segments

    # Merge continuous segments (i.e. connected segments in the same segment)
    def merge_segments (self):
        count = 0
        while (count < len(self.segments)):
            for current, nextone in pairwise(self.segments, retro=True):
                if current.same_line_as(nextone):
                    merged_segment = Segment(current.a, nextone.b)
                    self.segments = [ seg for seg in self.segments if seg not in [current, nextone] ]
                    self.segments.insert(count, merged_segment)
                    break
                count += 1

    # Check the perimeter to be closed
    # Try to close it if possible. Otherwise return an error, since not closed perimeters are not supported
    # DANI: En principio los perímetros no cerrados no tendrán soporte nunca porque no tienen mucho sentido o no les veo la utilidad
    # Check each perimeter segment to not be diagonal
    # Otherwise return an error, since perimeters with diagonal segments are not yet supported
    # DANI: En principio algun día se podría hacer esto
    # Merge continuous segments (i.e. connected segments in the same line)
    def check(self):
        # Check the perimeter to be closed
        if not self.is_closed():
            print('REMAKE!!')
            self.close_segments()
        # Check each perimeter segment to not be diagonal
        for segment in self.segments:
            if segment.is_diagonal():
                raise RuntimeError('The perimeter has diagonal segments, which are not supported')
        # Merge continuous segments (i.e. connected segments in the same line)
        self.merge_segments() 

    # Set the perimeter corners as points with additional stored values
    # The 'segments' variable defines the two segments of the perimeter which form the corner itself
    # The 'inside' variable defines if the corner is "pointing" to the inside of the perimeter
    # (e.g. rectangles have no inside corners while an 'L' has one inside corner)
    def set_corners(self):
        # For each segment, save the end point and if the direction of the next segment is left
        # Count how many corners are there in one direction (left in this case)
        precorners = []
        left_count = 0
        for current, nextone in pairwise(self.segments, retro=True):
            point = current.b
            # Given two continue segments which inevitably form a corner, set if the second segment goes "left" (true) or "right" (false) relative to the first segment
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
            raise RuntimeError('There is something wrong with the perimeter')

        # Check if are more corners in the counted direction (true) or the other (false)
        lefted_perimeter = difference > 0

        # Now set the real corners and save them in the global variable 'corners'
        # Corners with the minoritarian direction will be set as 'inside = True'
        corners = []
        for precorner in precorners:
            point, segments, lefted_corner = precorner
            inside = lefted_corner != lefted_perimeter
            corner = Corner(*point.coords, *segments, inside)
            corners.append(corner)
        return corners

    # Get all perimeter points
    # DANI: Obsoleto
    def get_points(self):
        return [ segment.a for segment in self.segments ]

    # Get a rectangle which contains the whole perimeter
    def get_box(self):
        points = self.get_points()
        x_coords = [ point.x for point in points ]
        y_coords = [ point.y for point in points ]
        x_min = min(x_coords)
        x_max = max(x_coords)
        y_min = min(y_coords)
        y_max = max(y_coords)
        pmin = Point(x_min, y_min)
        pmax = Point(x_max, y_max)
        return Rect(pmin, pmax)

    # Find out if a point is in the border of the perimeter (segments and corners)
    def in_border(self, point : Point) -> bool:
        for corner in self.corners:
            if point == corner:
                return True
        for segment in self.segments:
            if point in segment:
                return True
        return False

    # Get a specific perimeter corner or segment by specifying a point that matches this element
    def get_border_element (self, point : Point):
        for corner in self.corners:
            if point == corner:
                return corner
        for segment in self.segments:
            if point in segment:
                return segment
        return None

    # Return all points in the perimeter segments which intersect with a given segment
    # Points are sorted according to their distance with the 'a' point of the segment (from less to more distance)
    # WARNING: Intersection of paralel (overlapping) segments is not detected
    # WARNING: If the origen point is in a corner/segment it may return or not the origen as intersection point
    # DANI: No se ha probado a fondo
    # DANI: Actualmente NO está en uso
    def get_segment_intersection_points(self, segment : Segment) -> Optional[list]:
        # Get the intersection point of the specfied segment with each perimeter limit
        intersection_points = []
        for limit in self.segments:
            point = limit.get_intersection_point(segment)
            if point:
                intersection_points.append(point)
        # Find out also if the segments intersects any corner
        for corner in self.corners:
            if corner in segment:
                intersection_points.append(corner)
        # If no points are found return None
        if len(intersection_points) == 0:
            return None
        # Sort the points
        def by_distance(point):
            return segment.a.get_distance_to(point)
        sorted_points = sorted(intersection_points, key=by_distance)
        return sorted_points

    # Return all segments inside the perimeter which intersect with a given segment
    # segments are sorted according to their distance with the 'a' point of the segment (from less to more distance)
    # WARNING: Intersection of paralel (overlapping) segments is not detected
    # DANI: No se ha probado
    # DANI: No está acabado. No se contempla la posibilidad de que la segmenta que intersecta empieze dentro del perímetro
    # DANI: Actualmente NO está en uso
    def get_segment_intersection_segments(self, segment : Segment) -> Optional[list]:
        # Get the intersection points
        sorted_points = self.get_segment_intersection_points(segment)
        # If no points are found return None
        points_count = len(sorted_points)
        if points_count == 0:
            return None
        # Create new segments with the intersecting points
        intersecting_segments = []
        for a, b in pairwise(sorted_points):
            intersecting_segments.append(Segment(a,b))
        # If the last intersecting point has no pair then use the segment 'b' point a the end of the last intersecting segment
        if points_count % 2 == 1:
            last_intersection_point = sorted_points[-1]
            last_segment_point = segment.b
            last_intersection_segment = Segment(last_intersection_point, last_segment_point)
            intersecting_segments.append(last_intersection_segment)
        return intersecting_segments

    # Given a perimeter corner which is pointing inside,
    # Create two segments opposed to the corner segments but as long as required to cut the perimeter at onther segment or corner
    # This two segments will always be inside the perimeter
    def get_corner_insider_segments(self, corner : Corner, limit_points : list = [], limit_segments : list = []) -> list:

        # Get the length of the most large possible segment in the perimeter
        max_length = self.get_box().get_crossing_segment().length

        # Get the oppoiste segments to the corner segments but as long as the max_length
        # NEVER FORGET: The corner point is the entry segment 'b' point and the exit segment 'a' point
        entry_segment, exit_segment = corner.segments
        tracing1 = Segment(corner, corner + entry_segment.vector.normalized() * max_length, color='green')
        tracing2 = Segment(corner, corner + (-exit_segment.vector.normalized()) * max_length, color='green')

        # Set the limits if they were not passed
        if not limit_points or len(limit_points) == 0:
            limit_points = self.corners
        if not limit_segments or len(limit_segments) == 0:
            limit_segments = self.segments

        insider_segments = []
        for segment in [tracing1, tracing2]:

            # Get the intersection point of the specfied segment with each perimeter limit
            intersection_points = []
            for limit_segment in limit_segments:
                # Skip the segments of the main corner
                if limit_segment in corner.segments:
                    continue
                point = limit_segment.get_intersection_point(segment)
                if point:
                    # Intersection point may be the corner even ignoring corner segments
                    # This happens when one of the corner segments has been splited
                    if point == corner:
                        continue
                    intersection_points.append(point)
            # Find out also if the segment intersects any corner
            for limit_point in limit_points:
                # Skip the main corner
                if limit_point == corner:
                    continue
                if limit_point in segment:
                    intersection_points.append(limit_point)

            # There should be always at least 1 intersection
            if len(intersection_points) == 0:
                raise RuntimeError('An inside segment has no intersection point: ' + str(segment))

            # Sort the points by distance
            def by_distance(point):
                return segment.a.get_distance_to(point)
            sorted_points = sorted(intersection_points, key=by_distance)

            # The closest point will be the first point
            cut_point = sorted_points[0] # DANI: Me quedé aquí. Parece que una segmenta no corta, hay que chequear

            insider_segment = Segment(segment.a, cut_point, color='red')
            insider_segments.append(insider_segment)

        #return [tracing1, tracing2]
        return insider_segments

    # Split the current perimeter in a list of rectangles
    # The exclusion perimeters are those perimeters inside the current perimeter which are not considered
    # (e.g. children perimeters)
    def split_in_rectangles(self, exclusion_perimeters : list = []) -> list:

        # Inside corners are the parent perimeter inside corners and the children perimeters outside corners
        parent_inside_corners = [ corner for corner in self.corners if corner.inside == True ]
        children_inside_corners = []
        for perimeter in exclusion_perimeters:
            child_inside_corners = [ corner for corner in perimeter.corners if corner.inside == False ]
            children_inside_corners += child_inside_corners
        inside_corners = parent_inside_corners + children_inside_corners

        # Limits are the current perimeter and all exclusion perimeters
        limits = [ self, *exclusion_perimeters ]

        # Limit points are 'inside corners' which are not overlapped with any other perimeters
        limit_points = []
        for corner in inside_corners:
            count = 0
            for perimeter in limits:
                element = perimeter.get_border_element(corner)
                if element:
                    count += 1
            # Include this corner only if after all searches it was found only 1 time
            if count == 1:
                limit_points.append(corner)

        # Limit segments are both parent and children segments
        limit_segments = []
        for perimeter in limits:
            limit_segments += perimeter.segments

        # WARNING: If there are no inside corners at this point it means our perimeter/s are made only by single rectangles
        # However those rectangles may be formed by segments longer than the rectangle size, so they may require splitting

        # Get all inside separator segments
        inside_separators = []
        for corner in limit_points:
            for segment in self.get_corner_insider_segments(corner, limit_points, limit_segments):
                inside_separators.append(segment)

        # Remove duplicates
        inside_separators = list(set(inside_separators))

        # Join all inside segments with the perimeter limits in a single list
        # WARNING: inside segments MUST be before limit segments
        all_segments = [ *inside_separators, *limit_segments ]

        # DANI: Usa esto para ver los perímetros en negro y las segmentas interiores en rojo
        #add_frame(all_segments)

        # Find all points where the inside separator segments intersect each other
        all_intersections = []
        for segment1, segment2 in itertools.combinations(all_segments, 2):
            # All inside corners would be found as intersection points, since their 2 segments intersect
            # For this reason, we set 'in_extremis' as false
            intersection = segment1.get_intersection_point(segment2)
            if not intersection:
                continue
            all_intersections.append(intersection)

        # Remove duplicates
        all_intersections = list(set(all_intersections))

        # Split the inside separator segments at the intersection points
        all_splitted_segments = []
        for segment in all_segments:
            splitted_segments = segment.split_at_points(all_intersections)
            all_splitted_segments += list(splitted_segments)
        #print('All segments: ' + str(len(all_splitted_segments)))

        # Finally, for each segment, try to find 2 rectangles
        # Each segment will be connected to exactly 1 or 2 rectangles
        final_rectangles = []
        for segment in all_splitted_segments:
            this_segment_rects = [] # Max 2
            count = 0
            for other_segment in all_splitted_segments:
                # If 2 segments are connected making a corner there may be a rectangle
                if segment.makes_corner_with(other_segment):
                    # Create the rect that the two previous segments would make
                    new_rect = Rect.from_segments([ segment, other_segment ])
                    #print('New rect: ' + str(segment) + ' / ' + str(other_segment) + ' -> ' + str(count) + ': ' + str(new_rect))
                    # Must check that the other 2 segments which would close (complete) the rectangle do exist
                    # Actually they may be splitted is several segments, so we must check that they are covered
                    closing_segments = [ rect_segment for rect_segment in new_rect.segments if rect_segment not in [ segment, other_segment ] ]
                    if not all([ segment.is_covered_by(all_splitted_segments) for segment in closing_segments ]):
                        continue
                    if new_rect not in this_segment_rects:
                        count += 1
                        this_segment_rects.append(new_rect)
                    if count == 2:
                        break
            #print(str(segment) + ' -> ' + str(this_segment_rects[0]) + ' / ' + str(this_segment_rects[1]))
            final_rectangles += this_segment_rects

        # Remove duplicates
        final_rectangles = list(set(final_rectangles))

        # In some cases, some rectangles may be found inside exclusion perimeters
        # Find and discard those rectangles
        for perimeter in exclusion_perimeters:
            final_rectangles = [ rect for rect in final_rectangles if rect not in perimeter ]

        #print('Total rectangles: ' + str(len(final_rectangles)))

        return final_rectangles

    # Check if a rectangle fits somewhere in the perimeter
    # Iterate over all maximum rectangles and try to fit 
    def fit (self, x_fit_size : number, y_fit_size : number):
        for rect in self.mrects:
            x_size, y_size = rect.get_size()
            if x_fit_size <= x_size and y_fit_size <= y_size:
                return True
        return False
                

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

# One by one for each *available rectangle, where available rectangles are the splitted rectangles
# Get as many rectanges as possible which are connected horizontally to the current rectangle
# Get as many rows of rectangles as possible which are connected vertically to all previous rectangles
# Consider all previous rectangles as a single rectange
# Repeat in the inverse order (first vertically, then horizontally)
# Remove all previous rectangles from the *available rectangles list
def get_maximum_rectangles(splitted_rects : list = []) -> list:

    # First, some functions are defined to find colliding rects
    def get_left_rect (rect : Rect):
        pmax = rect.get_upper_left_point()
        for r in splitted_rects:
            if r.pmax == pmax:
                return r
        return None
    def get_right_rect (rect : Rect):
        pmin = rect.get_bottom_right_point()
        for r in splitted_rects:
            if r.pmin == pmin:
                return r
        return None
    def get_upper_rect (rect : Rect):
        pmin = rect.get_upper_left_point()
        for r in splitted_rects:
            if r.pmin == pmin:
                return r
        return None
    def get_bottom_rect (rect : Rect):
        pmax = rect.get_bottom_right_point()
        for r in splitted_rects:
            if r.pmax == pmax:
                return r
        return None

    # Set a function to merge multiple rects into a single big rect
    # Do it by finding the most maximum pmax and the most minimum pmin
    def sort_by_x(point):
        return point.x
    def sort_by_y(point):
        return point.y
    def merge_rectangles(rects):
        pmax_points = [ rect.pmax for rect in rects ]
        sorted_pmax_points = sorted( sorted( pmax_points, key=sort_by_x, reverse=True ), key=sort_by_y, reverse=True )
        maximum_pmax = sorted_pmax_points[0]
        pmin_points = [ rect.pmin for rect in rects ]
        sorted_pmin_points = sorted( sorted( pmin_points, key=sort_by_x), key=sort_by_y )
        minimum_pmin = sorted_pmin_points[0]
        return Rect(minimum_pmin, maximum_pmax)

    # Trace which rectangles have been already checked both horizontally and vertically to avoid repeating
    for splitted_rect in splitted_rects:
        splitted_rect.horizontal_check = False
        splitted_rect.vertical_check = False

    # Save all maximum rectangles in a list to be returned at the end
    maximum_rectangles = []

    # Group rectangles first horizontally and then vertically
    for splitted_rect in splitted_rects:
        # Skip already grouped rectangles horizontally
        if splitted_rect.horizontal_check == True:
            continue
        # Set the first row
        first_row = [ splitted_rect ]
        # Append all rectangles at right from current rectangle to the row
        rightest = splitted_rect
        while (True):
            rightest = get_right_rect(rightest)
            if not rightest:
                break
            first_row.append(rightest)
        # Append all rectangles at left from current rectangle to the row
        leftest = splitted_rect
        while (True):
            leftest = get_left_rect(leftest)
            if not leftest:
                break
            first_row.append(leftest)
        # Set all rectangles in the first row as checked horizontally
        for rect in first_row:
            rect.horizontal_check = True
        # Set the group of rectangles to be joined
        group = [ rect for rect in first_row ]
        # If all rectangles in the row have a botton rectangle then add all those new rects to a new row
        # This new row is then added to the whole group of rectanges and used to find the next row
        current_row = [ rect for rect in first_row ]
        while (True):
            new_row = [ get_bottom_rect(rect) for rect in current_row ]
            if not all(new_row):
                break
            group += new_row
            current_row = [ *new_row ]
        # Repeat the process but upperwards
        current_row = [ rect for rect in first_row ]
        while (True):
            new_row = [ get_upper_rect(rect) for rect in current_row ]
            if not all(new_row):
                break
            group += new_row
            current_row = [ *new_row ]
        # Create a new rect which contains all group rects
        maximum_rect = merge_rectangles(group)
        # Add the new maximum rectnagle to the list if it is not there already
        if maximum_rect not in maximum_rectangles:
            maximum_rectangles.append(maximum_rect)

    # Now repeat the process in the inverse order (first vertically, then horizontally)
    for splitted_rect in splitted_rects:
        # Skip already grouped rectangles vertically
        if splitted_rect.vertical_check == True:
            continue
        # Set the first column
        first_column = [ splitted_rect ]
        # Append all rectangles at the bottom from current rectangle to the column
        bottomest = splitted_rect
        while (True):
            bottomest = get_bottom_rect(bottomest)
            if not bottomest:
                break
            first_column.append(bottomest)
        # Append all rectangles upper from current rectangle to the column
        upperest = splitted_rect
        while (True):
            upperest = get_upper_rect(upperest)
            if not upperest:
                break
            first_column.append(upperest)
        # Set all rectangles in the first column as checked vertically
        for rect in first_column:
            rect.vertical_check = True
        # Set the group of rectangles to be joined
        group = [ rect for rect in first_column ]
        # If all rectangles in the column have a right rectangle then add all those new rects to a new column
        # This new column is then added to the whole group of rectanges and used to find the next column
        current_column = [ rect for rect in first_column ]
        while (True):
            new_column = [ get_right_rect(rect) for rect in current_column ]
            if not all(new_column):
                break
            group += new_column
            current_column = [ *new_column ]
        # Repeat the process but to the left
        current_column = [ rect for rect in first_column ]
        while (True):
            new_column = [ get_left_rect(rect) for rect in current_column ]
            if not all(new_column):
                break
            group += new_column
            current_column = [ *new_column ]
        # Create a new rect which contains all group rects
        maximum_rect = merge_rectangles(group)
        # Add the new maximum rectnagle to the list if it is not there already
        if maximum_rect not in maximum_rectangles:
            maximum_rectangles.append(maximum_rect)

    return maximum_rectangles

# One by one for each *available rectangle, where available rectangles are the splitted rectangles
# Get as many rectanges as possible which are connected horizontally to the current rectangle
# Consider all previous rectangles as a single rectange
# Remove all previous rectangles from the *available rectangles list
def get_row_rectangles(splitted_rects : list = []) -> list:

    # First, some functions are defined to find colliding rects
    def get_left_rect (rect : Rect):
        pmax = rect.get_upper_left_point()
        for r in splitted_rects:
            if r.pmax == pmax:
                return r
        return None
    def get_right_rect (rect : Rect):
        pmin = rect.get_bottom_right_point()
        for r in splitted_rects:
            if r.pmin == pmin:
                return r
        return None

    # Set a function to merge multiple rects into a single big rect
    # Do it by finding the most maximum pmax and the most minimum pmin
    def sort_by_x(point):
        return point.x
    def merge_rectangles(rects):
        pmax_points = [ rect.pmax for rect in rects ]
        sorted_pmax_points = sorted( pmax_points, key=sort_by_x, reverse=True )
        maximum_pmax = sorted_pmax_points[0]
        pmin_points = [ rect.pmin for rect in rects ]
        sorted_pmin_points = sorted( pmin_points, key=sort_by_x )
        minimum_pmin = sorted_pmin_points[0]
        return Rect(minimum_pmin, maximum_pmax)

    # Trace which rectangles have been already checked both horizontally and vertically to avoid repeating
    for splitted_rect in splitted_rects:
        splitted_rect.check = False

    # Save all maximum rectangles in a list to be returned at the end
    rows = []

    # Group rectangles first horizontally and then vertically
    for splitted_rect in splitted_rects:
        # Skip already grouped rectangles horizontally
        if splitted_rect.check == True:
            continue
        # Start the row
        row = [ splitted_rect ]
        # Append all rectangles at right from current rectangle to the row
        rightest = splitted_rect
        while (True):
            rightest = get_right_rect(rightest)
            if not rightest:
                break
            row.append(rightest)
        # Append all rectangles at left from current rectangle to the row
        leftest = splitted_rect
        while (True):
            leftest = get_left_rect(leftest)
            if not leftest:
                break
            row.append(leftest)
        # Set all rectangles in the row as checked
        for rect in row:
            rect.check = True
        # Merge all rectangles in row
        row_rectangle = merge_rectangles(rows)
        rows.append(row_rectangle)

    return rows

# One by one for each *available rectangle, where available rectangles are the splitted rectangles
# Get as many rectanges as possible which are connected vertically to the current rectangle
# Consider all previous rectangles as a single rectange
# Remove all previous rectangles from the *available rectangles list
def get_column_rectangles(splitted_rects : list = []) -> list:

    # First, some functions are defined to find colliding rects
    def get_upper_rect (rect : Rect):
        pmin = rect.get_upper_left_point()
        for r in splitted_rects:
            if r.pmin == pmin:
                return r
        return None
    def get_bottom_rect (rect : Rect):
        pmax = rect.get_bottom_right_point()
        for r in splitted_rects:
            if r.pmax == pmax:
                return r
        return None

    # Set a function to merge multiple rects into a single big rect
    # Do it by finding the most maximum pmax and the most minimum pmin
    def sort_by_y(point):
        return point.y
    def merge_rectangles(rects):
        pmax_points = [ rect.pmax for rect in rects ]
        sorted_pmax_points = sorted( pmax_points, key=sort_by_y, reverse=True )
        maximum_pmax = sorted_pmax_points[0]
        pmin_points = [ rect.pmin for rect in rects ]
        sorted_pmin_points = sorted( pmin_points, key=sort_by_y )
        minimum_pmin = sorted_pmin_points[0]
        return Rect(minimum_pmin, maximum_pmax)

    # Trace which rectangles have been already checked both horizontally and vertically to avoid repeating
    for splitted_rect in splitted_rects:
        splitted_rect.check = False

    # Save all maximum rectangles in a list to be returned at the end
    columns = []

    # Group rectangles first horizontally and then vertically
    for splitted_rect in splitted_rects:
        # Skip already grouped rectangles horizontally
        if splitted_rect.check == True:
            continue
        # Start the column
        column = [ splitted_rect ]
        # Append all rectangles at the bottom from current rectangle to the column
        bottomest = splitted_rect
        while (True):
            bottomest = get_bottom_rect(bottomest)
            if not bottomest:
                break
            column.append(bottomest)
        # Append all rectangles upper from current rectangle to the column
        upperest = splitted_rect
        while (True):
            upperest = get_upper_rect(upperest)
            if not upperest:
                break
            column.append(upperest)
        # Set all rectangles in the column as checked
        for rect in column:
            rect.check = True
        # Merge all rectangles in column
        column_rectangle = merge_rectangles(column)
        columns.append(column_rectangle)

    return columns