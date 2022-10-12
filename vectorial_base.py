from typing import Union, Optional, List, Tuple, Generator

from scheme_display import add_frame

from functools import reduce

import itertools

from math import sqrt, isclose, sin, cos, degrees, radians, acos

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
# When working with substractions (e.g. distances) the resolution should be decreased in 1 order
# When working with squared substractions (e.g. areas) the resolution should be decreased in 2 orders
base_resolution = 4
def resolute(num, base_offset = 0):
    base = base_resolution + base_offset
    resolution = 10**base
    return round(num * resolution) / resolution

# Calculate the minimum resolution
# i.e. changes below this resolution make no sense
minimum_resolution = 1 / 10 ** base_resolution
# Set functions to compare float according to the minimum resolution
def equal (a : float, b : float) -> bool:
    return a < b + minimum_resolution and a > b - minimum_resolution
def greater (a : float, b : float) -> bool:
    return a > b + minimum_resolution
def lower (a : float, b : float) -> bool:
    return a < b - minimum_resolution

# Set another function to check if two numbers are very close, out of resolution matters
# We do the double check because the 'isclose' function may fail for values close to 0
# see https://stackoverflow.com/questions/35324893/using-math-isclose-function-with-values-close-to-0
def same_number(a : number, b : number) -> bool:
    # DANI: Esto debería funcionar siempre pero son 4 operaciones
    #return isclose(a,b) or isclose(a+1,b+1)
    # Esto debería funcionar bien
    return isclose(a,b, abs_tol=minimum_resolution)

# An x,y coordinate
class Point:

    def __init__(self, x : number, y : number):
        # Save the coordinates in resoluted format applying the precision limit
        self.x = resolute(x)
        self.y = resolute(y)
        # Save both coordinates as a tuple
        self.coords = self.x, self.y

    def __str__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __repr__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

    def __eq__(self, other):
        if isinstance(other, self.__class__) or issubclass(self.__class__, other.__class__):
            return self.x == other.x and self.y == other.y
        return False

    def __hash__(self):
        return hash((self.x, self.y))
    
    # Point + Vector -> Point, Point + Point -> Vector
    def __add__(self, other):
        if isinstance(other, self.__class__) or issubclass(self.__class__, other.__class__):
            return Vector(other.x - self.x, other.y - self.y)
        if isinstance(other, Vector):
            return Point(self.x + other.x, self.y + other.y)
        raise ValueError('Point addition of ' + str(other.__class__) + ' is not supported')

    # Point - Vector -> Point, Point - Point -> Vector
    def __sub__(self, other):
        if isinstance(other, self.__class__) or issubclass(self.__class__, other.__class__):
            return Vector(self.x - other.x, self.y - other.y)
        if isinstance(other, Vector):
            return Point(self.x - other.x, self.y - other.y)
        raise ValueError('Point substraction of ' + str(other.__class__) + ' is not supported')

    # Get the distance from this point to other specified point
    def get_distance_to (self, other : 'Point') -> number:
        x_distance = self.x - other.x
        y_distance = self.y - other.y
        return sqrt( x_distance**2 + y_distance**2 )

class Vector:

    def __init__(self, x : number, y : number):
        # Save the coordinates as the whole float
        self.x = x
        self.y = y

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
            return same_number(self.x, other.x) and same_number(self.y, other.y)
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
        if same_number(self.x, 0):
            return None
        return self.y / self.x

    # Return a new vector with identical direction and sense but magnitude = 1
    def normalized(self) -> 'Vector':
        magnitude = self.get_magnitude()
        if magnitude == 0:
            raise ValueError('Can not normalize Vector(0,0)')
        return self / magnitude

    # Find out if the vector is totally vertical
    def is_vertical(self) -> bool:
        return same_number(self.x, 0)
    # Find out if the segment is totally horizontal
    def is_horizontal(self) -> bool:
        return same_number(self.y, 0)

    # Find out if the segment is diagonal
    def is_diagonal(self) -> bool:
        return not self.is_vertical() and not self.is_horizontal()

    # Find out if two vectors are equivalent
    # i.e. they have the same direction and magnitude, no matter the sense
    def is_equivalent_to (self, other : 'Vector') -> bool:
        same_slope = same_number(self.get_slope(), other.get_slope())
        same_magnitude = same_number(self.get_magnitude(), other.get_magnitude())
        return same_slope and same_magnitude

    # Get the angle between this vector and other vector
    def get_angle_with (self, other : 'Vector') -> number:
        # By normalizing vectors we skip having to divide the dot product by the vector magnitudes product
        a = self.normalized()
        b = other.normalized()
        dot_product = a.x * b.x + a.y * b.y
        angle = degrees(acos(dot_product))
        # If the second segment is pointing counterclockwise respect to the first vector then turn the angle negative
        cz = a.x * b.y - a.y * b.x
        if cz > 0:
            angle = -angle
        return angle

    # Given a vector and an angle get a new clockwise rotated vector
    # Source: https://en.wikipedia.org/wiki/Rotation_matrix
    def rotate (self, angle) -> 'Vector':
        angle_sin = sin(radians(angle))
        angle_cos = cos(radians(angle))
        x = + angle_cos * self.x + angle_sin * self.y
        y = - angle_sin * self.x + angle_cos * self.y
        #print(str(self) + ' (' + str(angle) + ') ' + str(Vector(x,y)))
        return Vector(x,y)

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

    def __hash__(self):
        return hash(self.get_hash())

    def __contains__(self, other):
        if isinstance(other, Point):
            return self.line_contains_point(other)
        if isinstance(other, Segment):
            return self.same_line_as(other)
        return False

    # Get a value which is always the same for equal lines, no matter the specific point and vector
    def get_hash(self) -> tuple:
        slope = self.get_slope()
        intercept = self.get_intercept_point()
        return slope, intercept

    # Get the vector from the line
    # Return None in case the line is vertical
    def get_slope(self) -> Optional[number]:
        return self.vector.get_slope()

    # Get the line y intercept
    # i.e. the y coordinate from the point in the vertical line x = 0 where this line cuts
    # Return None in case the line is vertical
    def get_y_intercept (self) -> number:
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
    def get_x_intercept(self) -> number:
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

    # Check if a point is inside the line
    # WARNING: This function is independent from __contain__ since it must be inherited by the Segment class
    def line_contains_point (self, point : Point) -> bool:
        if self.is_vertical():
            x = self.get_x_intercept()
            return resolute(x) == point.x
        slope = self.get_slope()
        intercept = self.get_y_intercept()
        y = point.x * slope + intercept
        return resolute(y) == point.y

    # Get the intersection point between two lines
    # DANI: Esto está hecho en la libreta
    def get_line_intersection_point (self, line : 'Line') -> Optional[Point]:
        if self.is_paralel_to(line):
            return None
        # Asuming 2 vector equation of the line, which can be equaled since they are intersecting, we have isolated one of their "alphas"
        # In this case it is the alpha of 'line', not 'self'
        a1 = self.point.x
        b1 = self.point.y
        x1 = self.vector.x
        y1 = self.vector.y
        a2 = line.point.x
        b2 = line.point.y
        x2 = line.vector.x
        y2 = line.vector.y
        alpha = ( a2*y1 - a1*y1 - b2*x1 + b1*x1 ) / ( y2*x1 - x2*y1 )
        # Now solve the corresponding point in the line
        x = a2 + x2 * alpha
        y = b2 + y2 * alpha
        return Point(x, y)

    # Get the intersection point between this line and a segment
    def get_segment_intersection_point (self, segment : 'Segment') -> Optional[Point]:
        line_intersection_point = self.get_line_intersection_point(segment.line)
        if line_intersection_point in segment:
            return line_intersection_point
        return None

    # Get the angle between this line and other line
    def get_angle_with (self, other : 'Line') -> number:
        return self.vector.get_angle_with(other.vector)

    # Get the closer point in the line to other specified point 
    def get_closer_point (self, other : 'Point') -> 'Point':
        perpendicular_vector = self.vector.rotate(90)
        perpendicular_line = Line(other, perpendicular_vector)
        return self.get_line_intersection_point(perpendicular_line)

    # Get the distance from this line to other specified point (i.e. from the closer point in the line)
    def get_distance_to (self, other : 'Point') -> number:
        closer_point = self.get_closer_point(other)
        return closer_point.get_distance_to(other)
        
# A segment defined by 2 coordinates (Points): 'a' and 'b'
class Segment(Line):

    def __init__(self, a : Point, b : Point, color : str = 'black'):
        if a == b:
            raise ValueError('Inexistent segment. Points "a" and "b" must be different: ' + str(a))
        self.a = a
        self.b = b
        self.color = color
        vector = a + b
        super().__init__(a, vector)
        self.line = Line(a, vector)
        self.direction = vector.normalized()
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

    def __neg__(self):
        return self.inverted

    def __contains__(self, other):
        if isinstance(other, Point):
            # First of all check that the point is in the line
            # WARNING: This may seem redundant if then we use the distances method, but it is not
            # WARNING: The distances method may return True when a point is just close to the line
            # WARNING: The more far the segment points are the less resolution this method has
            if not self.line_contains_point(other):
                return False
            # Now that we know the point is in the line, check if it is between both segment points
            distance1 = self.a.get_distance_to(other)
            distance2 = self.b.get_distance_to(other)
            return same_number(distance1 + distance2, self.length)
        if isinstance(other, self.__class__):
            return other.a in self and other.b in self
        return False

    # Get a value which is always the same for equal segments no matter the order of a and b
    def get_hash(self) -> tuple:
        points = [self.a, self.b]
        sorted_points = sort_points(points)
        return tuple(sorted_points)

    # Get the inverted segment (i.e. invert a and b points)
    def inverted (self):
        return Segment(self.b, self.a, self.color)

    # Create a new segment identical to this segment but with a specified color
    def get_colored_segment (self, color : str):
        return Segment(self.a, self.b, color)

    # Check if one of the segment points matches a point
    def has_point (self, point : Point) -> bool:
        return point == self.a or point == self.b

    # Given a segment point (a or b) return the other one
    def get_other_point (self, point : Point) -> Point:
        if point == self.a:
            return self.b
        elif point == self.b:
            return self.a
        else:
            raise ValueError('Input point is not one of the segment points')

    # Check if two segments have a point in common
    # WARNING: They may be overlapped
    def is_connected_with (self, other : 'Segment') -> bool:
        return self.has_point(other.a) or self.has_point(other.b)

    # Check if two segments form a corner
    # i.e. they have a point in common and both segments are not paralel
    def makes_corner_with (self, other : 'Segment') -> bool:
        if self.is_paralel_to(other):
            return False
        return self.is_connected_with(other)

    # Check if a segment is fully covered by a list of segments
    # i.e. substract all segments from self and if there is still a remaining segment it is not covered
    def is_covered_by (self, segments : List['Segment']) -> bool:
        remaining_segments = self.substract_segments(segments)
        return len(remaining_segments) == 0

    # Split the segment at multiple points and return the splitted segments
    # Duplicated or non-cutting points are discarded
    def split_at_points (self, points : List[Point]) -> Generator['Segment', None, None]:
        # Get only those points which are cutting the segment
        cutting_points = [ point for point in set(points) if point in self and point != self.a and point != self.b ]
        # If no points are cutting the segment then return only the intact segment
        if len(cutting_points) == 0:
            yield self
            return
        # Sort points by distance with the segment 'a' point
        def by_distance (point):
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
    # in_extremis = 2 -> All intersections are considered
    def get_intersection_point (self, segment, in_extremis : int = 2) -> Optional[Point]:
        xdiff = Vector(self.a.x - self.b.x, segment.a.x - segment.b.x)
        ydiff = Vector(self.a.y - self.b.y, segment.a.y - segment.b.y)

        def det(a, b) -> number:
            return a.x * b.y - a.y * b.x

        div = det(xdiff, ydiff)
        # segments are paralel
        if div == 0:
            return None

        self_det = det( Vector(self.a.x, self.a.y), Vector(self.b.x, self.b.y) )
        segment_det = det( Vector(segment.a.x, segment.a.y), Vector(segment.b.x, segment.b.y))
        d = Vector(self_det, segment_det)
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
    # CRITICAL WARNING: Relying on the result of this function is risky and it may lead to unaccuracies and unstability
    # CRITICAL WARNING: Do never use this function to solve grid or polygon split calculations
    # According to the resolutution of 0.0001, imagine we have a segment where self.length % 10000 is not 0
    # Then the middle point of a -> b will be different from the middle point of b -> a although they should be the same
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
        points = unique(self_points + other_points)
        sorted_points = sort_points(points)
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
    def substract_segments (self, others : List['Segment']) -> List['Segment']:
        # Filter others to be in the same line that self segment
        inline_segments = [ segment for segment in others if self.same_line_as(segment) ]
        if len(inline_segments) == 0:
            return [self]
        # Order segment points and check how they alternate
        self_points = [self.a, self.b]
        other_points = []
        for segment in inline_segments:
            other_points += [segment.a, segment.b]
        points = unique(self_points + other_points)
        sorted_points = sort_points(points)
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

    # Given two segments which have a point in common, return a segment from the non-common points
    def combine_segment (self, other : 'Segment') -> 'Segment':
        # Find the non-common point in self, which must be only one
        different_self_points = [ point for point in self.points if not other.has_point(point) ]
        if len(different_self_points) != 1:
            raise ValueError('Segments can not be combined')
        a = different_self_points[0]
        # Find the non-common point in other and build the new segment
        b = next( point for point in other.points if not self.has_point(point) )
        return Segment(a,b)

    # Get a point inside the current segment with a specific margin length
    def fit_point (self, margin : number) -> 'Point':
        # Check the margin width to be suitable for self segment length
        if self.length < 2 * margin:
            raise ValueError('Margins are not suitable for this segment')
        # Get a random point inside the segment where the a point of the new segment could be
        new_limit_a = self.a + self.direction * margin
        new_limit_b = self.b - self.direction * margin
        # If margins fit perfectly in the segment then the segment center is the only available point
        if new_limit_a == new_limit_b:
            return new_limit_a
        new_region = Segment(new_limit_a, new_limit_b)
        return new_region.get_random_point()

    # Get a segment inside the current segment with a specific length and margin length
    def fit_segment (self, width : number, margin : number) -> 'Segment':
        # Check this segment and the margins widths to be suitable for self segment length
        if self.length < width + 2 * margin:
            raise ValueError('Conditions are not suitable for this segment')
        # Get a random point inside the segment where the a point of the new segment could be
        new_a_limit_a = self.a + self.direction * margin
        new_a_limit_b = self.b - self.direction * (width + margin)
        new_a_region = Segment(new_a_limit_a, new_a_limit_b)
        new_a = new_a_region.get_random_point()
        new_b = new_a + self.direction * width
        return Segment(new_a, new_b)

    # Get a random point inside this segment
    def get_random_point (self) -> 'Point':
        return self.a + self.direction * (self.length * random.random())

    # Get a segment after substracting a margin at both ends
    def get_margined_segment (self, margin : number) -> 'Segment':
        # Check the margin width to be suitable for self segment length
        if self.length < 2 * margin:
            raise ValueError('Margins (' + str(margin) + ') would consume this segment totally: ' + str(self))
        if self.length == 2 * margin:
            raise ValueError('Margins (' + str(margin) + ') would consume this segment exact length: ' + str(self))
        # Build the margined segment
        new_a = self.a + self.direction * margin
        new_b = self.b - self.direction * margin
        return Segment(new_a, new_b)

    # Get a segment after translating the current segment
    def translate (self, direction : 'Vector') -> 'Segment':
        new_a = self.a + direction
        new_b = self.b + direction
        return Segment(new_a, new_b)

    # Get the closer point in the segment to other specified point 
    def get_closer_point (self, other : 'Point') -> 'Point':
        line_closer_point = self.line.get_closer_point(other)
        if line_closer_point in self:
            return line_closer_point
        point_distances = [ point.get_distance_to(other) for point in self.points ]
        if point_distances[0] >= point_distances[1]:
            return self.a
        return self.b

    # Get the distance from this segment to other specified point (i.e. from the closer point in the segment)
    # DANI: Ahora mismo no se usa
    def get_distance_to (self, other : 'Point') -> number:
        closer_point = self.get_closer_point(other)
        return closer_point.get_distance_to(other)

    # Given another segment, get the part of it which overlaps the "perpendicular range" of this segment
    # Note that they may overlap through a segment, just a point or do not overlap at all
    def get_segment_perpendicular_range_overlap (self, other : 'Segment') -> Optional[Union['Point', 'Segment']]:
        # Throw a perpendicular line from each segment point and check if they intersect with the other segment
        intersections = []
        perpendicular_vector = self.direction.rotate(90)
        for point in self.points:
            line = Line(point, perpendicular_vector)
            intersection = line.get_segment_intersection_point(other)
            if intersection:
                intersections.append(intersection)
        # If there are two intersections then it means the other segment is crossing the whole perpendicular range
        if len(intersections) == 2:
            return Segment(intersections[0], intersections[1])
        # Get the self line perpendicular point of one of the other segment points and check if it is in the perpendicular range
        other_perpendicular_point = self.line.get_closer_point(other.points[0])
        is_point_in_range = other_perpendicular_point in self
        # If there is only one intersection then it means part of the other segment is inside the range and part is outside
        if len(intersections) == 1:
            point_in_range = other.points[0] if is_point_in_range else other.points[1]
            if intersections[0] == point_in_range:
                return point_in_range
            return Segment(intersections[0], point_in_range)
        # If there are no intersections then it means the other segment is fully inside or outside the perpendicular range
        if is_point_in_range:
            return other
        return None

# A corner is a point where 2 non-paralel segments are connected
# i.e. both segments have this point ('a' or 'b') in common
# IMPORTANT: It is a standard that the first segment enters the corner while the second exits the corner
# i.e. first segment is 'a' -> (corner) and second segment is (corner) -> 'b'
# The inside argument is refered to if this corner points towards the 'inside' of the polygon which belongs
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

# A rectangular area defined by 2 coordinates (Points): 'max' and 'min'
class Rect:

    def __init__(self, x_min : number, y_min : number, x_max : number, y_max : number, segments_color : str = 'black', fill_color : str = 'white'):
        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max
        self.x = (x_min, x_max)
        self.y = (y_min, y_max)
        self.x_size, self.y_size = self.get_size()
        if resolute(self.x_size) < minimum_resolution or resolute(self.y_size) < minimum_resolution:
            raise ValueError('The rectangle has not 2 dimensions: ' + str(self))
        self.segments_color = segments_color
        self.fill_color = fill_color
        self.segments = self.get_segments()
        self._area = None

    # Set the rect from segments
    # The new rect will contain all segments
    @classmethod
    def from_segments(cls, segments : List[Segment], segments_color : str = 'black', fill_color : str = 'white'):
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
        return cls(x_min, y_min, x_max, y_max, segments_color, fill_color)

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
        return 'X: ' + str(self.x_min) + '/' + str(self.x_max) + ', Y: ' + str(self.y_min) + '/' + str(self.y_max)

    def __repr__(self):
        return 'X: ' + str(self.x_min) + '/' + str(self.x_max) + ', Y: ' + str(self.y_min) + '/' + str(self.y_max)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x_min == other.x_min and self.y_min == other.y_min and self.x_max == other.x_max and self.y_max == other.y_max
        return False

    def __hash__(self):
        return hash((self.x_min, self.y_min, self.x_max, self.y_max))

    def __contains__(self, other):
        if isinstance(other, Point):
            in_x = other.x >= self.x_min and other.x <= self.x_max
            in_y = other.y >= self.y_min and other.y <= self.y_max
            return in_x and in_y
        if isinstance(other, Segment):
            in_a = other.a in self
            in_b = other.b in self
            return in_a and in_b
        if isinstance(other, Rect):
            in_x = other.x_min >= self.x_min and other.x_max <= self.x_max
            in_y = other.y_min >= self.y_min and other.y_max <= self.y_max
            return in_x and in_y
        return False

    # Check if self rect is colliding with other rect
    # i.e. one of their segments is totally or partially overlapping
    # Note that one rectangle may be inside of the other
    def is_colliding_with (self, other : 'Rect') -> bool:
        for self_segment in self.segments:
            for other_segment in other.segments:
                overlap = self_segment.get_overlap_segment(other_segment)
                if overlap:
                    return True
        return False

    # Check if self rect is connected with other rect
    # i.e. one of their segments is totally overlapping
    # Note that rectangles cannot overlap
    def is_connected_with (self, other : 'Rect') -> bool:
        # They are aligned horizontally
        if self.x_min == other.x_min and self.x_max == other.x_max:
            # The other rect is below self rect
            if self.y_min == other.y_max:
                return True
            # The other rect is above self rect
            if self.y_max == other.y_min:
                return True
        # They are aligned vertically
        if self.y_min == other.y_min and self.y_max == other.y_max:
            # The other rect is at the left of self rect
            if self.x_min == other.x_max:
                return True
            # The other rect is at the right self rect
            if self.x_max == other.x_min:
                return True
        # They are not connected
        return False

    # Get rectangle points
    def get_bottom_left_point (self) -> Point:
        return Point(self.x_min, self.y_min)
    def get_upper_left_point (self) -> Point:
        return Point(self.x_min, self.y_max)
    def get_upper_right_point (self) -> Point:
        return Point(self.x_max, self.y_max)
    def get_bottom_right_point (self) -> Point:
        return Point(self.x_max, self.y_min)

    # Return all rectangle points in a 'polygon-friendly' order
    def get_points (self) -> List[Point]:
        point1 = self.get_bottom_left_point()
        point2 = self.get_upper_left_point()
        point3 = self.get_upper_right_point()
        point4 = self.get_bottom_right_point()
        return [ point1, point2, point3, point4 ]

    # Return all rectangle segments in a 'polygon-friendly' order
    def get_segments (self, color : Optional[str] = None) -> List[Segment]:
        points = self.get_points()
        segment_color = color if color else self.segments_color
        segments = [ Segment(a, b, segment_color) for a, b in pairwise(points, retro=True) ]
        return segments

    # Return all rectangle points in a 'polygon-friendly' order
    # Each point contains its two adjacent segments
    def get_corners (self) -> List['Corner']:
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

    # Get a corner from its point
    def get_corner (self, point : 'Point') -> Optional['Corner']:
        corners = self.get_corners()
        corner = next((corner for corner in corners if corner == point), None)
        if not corner:
            return None
        return corner

    # Get horizontal size
    def get_x_size (self) -> number:
        return self.x_max - self.x_min

    # Get vertical size
    def get_y_size (self) -> number:
        return self.y_max - self.y_min

    # Get both dimension sizes
    def get_size (self) -> tuple:
        x_size = self.get_x_size()
        y_size = self.get_y_size()
        return x_size, y_size

    # Get both dimension sizes
    def get_diagonal_size (self) -> number:
        x_size, y_size = self.get_size()
        return sqrt(x_size**2 + y_size**2)

    # Calculate the rectangle area
    def get_area (self) -> number:
        if self._area:
            return self._area
        x_size, y_size = self.get_size()
        self._area = resolute(x_size * y_size)
        return self._area

    # The rectangle area (read only)
    area = property(get_area, None, None, "The rectangle area")

    # Given a segment of the current rectangle, get the normalized direction from this segment to the rectangle center
    def get_side_direction_to_center (self, segment : Segment) -> Vector:
        segments = self.get_segments()
        # Left side
        if segment == segments[0]:
            return Vector(1,0)
        # Upper side
        if segment == segments[1]:
            return Vector(0,-1)
        # Right side
        if segment == segments[2]:
            return Vector(-1,0)
        # Buttom side
        if segment == segments[3]:
            return Vector(0,1)
        raise ValueError('Segment "' + str(segment) + '" is not part of rectangle "' + str(self) + '"')

    # Return a segment which crosses the rectangle in diagonal from bottom left to upper right points
    def get_crossing_segment (self) -> Segment:
        return Segment(self.get_bottom_left_point(), self.get_upper_right_point(), self.segments_color)

    # Get the point in the middle of the rect
    # WARNING: Avoid using this function since it uses the unstable 'Segment.get_middle_point'
    def get_middle_point (self) -> Point:
        # Get the point in the middle of the crossing line
        crossing_segment = self.get_crossing_segment()
        return crossing_segment.get_middle_point()

    # Create a new rectangle identical to this rectangle but with a specified color
    def get_colored_rect(self, segments_color : str = 'black', fill_color : str = 'white') -> 'Rect':
        return Rect(self.x_min, self.y_min, self.x_max, self.y_max, segments_color, fill_color)

    # Split the rect in as many rects as specified by the 'x' and 'y' cuts
    # DANI: No lo he provado desde que lo moví de abajo
    def split (self, x_splits : list = [], y_splits : list = []):
        # Sort the split values and discard those values out of range
        def x_in_range (x):
            return x > self.x_min and x < self.x_max
        formatted_x_splits = list(filter(x_in_range, sorted(x_splits)))
        def y_in_range (y):
            return y > self.y_min and y < self.y_max
        formatted_y_splits = list(filter(y_in_range, sorted(y_splits)))
        # Set the steps to form rectangles
        # WARNING: The doble 'for' loop with generators seems to work wrong
        x_steps = list(pairwise([ self.x_min, *formatted_x_splits, self.x_max ]))
        y_steps = list(pairwise([ self.y_min, *formatted_y_splits, self.y_max ]))
        # Remove steps where both values are identical
        x_steps = [ step for step in x_steps if step[0] != step[1] ]
        y_steps = [ step for step in y_steps if step[0] != step[1] ]
        # Create as many rectangles as required
        for x_min, x_max in x_steps:
            for y_min, y_max in y_steps:
                yield Rect(x_min, y_min, x_max, y_max)

    # Wrap two rectangles by creating a rectangle which contains both
    def wrap (self, rect) -> 'Rect':
        x_min = min(self.x_min, rect.x_min)
        y_min = min(self.y_min, rect.y_min)
        x_max = max(self.x_max, rect.x_max)
        y_max = max(self.y_max, rect.y_max)
        return Rect(x_min, y_min, x_max, y_max)

    # Given another segment, it returns the overlapping segment with this rectangle if exists, as a new segment
    def get_overlap_segment (self, segment : 'Segment') -> Optional['Segment']:
        if segment.is_horizontal():
            y = segment.a.y
            if y < self.y_min or y > self.y_max:
                return None
            # Find the new a and b x coords
            new_a_x = segment.a.x
            if new_a_x < self.x_min:
                new_a_x = self.x_min
            elif new_a_x > self.x_max:
                new_a_x = self.x_max
            new_b_x = segment.b.x
            if new_b_x < self.x_min:
                new_b_x = self.x_min
            elif new_b_x > self.x_max:
                new_b_x = self.x_max
            # If the coords are the same then it means the segment was out of the rectangle
            if new_a_x == new_b_x:
                return None
            new_a = Point(new_a_x, y)
            new_b = Point(new_b_x, y)
            return Segment(new_a, new_b)
        elif segment.is_vertical():
            x = segment.a.x
            if x < self.x_min or x > self.x_max:
                return None
            # Find the new a and b y coords
            new_a_y = segment.a.y
            if new_a_y < self.y_min:
                new_a_y = self.y_min
            elif new_a_y > self.y_max:
                new_a_y = self.y_max
            new_b_y = segment.b.y
            if new_b_y < self.y_min:
                new_b_y = self.y_min
            elif new_b_y > self.y_max:
                new_b_y = self.y_max
            # If the coords are the same then it means the segment was out of the rectangle
            if new_a_y == new_b_y:
                return None
            new_a = Point(x, new_a_y)
            new_b = Point(x, new_b_y)
            return Segment(new_a, new_b)
        else:
            segment_rect = Rect.from_segments([segment])
            overlap_rect = self.get_overlap_rect(segment_rect, borders = False)
            if not overlap_rect:
                return None
            if overlap_rect.get_bottom_left_point() in segment:
                new_a = overlap_rect.get_bottom_left_point()
                new_b = overlap_rect.get_upper_right_point()
                return Segment(new_a, new_b)
            else:
                new_a = overlap_rect.get_bottom_right_point()
                new_b = overlap_rect.get_upper_left_point()
                return Segment(new_a, new_b)

    # Given another rectangle, it returns the overlapping region with this rectangle, if exists, as a new rectangle
    # DANI: No lo he provado desde que lo moví de abajo
    def get_overlap_rect (self, rect : 'Rect', borders : bool = True) -> Optional[ Union[ 'Point', 'Rect', 'Segment'] ]:
        # Find the overlap in the 'x' dimension
        # Get the maximum of the minimums
        x_minimum = max(self.x_min, rect.x_min)
        # Get the minimum of the maximums
        x_maximum = min(self.x_max, rect.x_max)
        # Check that the overlap range exists
        if x_minimum > x_maximum or (x_minimum == x_maximum and not borders):
            return None
        # Find the overlap in the 'y' dimension
        # Get the maximum of the minimums
        y_minimum = max(self.y_min, rect.y_min)
        # Get the minimum of the maximums
        y_maximum = min(self.y_max, rect.y_max)
        # Check that the overlap range exists
        if y_minimum > y_maximum or (y_minimum == y_maximum and not borders):
            return None
        # In case the overlap is happening in a border
        if (x_minimum == x_maximum and y_minimum != y_maximum) or (x_minimum != x_maximum and y_minimum == y_maximum):
            a = Point(x_minimum, y_minimum)
            b = Point(x_maximum, y_maximum)
            return Segment(a,b)
        # In case the overlap is happening in a corner
        if x_minimum == x_maximum and y_minimum == y_maximum:
            return Point(x_minimum, y_minimum)
        # Otherwise, the overlap is a whole rectangle
        # Build a rectangle with both dimensional overlaps
        return Rect(x_minimum, y_minimum, x_maximum, y_maximum)

    # Join 2 rectangles by returning a list of all rectangles which define the resulting grid
    # The overlapped region, if exists is transformed to a single rectangle
    def join_rect (self, rect : 'Rect') -> List['Rect']:
        # Find the overlap between these two rectangles
        overlap = self.get_overlap_rect(rect, borders = False)
        # If there is no overlap then just return both input rectangles
        if not overlap:
            return [self, rect]
        # Split both input rectangles using the overlap maximum and minimum points as split points
        x_splits = [ overlap.x_min, overlap.x_max ]
        y_splits = [ overlap.y_min, overlap.y_max ]
        split1 = self.split(x_splits, y_splits)
        split2 = rect.split(x_splits, y_splits)
        rects = list(split1) + list(split2)
        # Return only unique rectangles
        return set(rects)

    # Substract the second rectangle form the first rectangle and return the rectangles which define the resulting grid
    def subtract_rect (self, rect : 'Rect') -> List['Rect']:
        # Find the overlap between these two rectangles
        overlap = self.get_overlap_rect(rect, borders = False)
        # If there is no overlap then just return the first rectangle intact
        if not overlap:
            return [self]
        # Split the first input rectangle using the overlap maximum and minimum points as split points
        x_splits = [ overlap.x_min, overlap.x_max ]
        y_splits = [ overlap.y_min, overlap.y_max ]
        split = self.split(x_splits, y_splits)
        # Return all splits but the overlapped rectangle
        rects = [ rect for rect in list(split) if rect != overlap ]
        return rects

    # Substract several rectangles from self rectangle and return the rectangles which define the resulting grid
    def subtract_rects (self, rects : List['Rect']) -> List['Rect']:
        # Find the overlap rects
        overlaps = []
        for rect in rects:
            overlap = self.get_overlap_rect(rect, borders = False)
            if overlap:
                overlaps.append(overlap)
        # If there is no overlap then just return the first rectangle intact
        if len(overlaps) == 0:
            return [self]
        return self.subtract_overlaps(overlaps)

    # Substract several rectangles which MUST overlap from self rectangle and return the rectangles which define the resulting grid
    def subtract_overlaps (self, overlaps : List['Rect']) -> List['Rect']:
        # Split the first input rectangle using the overlap maximum and minimum points as split points
        x_splits = [ x_coord for overlap in overlaps for x_coord in [ overlap.x_min, overlap.x_max ] ]
        y_splits = [ y_coord for overlap in overlaps for y_coord in [ overlap.y_min, overlap.y_max ] ]
        split = self.split(x_splits, y_splits)
        # Return all splits but the overlapped rectangle
        def get_overlapped (rect):
            return next((overlap for overlap in overlaps if rect in overlap), None)
        rects = [ rect for rect in list(split) if not get_overlapped(rect) ]
        # Test the grid
        return rects

    # Generate a new rectangle with expanded margins
    def expand_margins (self, margin_size : number) -> 'Rect':
        new_x_min = self.x_min - margin_size
        new_y_min = self.y_min - margin_size
        new_x_max = self.x_max + margin_size
        new_y_max = self.y_max + margin_size
        return Rect(new_x_min, new_y_min, new_x_max, new_y_max)

# A polygon is a list of connected segments which is closed
# IMPORTANT: Segments must follow some standards:
# - All segments follow an order: each segment 'b' point is the next segment 'a' point
# - The polygon is closed: the last segment 'b' point is the first segment 'a' point
# - There are no splitted segments: connected segments are always in different lines (i.e. making a corner)
class Polygon:

    def __init__(self, segments : list, color = 'black'):
        self.segments = segments
        # Check the polygon is closed and segments are not diagonal
        self.check()
        # Save display parameters
        self.color = color
        # Color at this moment all segments
        for segment in segments:
            segment.color = color
        # The size is the length to cross the whole boundary
        self.size = self.get_box().get_crossing_segment().length
        # Set internar variables
        self._clockwise = None
        self._corners = None
        self._grid = None


    # Set a class constant error
    open_polygon_error = ValueError('The polygon is not closed')

    # Set the polygon from segments in a non-canonical format
    # Segments will be sorted, flipped and merged as necessary to accomplish the polygon standard
    @classmethod
    def non_canonical(cls, segments : List[Segment]):

        # Check each polygon segment to not be diagonal
        for segment in segments:
            if segment.is_diagonal():
                raise RuntimeError('The polygon has diagonal segments, which are not supported')

        # ---------------------------------------------------------------------------------
        # Format segments in a way that segments and their points are ordered
        # i.e. each segment 'b' point is the 'a' point of the next segment
        # The first segment 'a' point must be the final segment 'b' point
        # If the list of segments cannot be formatted like this then return an error
        # ---------------------------------------------------------------------------------
        # Check that each segment ends in the same point that the next segment starts
        ordered_segments = [ segments[0] ]
        available_segments = segments[1:]
        while len(ordered_segments) < len(segments):
            # If there are no more available segments return error
            # This may happen in case there was a duplicated segment, which means the polygon is wrong
            if len(available_segments) == 0:
                add_frame(segments)
                raise cls.open_polygon_error
            # Get the last ordered segment to find which is the next connected segment
            last_ordered_segment = ordered_segments[-1]
            last_point = last_ordered_segment.b
            # Get the segment which is connected with the previous segment
            connected_segment = next((segment for segment in available_segments if last_point in segment.points), None)
            if not connected_segment:
                add_frame(segments)
                raise cls.open_polygon_error
            # Remove the connected segment from the available 
            available_segments = [ segment for segment in available_segments if segment != connected_segment ]
            # The connected segment must be connected by the 'a' point
            # If it is connected by the 'b' point then get the inverted segment
            if connected_segment.b == last_point:
                connected_segment = connected_segment.inverted()
            ordered_segments.append(connected_segment)
        # Finally, check that the first segment and the last segment are also connected
        if ordered_segments[0].a != ordered_segments[-1].b:
            add_frame(segments)
            raise cls.open_polygon_error
        segments = ordered_segments
            
        # Merge continuous segments (i.e. connected segments in the same line)
        count = 0
        while count < len(segments):
            for current, nextone in pairwise(segments, retro=True):
                if current.same_line_as(nextone):
                    merged_segment = Segment(current.a, nextone.b)
                    segments = [ seg for seg in segments if seg not in [current, nextone] ]
                    segments.insert(count, merged_segment)
                    count = 0
                    break
                count += 1

        # Build the canonical polygon
        return cls(segments)

    # Set the polygon from its corners in order
    @classmethod
    def from_corners (cls, corners : List[Corner]):
        # Check that there are at least 4 points
        if len(corners) < 4:
            raise RuntimeError('It is required at least 4 points')
        # Set the segment between each pair of points
        segments = []
        for a,b in pairwise(corners, retro=True):
            new_segment = Segment(a,b)
            segments.append(new_segment)
        # Build the polygon
        return cls(segments)

    # Set the polygon from a rectangle
    @classmethod
    def from_rect (cls, rect : Rect):
        return cls(rect.segments)

    def __str__ (self):
        return ', '.join([str(segment) for segment in self.segments])

    def __repr__ (self):
        return ', '.join([str(segment) for segment in self.segments])

    def __lt__(self, other):
        return self.get_box().area < other.get_box().area

    def __eq__ (self, other):
        if isinstance(other, self.__class__) or isinstance(other, Rect):
            return all([ segment in self.segments for segment in other.segments ])
        return False

    def __contains__ (self, other):
        for segment in self.segments:
            if other in segment:
                return True
        return False

    # Check if the poligon is clockwise or counterclockwise
    def is_clockwise (self):
        # If clockwise was previously calculated then return it
        if self._clockwise:
            return self._clockwise
        # Find out if the polygon is clockwise and set is corners
        clockwise, corners = self.setup()
        self._clockwise = clockwise
        self._corners = corners
        return self._clockwise

    # Polygon clockwise (read only)
    clockwise = property(is_clockwise, None, None, "Get if the polygon is clockwise or counterclockwise")

    # Get the polygon corners
    def get_corners (self):
        # If corners were previously calculated then return them
        if self._corners:
            return self._corners
        # Find out if the polygon is clockwise and set is corners
        clockwise, corners = self.setup()
        self._clockwise = clockwise
        self._corners = corners
        return self._corners

    # Polygon corners (read only)
    corners = property(get_corners, None, None, "The polygon corners")

    # Get a corner from its point
    def get_corner (self, point : 'Point') -> Optional['Corner']:
        corners = self.get_corners()
        corner = next((corner for corner in corners if corner == point), None)
        if not corner:
            return None
        return corner

    # Get the polygon grid
    # If grid was previously calculated then return it
    # Otherwise, get all rectangles which form the polygon and set a new grid with them
    def get_grid (self):
        if not self._grid:
            rects = self.split_in_rectangles()
            self._grid = Grid(rects)
        return self._grid

    # Polygon grid (read only)
    grid = property(get_grid, None, None, "The polygon grid")

    # Get the polygon area
    def get_area (self):
        return self.grid.area

    # Polygon area (read only)
    area = property(get_area, None, None, "The polygon area")

    # Get the inside corners
    def get_inside_corners (self):
        return [ corner for corner in self.corners if corner.inside ]

    # Check if the current segments create a closed polygon or if it is open
    def is_closed (self):
        # Check that each segment ends in the same point that the next segment starts
        for current, nextone in pairwise(self.segments, retro=True):
            if current.b != nextone.a:
                return False
        return True

    # Check the polygon to be closed
    # DANI: En principio los perímetros no cerrados no tendrán soporte nunca porque no tienen mucho sentido o no les veo la utilidad
    # Check each polygon segment to not be diagonal
    # Otherwise return an error, since polygons with diagonal segments are not yet supported
    # DANI: En principio algun día se podría hacer esto
    def check (self):
        # Check the polygon to be closed
        if not self.is_closed():
            raise self.open_polygon_error
        # Check each polygon segment to not be diagonal
        for segment in self.segments:
            if segment.is_diagonal():
                raise RuntimeError('The polygon has diagonal segments, which are not supported')

    # Set the polygon corners as points with additional stored values
    # The 'segments' variable defines the two segments of the polygon which form the corner itself
    # The 'inside' variable defines if the corner is "pointing" to the inside of the polygon
    # (e.g. rectangles have no inside corners while an 'L' has one inside corner)
    # DANI: OBSOLETO -> Ahora se usa la función 'setup'
    # DANI: Este sistema no soporta lineas diagonales aunque es más rápido
    def set_corners_obsolete (self):
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
        if abs(difference) != 4:
            # If you see this error there may be splitted segments in your polygon
            # Use the non-canonical class method to set your polygon
            raise RuntimeError('There is something wrong with the polygon')

        # Check if are more corners in the counted direction (true) or the other (false)
        lefted_polygon = difference > 0

        # Now set the real corners and save them in the global variable 'corners'
        # Corners with the minoritarian direction will be set as 'inside = True'
        corners = []
        for precorner in precorners:
            point, segments, lefted_corner = precorner
            inside = lefted_corner != lefted_polygon
            corner = Corner(*point.coords, *segments, inside)
            corners.append(corner)
        return corners

    # Understand which is the inside and which is the outside of the polygon
    # Set the polygon corners as points with additional stored values
    #   - The 'segments' variable defines the two segments of the polygon which form the corner itself
    #   - The 'inside' variable defines if the corner is "pointing" to the inside of the polygon
    #     (e.g. rectangles have no inside corners while an 'L' has one inside corner)
    # Set the polygon sides as segments with additional stored values
    #   - The 'inside_left'
    # In order to know which is the inside and which is the outside we do the following:
    # First we count the angle difference between each pair of segments
    # When the further segment is rotated clockwise respect to the previous segment the difference is positive
    # When the further segment is rotated counterclockwise respect to the previous segment the difference is negative
    # If all differences are added then they will always sum -360 (counterclockwise polygon) or 360 (clockwise polygon)
    #   * Note that the "clockwise" or "counterclockwise" polygon term referes only to the sense of the segments it is made of
    #     It could be inverted easily by exchanging the a and b points of every segment, but this is useless
    # If we have a clockwise polygon then all segments clockwise direction points to the inside of the polygon and vice versa
    # If we have a counterclockwise polygon then all segments counterclockwise direction points to the inside of the polygon and vice versa
    def setup (self) -> Tuple[bool, List['Corner']]:
        # For each segment, save the end point and if the direction of the next segment is left
        # Count how many corners are there in one direction (left in this case)
        precorners = []
        angle_count = 0
        for current, nextone in pairwise(self.segments, retro=True):
            # Given two continue segment vectors, find out the angle between them
            angle = current.get_angle_with(nextone)
            if angle == 0:
                # If you see this error there may be splitted segments in your polygon
                # Use the non-canonical class method to set your polygon
                raise ValueError('There are 2 aligned segments: ' + str(current) + ' and ' + str(nextone))
            clockwise = angle > 0
            angle_count += angle
            point = current.b
            precorners.append((point, [current, nextone], clockwise))

        # There should be always 360 more grades in one direction than in the other (may be clockwise or counterclockwise)
        if abs(angle_count) != 360:
            # If you see this error there may be splitted segments in your polygon
            # Use the non-canonical class method to set your polygon
            add_frame(self.segments)
            raise RuntimeError('There is something wrong with the polygon')

        # Check if are more corners in the counted direction (true) or the other (false)
        clockwise_polygon = angle_count > 0

        # Now set the real corners and save them in the global variable 'corners'
        # Corners with the minoritarian direction will be set as 'inside = True'
        corners = []
        for precorner in precorners:
            point, segments, clockwise = precorner
            inside = clockwise != clockwise_polygon
            corner = Corner(*point.coords, *segments, inside)
            corners.append(corner)
        return clockwise_polygon, corners

    # Get all polygon points
    # DANI: Obsoleto -> Usa las corners
    def get_points (self) -> List['Point']:
        return [ segment.a for segment in self.segments ]

    # Get a rectangle which contains the whole polygon
    def get_box (self) -> 'Rect':
        points = self.get_points()
        x_coords = [ point.x for point in points ]
        y_coords = [ point.y for point in points ]
        x_min = min(x_coords)
        x_max = max(x_coords)
        y_min = min(y_coords)
        y_max = max(y_coords)
        return Rect(x_min, y_min, x_max, y_max)

    # Get the perimeter of the polygon, which is the sum of all its segment lengths
    def get_perimeter (self) -> number:
        return sum([ segment.length for segment in self.segments ])

    # Find out if a point is in the border of the polygon (segments and corners)
    def in_border (self, point : Point) -> bool:
        for corner in self.corners:
            if point == corner:
                return True
        for segment in self.segments:
            if point in segment:
                return True
        return False

    # Get a specific polygon corner or segment by specifying a point that matches this element
    def get_border_element (self, point : Point) -> Optional[Union[Point, Segment]]:
        for corner in self.corners:
            if point == corner:
                return corner
        for segment in self.segments:
            if point in segment:
                return segment
        return None

    # Get a specific polygon segment by specifying a segment that overlaps partial or totally this segment
    def get_segment_from_segment (self, segment : Segment) -> Optional[Segment]:
        for polygon_segment in self.segments:
            overlap = polygon_segment.get_overlap_segment(segment)
            if overlap:
                return polygon_segment
        return None

    # Given a segment which overlaps the polygon, get a vector which is perpendicular to this segment and points into de inside of the polygon
    def get_border_inside (self, segment : Segment) -> Vector:
        polygon_segment = self.get_segment_from_segment(segment)
        if not polygon_segment:
            raise ValueError('The segment is not in the polyigon')
        angle = 90 if self.clockwise else -90
        return polygon_segment.direction.rotate(angle)

    # Return all points in the polygon segments which intersect with a given segment
    # Points are sorted according to their distance with the 'a' point of the segment (from less to more distance)
    # WARNING: Intersection of paralel (overlapping) segments is not detected
    # WARNING: If the origen point is in a corner/segment it may return or not the origen as intersection point
    # DANI: No se ha probado a fondo
    # DANI: Actualmente NO está en uso
    def get_segment_intersection_points (self, segment : Segment) -> Optional[list]:
        # Get the intersection point of the specfied segment with each polygon limit
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

    # Return all segments in the polygon which intersect with a given segment
    # segments are sorted according to their distance with the 'a' point of the segment (from less to more distance)
    # WARNING: Intersection of paralel (overlapping) segments is not detected
    # DANI: No se ha probado
    # DANI: No está acabado. No se contempla la posibilidad de que el segmento que intersecta empieze dentro del perímetro
    # DANI: Actualmente NO está en uso
    def get_segment_intersection_segments (self, segment : Segment) -> Optional[list]:
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

    # Find overlap segments between self polygon segments and other segment
    # If the other segment is inside the area of the polygon it will not be considered
    def get_segment_overlap_segments (self, other : 'Segment') -> Generator[Segment, None, None]:
        for segment in self.segments:
            overlap_segment = segment.get_overlap_segment(other)
            if overlap_segment:
                yield overlap_segment

    # Given a group of segments, find overlaps with the polygon segments
    def get_segments_overlap_segments (self, segments : List[Segment]) -> List[Segment]:
        overall_overlap_segments = []
        for segment in self.segments:
            for other_segment in segments:
                overlap_segments = list(self.get_segment_overlap_segments(other_segment))
                overall_overlap_segments += overlap_segments
        return list(set(overall_overlap_segments))

    # Get polygon segments after substracting the overlap region with a list of segments
    def get_non_overlap_segments (self, segments : List[Segment]) -> List[Segment]:
        overall_non_overlap_segments = []
        for segment in self.segments:
            non_overlap_segments = segment.substract_segments(segments)
            overall_non_overlap_segments += non_overlap_segments
        return overall_non_overlap_segments

    # Get overlapping regions between polygon segments
    # Note that segments inside the area of the polygon will not be considered, just the perimeter
    def get_polygon_overlap_segments (self, other : 'Polygon') -> List[Segment]:
        return self.get_segments_overlap_segments(other.segments)

    # Get non overlapping regions between polygon segments
    # Note that segments inside the area of the polygon will not be considered, just the perimeter
    def get_polygon_non_overlap_segments (self, other : 'Polygon') -> List[Segment]:
        return self.get_non_overlap_segments(other.segments)

    # Check if self polygon is colliding with other polygon
    # i.e. one of their segments is totally or partially overlapping
    # Note that one polygon may be inside of the other
    def is_colliding_with (self, other : 'Polygon') -> bool:
        for self_segment in self.segments:
            for other_segment in other.segments:
                overlap = self_segment.get_overlap_segment(other_segment)
                if overlap:
                    return True
        return False

    # Split the space inside the perimeter in a list of rectangles in a grid friendly format
    def split_in_rectangles (self) -> List[Rect]:

        # Limit points are inside corners
        limit_points = self.get_inside_corners()

        # Limit segments are perimeter segments
        limit_segments = self.segments

        # Given a polygon corner which is pointing inside,
        # Create two segments opposed to the corner segments but as long as required to cut the polygon at onther segment or corner
        # This two segments will always be inside the polygon
        def get_corner_insider_segments (corner : Corner) -> list:

            # Get the oppoiste segments to the corner segments but as long as the boundary size
            # NEVER FORGET: The corner point is the entry segment 'b' point and the exit segment 'a' point
            entry_segment, exit_segment = corner.segments
            tracing1 = Segment(corner, corner + entry_segment.vector.normalized() * self.size, color='green')
            tracing2 = Segment(corner, corner + (-exit_segment.vector.normalized()) * self.size, color='green')

            insider_segments = []
            for segment in [tracing1, tracing2]:

                # Get the intersection point of the specfied segment with each polygon limit
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
                # This may fail, for example, if interior polygons are overlapped
                if len(intersection_points) == 0:
                    # Represent polygons in the problematic boundary
                    # for limit_segment in limit_segments:
                    #     limit_segment.color = 'black'
                    # interior_segments = [ segment for polygon in self.interior_polygons for segment in polygon.segments ]
                    # for interior_segment in interior_segments:
                    #     interior_segment.color = 'blue'
                    # segment.color = 'green'
                    # add_frame(limit_segments + interior_segments + list(mark_point(corner, 'red')) + [segment])
                    raise RuntimeError('An inside segment has no intersection point: ' + str(segment))

                # Sort the points by distance
                def by_distance(point):
                    return segment.a.get_distance_to(point)
                sorted_points = sorted(intersection_points, key=by_distance)

                # The closest point will be the first point
                cut_point = sorted_points[0]

                insider_segment = Segment(segment.a, cut_point, color='red')
                if insider_segment.is_diagonal():
                    # This may happen due to a resolution problem
                    raise RuntimeError('An insider segment should never be diagonal: ' + str(insider_segment))
                insider_segments.append(insider_segment)

            #return [tracing1, tracing2]
            return insider_segments

        # Get all inside separator segments
        inside_separators = []
        for corner in limit_points:
            for segment in get_corner_insider_segments(corner):
                inside_separators.append(segment)

        # Remove duplicates
        inside_separators = unique(inside_separators)

        # Join all inside segments with the polygon limits in a single list
        # WARNING: inside segments MUST be before limit segments
        all_segments = [ *inside_separators, *limit_segments ]

        # DANI: Usa esto para ver los perímetros en negro y los segmentos interiores en rojo
        #add_frame(all_segments)

        # Find all points where the inside separator segments intersect each other
        all_intersections = []
        for segment1, segment2 in itertools.combinations(all_segments, 2):
            intersection = segment1.get_intersection_point(segment2)
            if not intersection:
                continue
            all_intersections.append(intersection)

        # Remove duplicates
        all_intersections = unique(all_intersections)

        # Split the inside separator segments at the intersection points
        all_splitted_segments = []
        for segment in all_segments:
            splitted_segments = segment.split_at_points(all_intersections)
            all_splitted_segments += list(splitted_segments)
        #print('All segments: ' + str(len(all_splitted_segments)))

        # DANI: Usa para ver los segmentos ya cortados
        # for segment in all_splitted_segments:
        #     segment.color = generate_random_color()
        # add_frame(all_splitted_segments)

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
        final_rectangles = unique(final_rectangles)

        #print('Total rectangles: ' + str(len(final_rectangles)))

        return final_rectangles

# A boundary is made of an external polygon and any number of internal polygons which mark the limits of an area
class Boundary:
    def __init__(self, exterior_polygon : Polygon, interior_polygons : List[Polygon] = [], color = 'black', fill_color = 'white'):
        self.exterior_polygon = exterior_polygon
        self.interior_polygons = interior_polygons
        self.polygons = [ exterior_polygon, *interior_polygons ]
        self.segments = [ segment for polygon in self.polygons for segment in polygon.segments ]
        # Calculate some internal values -------------------------------------------------------------------------
        # The size is the length to cross the whole boundary
        self.size = exterior_polygon.get_box().get_crossing_segment().length
        # Inside corners are the exterior polygon inside corners and the interior polygons outside corners
        self.exterior_inside_corners = [ corner for corner in exterior_polygon.corners if corner.inside == True ]
        self.interior_inside_corners = []
        for polygon in interior_polygons:
            self.interior_inside_corners += [ corner for corner in polygon.corners if corner.inside == False ]
        self.inside_corners = self.exterior_inside_corners + self.interior_inside_corners
        # --------------------------------------------------------------------------------------------------------
        # Set some display parameters
        self.color = color
        self.fill_color = fill_color
        # --------------------------------------------------------------------------------------------------------
        self._grid = None
        self._area = None
        # Check interior polygons to be inside the exterior polygon
        self.simple = len(interior_polygons) == 0
        if not self.simple:
            self.check()

    def __contains__(self, other):
        return other in self.grid

    def __eq__(self, other):
        if not other:
            return False
        if self.exterior_polygon != other.exterior_polygon:
            return False
        if len(self.interior_polygons) != len(other.interior_polygons):
            return False
        for self_interior_polygon in self.interior_polygons:
            if self_interior_polygon not in other.interior_polygons:
                return False
        return True

    # Get the boundary grid
    def get_grid (self) -> Optional['Grid']:
        # If grid already exists then return it
        if self._grid:
            return self._grid
        # Otherwise, set the grid
        # Split the boundary in rectangles in a grid-friendly format
        rects = self.split_in_rectangles()
        # Return None when the interior polygons fully consume the exterior polygon
        if len(rects) == 0:
            return None
        # Apply the boundary color to each rectangle
        for rect in rects:
            rect.color = self.fill_color
        # Set, save and return the grid
        grid = Grid(rects)
        grid._boundaries = [self]
        self._grid = grid
        return grid

    # The boundary grid (read only)
    grid = property(get_grid, None, None, "The boundary grid")

    # Get the area inside the polygon (read only)
    def get_area(self):
        return self.grid.area
    area = property(get_area, None, None, "The area inside the boundary")

    # Get a boundary made with the exterior polygon only
    def get_simple_boundary(self):
        if self.simple:
            return self
        return Boundary(self.exterior_polygon)

    # Check all interior polygons to be inside the exterior polygon
    # Note that this function is called when the self grid is not made yet, so checking this is not that easy
    # To do so, we calculate the grid of the exterior polygon alone and then we check if each interior polygon is inside
    def check (self):
        simple_boundary = self.get_simple_boundary()
        for polygon in self.interior_polygons:
            if polygon not in simple_boundary:
                # Represent polygons in the problematic boundary
                exterior_segments = self.exterior_polygon.segments
                for segment in exterior_segments:
                    segment.color = 'black'
                interior_segments = polygon.segments
                for segment in interior_segments:
                    segment.color = 'red'
                add_frame(exterior_segments + interior_segments)
                raise ValueError('At least one interior polygon is not fully inside the exterior polygon')

    # Split the space inside the boundary in a list of rectangles in a grid friendly format
    # It works even if the interior polygons are overlapped
    def split_in_rectangles (self) -> List[Rect]:
        # Calculate the grid from the exterior polygon
        final_grid = self.exterior_polygon.grid
        # Calculate the grid for each interior polygon
        interior_grids = [ polygon.grid for polygon in self.interior_polygons ]
        for interior_grid in interior_grids:
            final_grid = final_grid.get_substract_grid(interior_grid)
        # Return an empty list when the interior polygons fully consume the exterior polygon
        if not final_grid:
            return []
        return final_grid.rects
    
    # Get the overlap boundaries between 2 boundaries
    def get_overlap_boundaries (self, other : 'Boundary') -> List['Boundary']:
        # Get the overlap grid of self boundary grid and the other boundary grid
        overlap_grid = self.grid.get_overlap_grid(other.grid)
        return overlap_grid.boundaries

    # Fuse other boundary to self boundary
    # Check that both boundaries can be joined as a single boundary
    # i.e. both exterior polygons must be colliding and the colliding region/s must be as wide as the minimum size or more
    def merge_boundary (self, other : 'Boundary', min_size = None) -> 'Boundary':
        # Check if both boundaries are colliding
        # i.e. both external polygons have overlapping segments
        colliding_segments = self.exterior_polygon.get_polygon_overlap_segments(other.exterior_polygon)
        if len(colliding_segments) == 0:
            raise ValueError('Boundaries are not colliding')
        # In case there is a minimum size restriction check that the colliding segment is as long as required
        if min_size and all([ segment.length < min_size for segment in colliding_segments ]):
            raise ValueError('Some colliding region is not wide enough')
        # Check that there are not overlapping areas between both boundaries
        overlaps = self.get_overlap_boundaries(other)
        if len(overlaps) == 0:
            raise ValueError('Boundaries are overlapping')
        # Finally, merge both boundary grids into a single grid and get its boundary
        grid = self.grid.get_merge_grid(other.grid, check_overlaps=False)
        boundaries = grid.boundaries
        if len(boundaries) != 1:
            raise ValueError('Something went wrong with the boundary')
        return boundaries[0]

    # Check if self boundary exterior polygon is colliding with other boundary exterior polygon/s
    # i.e. one of their segments is totally or partially overlapping
    # Note that one polygon may be inside of the other
    def is_colliding_with (self, other : 'Boundary') -> bool:
        return self.exterior_polygon.is_colliding_with(other.exterior_polygon)
    def is_colliding_with_any (self, others : List['Boundary']) -> bool:
        for other in others:
            if self.exterior_polygon.is_colliding_with(other.exterior_polygon):
                return True
        return False

# A grid is a group or groups of rectangles connected according to a standard:
# Two connected rectangles have the whole segment connected
# i.e. both rectangles have an identical segment
# This class is used to handle space inside boundaries
class Grid:
    def __init__(self, rects : List[Rect]):
        self._rects = rects
        # Check rectangles to match the grid format requirements
        self.check()
        self._max_rects = None
        self._rows = None
        self._columns = None
        self._area = None
        self._boundaries = None

    def __str__(self):
        return str(self.rects)

    def __repr__(self):
        return str(self.rects)

    def __contains__(self, other):
        if isinstance(other, Point):
            for rect in self.rects:
                if other in rect:
                    return True
            return False
        if isinstance(other, Segment):
            in_a = other.a in self
            cross_any_segment = any([ other.get_intersection_point(segment, in_extremis=0) for segment in self.get_boundary_segments() ])
            return in_a and not cross_any_segment
        if isinstance(other, Rect):
            return all([ segment in self for segment in other.get_segments() ])
        if isinstance(other, Polygon):
            return all([ segment in self for segment in other.segments ])
        if isinstance(other, Boundary):
            return all([ segment in self for segment in other.exterior_polygon.segments ])
        if isinstance(other, self.__class__):
            return all([ boundary in self for boundary in other.boundaries ])
        return False

    def __add__(self, other : 'Grid') -> 'Grid':
        return self.get_merge_grid(other)

    def __sub__(self, other : 'Grid') -> 'Grid':
        return self.get_substract_grid(other)

    # Check for each corner on each rect that, if it is inside any other rect, it is also a corner in this rect
    # Check also that rects to do not overlap, since this may happen even with no corners inside each other
    def check (self):
        for rect, other_rects in otherwise(self.rects):
            rect_corners = rect.get_points()
            for other_rect in other_rects:
                other_corners = other_rect.get_points()
                for corner in rect_corners:
                    if corner in other_rect and corner not in other_corners:
                        print('WARNING: Conflict rects ' + str(rect) + ' and ' + str(other_rect))
                        raise RuntimeError('Grid rects are wrong')
                if rect.get_overlap_rect(other_rect, borders = False):
                    #add_frame(self.rects)
                    print('WARNING: Overlapping rects ' + str(rect) + ' and ' + str(other_rect))
                    raise RuntimeError('Grid rects are wrong')

    # Find redundant columns/rows of rectangles and merge them to make the grid more efficient
    # This should be necessary only when the grid is built in a non canonical way
    # Note that this function modifies the current grid rects
    def compact (self):
        # Start by checking the rows
        new_rows = []
        row_pool = [ *self.rows ]
        while len(row_pool) > 0:
            # Get the next first element in the pool
            current_row = row_pool[0]
            del row_pool[0]
            # Find another element in the pool which is connected to the current element
            while True:
                current_row_overall_rectangle, current_row_contained_rects = current_row
                next_row = next((row for row in row_pool if row[0].is_connected_with(current_row_overall_rectangle)), None)
                # If there are no more rows connected to the current next row then we are done
                if not next_row:
                    break
                row_pool.remove(next_row)
                # Join the next row to the current row
                next_row_overall_rectangle, next_row_contained_rects = next_row
                new_row_overall_rectangle = current_row_overall_rectangle.wrap(next_row_overall_rectangle)
                new_row_contained_rects = []
                for current_row_contained_rect in current_row_contained_rects:
                    connected_contained_rect = next(rect for rect in next_row_contained_rects if rect.is_connected_with(current_row_contained_rect))
                    merged_rect = current_row_contained_rect.wrap(connected_contained_rect)
                    new_row_contained_rects.append(merged_rect)
                current_row = new_row_overall_rectangle, new_row_contained_rects
            # Add the current row, including all the merges (if any) to the list of new rows
            new_rows.append(current_row)
        # Update the grid rects with the new merged rects
        self._rects = sum([ new_row[1] for new_row in new_rows ], [])
        # Reset the columns just in case they were previously calculated, since they will be not valid anymore
        self._columns = None
        # Now repeat the whole process with columns
        new_columns = []
        column_pool = [ *self.columns ]
        while len(column_pool) > 0:
            # Get the next first element in the pool
            current_column = column_pool[0]
            del column_pool[0]
            # Find another element in the pool which is connected to the current element
            while True:
                current_column_overall_rectangle, current_column_contained_rects = current_column
                next_column = next((column for column in column_pool if column[0].is_connected_with(current_column_overall_rectangle)), None)
                # If there are no more columns connected to the current next column then we are done
                if not next_column:
                    break
                column_pool.remove(next_column)
                # Join the next column to the current column
                next_column_overall_rectangle, next_column_contained_rects = next_column
                new_column_overall_rectangle = current_column_overall_rectangle.wrap(next_column_overall_rectangle)
                new_column_contained_rects = []
                for current_column_contained_rect in current_column_contained_rects:
                    connected_contained_rect = next(rect for rect in next_column_contained_rects if rect.is_connected_with(current_column_contained_rect))
                    merged_rect = current_column_contained_rect.wrap(connected_contained_rect)
                    new_column_contained_rects.append(merged_rect)
                current_column = new_column_overall_rectangle, new_column_contained_rects
            # Add the current column, including all the merges (if any) to the list of new columns
            new_columns.append(current_column)
        # Update the grid rects with the new merged rects
        self._rects = sum([ new_column[1] for new_column in new_columns ], [])
        # Now we can keep the previously calculated columns as the current grid columns
        # However we can not keep the previous calculated rows, since their rows may have been merged
        # So we just clean the curren grid rows, so they will be recalcualted if needed
        self._rows = None
        self._columns = new_columns

    # Set the grid from rectangles which do not follow the standards
    # i.e. they are not connected each other by the whole segment
    # However these rectangles must never overlap each other
    # These rectangles will be used to build a new grid with different rectangles which do follow the standards
    @classmethod
    def non_canonical (cls, rects : List[Rect]):
        # There must be at least one rectangle
        if len(rects) == 0:
            raise ValueError('There must be at least one rect to build a grid')
        # First of all find rectangle groups
        # This is the non-canonical version of find_connected_rect_groups method
        rects_to_group = [ *rects ]
        rect_groups = []
        while len(rects_to_group) > 0:
            first_rect = rects_to_group[0]
            new_group = [ first_rect ]
            for current_rect in new_group:
                connected_rects = [ rect for rect in rects if rect.is_colliding_with(current_rect) ]
                new_group += [ rect for rect in connected_rects if rect not in new_group ]
            rects_to_group = [ rect for rect in rects_to_group if rect not in new_group ]
            rect_groups.append(new_group)
        # Now split all rectangles in every x and y limits of all rectangles
        canonical_rects = []
        for rect_group in rect_groups:
            x_splits = sum([ [ rect.x_min, rect.x_max ] for rect in rect_group ], [])
            x_splits = list(set(x_splits))
            y_splits = sum([ [ rect.y_min, rect.y_max ] for rect in rect_group ], [])
            y_splits = list(set(y_splits))
            for rect in rect_group:
                splitted_rects = rect.split(x_splits, y_splits)
                canonical_rects += splitted_rects
        # Now that we have canonical rects, set the grid
        grid = cls(canonical_rects)
        # Now rectangles should be canonical, but we may have redundant "cuts" where several rows/columns may be merged to optimize the grid
        # DANI: Igual hay replantearse esto porque lo mismo es más ineficiente el chekeo + arreglo de la ineficiencia que la propia ineficiencia
        grid.compact()
        return grid

    # Grid rectangles are read only
    def get_rects (self) -> List[Rect]:
        return self._rects
    rects = property(get_rects, None, None, "Grid rectangles")

    # Get maximum possible rectangles in the grid
    # This variable is treated appart since its calculation may be an expensive calculation
    def get_max_rects (self) -> List[ Tuple [ Rect, List[Rect] ] ]:
        # If maximum rectangles are previously calculated then return them
        if self._max_rects:
            return self._max_rects
        # Calculate all maximum rectangles
        self._max_rects = self.find_maximum_rectangles()
        return self._max_rects
    # Grid maximum rectangles (read only)
    max_rects = property(get_max_rects, None, None, "Maximum possible rectangles in the grid")

    # Get row rectangles in the grid
    # This variable is treated appart since its calculation may be an expensive calculation
    def get_rows (self) -> List[ Tuple [ Rect, List[Rect] ] ]:
        # If row rectangles are previously calculated then return them
        if self._rows:
            return self._rows
        # Calculate all row rectangles
        self._rows = self.find_rows()
        return self._rows
    # Grid row rectangles (read only)
    rows = property(get_rows, None, None, "Row rectangles in the grid")

    # Get column rectangles in the grid
    # This variable is treated appart since its calculation may be an expensive calculation
    def get_columns (self) -> List[ Tuple [ Rect, List[Rect] ] ]:
        # If column rectangles are previously calculated then return them
        if self._columns:
            return self._columns
        # Calculate all column rectangles
        self._columns = self.find_columns()
        return self._columns
    # Grid column rectangles (read only)
    columns = property(get_columns, None, None, "Column rectangles in the grid")

    # Get the area of the whole grid
    # This variable is treated appart since its calculation may be an expensive calculation
    def get_area (self) -> number:
        # If the area is previously calculated then return it
        if self._area:
            return self._area
        # Otherwise, calculate the area
        # Add the area of all grid rectangle
        rect_areas = [ rect.get_area() for rect in self.rects ]
        self._area = resolute(sum(rect_areas))
        return self._area
    # Grid area (read only)
    area = property(get_area, None, None, "The area of the whole grid")

    # Find all grid boundaries
    # i.e. find the enveloping boundaries for each group of connected rects
    def find_boundaries (self) -> List[Boundary]:
        boundaries = []
        # First of all isolate groups of connected rectangles
        rect_groups = self.find_connected_rect_groups()
        # Now find the boundary of each group
        for group in rect_groups:
            # Get all rectangle segments and find those which are not duplicated
            # Each non duplicated segment is an outter segment and thus it is part of the boundary
            segments = []
            for rect in group:
                segments += rect.segments
            boundary_segments = [ segm for segm in segments if segments.count(segm) == 1 ]
            # Now get groups of connected segments (i.e. polygons)
            # A group must contain an exterior polygon and it may contain also several interior polygons
            polygons = list(connect_segments(boundary_segments))
            # The external polygon is the polygon with the biggest box
            sorted_polygons = sorted(polygons)
            exterior_polygon = sorted_polygons[-1]
            interior_polygons = sorted_polygons[0:-1]
            # Now create the boundary
            boundary = Boundary(exterior_polygon, interior_polygons)
            boundaries.append(boundary)
        return boundaries

    # Get grid rectangle groups boundaries
    def get_boundaries (self) -> List[Boundary]:
        if not self._boundaries:
            self._boundaries = self.find_boundaries()
        return self._boundaries

    # Grid boundaries (read only)
    boundaries = property(get_boundaries, None, None, "Boundaries which wrap the whole grid")

    # Get all boundaries in all segments
    def get_boundary_segments (self):
        return [ segment for boundary in self.boundaries for segment in boundary.segments ]

    # Check if a rect is overlapping at some rect in this grid
    def is_rect_overlapping (self, rect : Rect) -> bool:
        for self_rect in self.rects:
            overlap = self_rect.get_overlap_rect(rect, borders = False)
            if overlap:
                return True
        return False

    # Return the overlap space between two grids splitted in rectangles
    # Note that these rectangles will not follow the grid standard rules
    def get_overlap_rects (self, grid : 'Grid') -> List[Rect]:
        overlap_rects = []
        for rect in self.rects:
            for other in grid.rects:
                overlap_rect = rect.get_overlap_rect(other, borders = False)
                if overlap_rect:
                    overlap_rects.append(overlap_rect)
        return overlap_rects

    # Return the resulting grid after overlapping self grid to other grid
    def get_overlap_grid (self, grid : 'Grid') -> Optional['Grid']:
        overlap_rects = self.get_overlap_rects(grid)
        if len(overlap_rects) == 0:
            return None
        return Grid.non_canonical(overlap_rects)

    # Return the merge space between two grids splitted in rectangles
    # Note that these rectangles will not follow the grid standard rules
    # DANI: Esto antes se usaba en el get_merge_grid pero ahora ya no se usa en ningún lado
    def get_merge_rects (self, grid : 'Grid') -> List[Rect]:
        merge_rects = []
        for rect in self.rects:
            for other in grid.rects:
                join_rects = rect.join_rect(other)
                merge_rects += join_rects
        return unique(merge_rects)

    # Return the resulting grid after merging self grid to other grid
    # If you are sure grids do not overlap then set check_overlaps as False to avoid a lot of calculus
    def get_merge_grid (self, grid : 'Grid', check_overlaps : bool = True) -> 'Grid':
        overlap_rects = None
        if check_overlaps:
            overlap_rects = self.get_overlap_rects(grid)
        # Get the rects of both grids together
        # Note that rects will not follow the grid standards
        if overlap_rects:
            difference = grid - self
            merge_rects = self.rects + difference.rects
        else:
            merge_rects = self.rects + grid.rects
        return Grid.non_canonical(merge_rects)

    # Return the overlaps between this grid and another grid in a specific format
    # Return a dict where keys are self rectangles and values are the overlaps in self rectangles
    def get_overlaps (self, grid : 'Grid') -> dict:
        overlaps = {}
        for self_rect in self.rects:
            # Find the overlap rects
            overlap_rects = []
            for grid_rect in grid.rects:
                overlap_rect = self_rect.get_overlap_rect(grid_rect, borders = False)
                if overlap_rect:
                    overlap_rects.append(overlap_rect)
            # If there is no overlap then just return the first rectangle intact
            if len(overlap_rects) == 0:
                continue
            # Save the overlaps
            overlaps[self_rect] = overlap_rects
        return overlaps

    # Return self grid space after substracting other grid space splitted in rectangles
    # Note that these rectangles will not follow the grid standard rules
    def get_substract_overlaps (self, overlaps : dict) -> List[Rect]:
        substract_rects = []
        for rect in self.rects:
            overlap = overlaps.get(rect, None)
            if overlap:
                new_rects = rect.subtract_overlaps(overlap)
                substract_rects += new_rects
                continue
            substract_rects.append(rect)
        return substract_rects

    # Return self grid space after substracting other grid space splitted in rectangles
    # Note that these rectangles will not follow the grid standard rules
    # DEPRECATED
    def get_substract_rects (self, grid : 'Grid') -> List[Rect]:
        substract_rects = []
        for rect in self.rects:
            new_rects = rect.subtract_rects(grid.rects)
            substract_rects += new_rects
        return substract_rects

    # Return the resulting grid after substracting other grid to self grid
    # Return None if the other grid fully consumes self grid
    def get_substract_grid (self, grid : 'Grid') -> Optional['Grid']:
        #substract_rects = self.get_substract_rects(grid)
        overlaps = self.get_overlaps(grid)
        if not overlaps:
            return self
        substract_rects = self.get_substract_overlaps(overlaps)
        if len(substract_rects) == 0:
            return None
        return Grid.non_canonical(substract_rects)

    # Search all maximum rectangles which fulfill the specified minimum x and y sizes
    # Fitting rects are returned as a generator
    # If only the x size parameter is passed it is assumed to be both x and y size
    def get_fitting_space (self, x_fit_size : number, y_fit_size : Optional[number] = None) -> Generator[Rect, None, None]:
        if not y_fit_size:
            y_fit_size = x_fit_size
        for max_rect in self.max_rects:
            rect = max_rect[0]
            x_size, y_size = rect.get_size()
            if x_fit_size <= x_size and y_fit_size <= y_size:
                yield rect

    # Some functions are defined to find colliding rects
    # Note that these functions rely on the fact that grids must never have overlapping rects
    def get_left_rect (self, rect : Rect) -> Optional[Rect]:
        x_max = rect.x_min
        y_max = rect.y_max
        for r in self.rects:
            if r.x_max == x_max and r.y_max == y_max:
                return r
        return None
    def get_right_rect (self, rect : Rect) -> Optional[Rect]:
        x_min = rect.x_max
        y_min = rect.y_min
        for r in self.rects:
            if r.x_min == x_min and r.y_min == y_min:
                return r
        return None
    def get_upper_rect (self, rect : Rect) -> Optional[Rect]:
        x_min = rect.x_min
        y_min = rect.y_max
        for r in self.rects:
            if r.x_min == x_min and r.y_min == y_min:
                return r
        return None
    def get_bottom_rect (self, rect : Rect) -> Optional[Rect]:
        x_max = rect.x_max
        y_max = rect.y_min
        for r in self.rects:
            if r.x_max == x_max and r.y_max == y_max:
                return r
        return None
    def get_connected_rects (self, rect : Rect) -> List[Rect]:
        left_rect = self.get_left_rect(rect)
        right_rect = self.get_right_rect(rect)
        upper_rect = self.get_upper_rect(rect)
        bottom_rect = self.get_bottom_rect(rect)
        return [ rect for rect in [ left_rect, right_rect, upper_rect, bottom_rect ] if rect != None ]

    # Check the grid to respect minimum size in both x and y dimensions
    def check_minimum (self, minimum : number) -> bool:
        # Check all maximum rectangles
        for max_rect in self.max_rects:
            rect = max_rect[0]
            x_size, y_size = rect.get_size()
            # If there is at least one maximum rectangle which does not respect minimum size in both x and y dimensions then it is wrong
            if lower(x_size, minimum) and lower(y_size, minimum):
                return False
            # Ignore rectangles which respect size only in one dimension
            if lower(x_size, minimum) or lower(y_size, minimum):
                continue
            # For all rectangles which respect size in both dimensions,
            # Find their free borders: border regions which are not overlapping with the grid boundaries
            limits = [ segment for boundary in self.boundaries for segment in boundary.segments ]
            free_borders = []
            for segment in rect.segments:
                free_borders += segment.substract_segments(limits)            
            # Check free borders to respect the minimum size
            for segment, other_segments in otherwise(free_borders):
                # If the segment is wide enough then skip it
                if not lower(segment.length, minimum):
                    continue
                # If the segment is not wide enough it may make a corner with other free regions thus beeing correct
                # Find out if the free border is making a corner
                corner_segment = next((other for other in other_segments if other.makes_corner_with(segment)), None)
                # In case it is, use the hipotenuse of the border to check the minimum size
                if corner_segment:
                    size = sqrt( segment.length**2 + corner_segment.length**2 )
                    if size >= minimum:
                        continue
                # Otherwise, we are not respecting the minimum size
                return False
        return True

    # This function is used to generate a new grid without all regions which do not respect a minimum size
    def keep_minimum (self, minimum : number) -> Optional['Grid']:
        # Find the maximum size both in x and y for each rect in the grid
        # To do so, find all maximum rects which include a specific rect and get the maximum sizes among them
        respecting_rects = []
        for rect in self.rects:
            max_rects = [ max_rect[0] for max_rect in self.max_rects if rect in max_rect[1] ]
            max_x_size = max([ rect.get_x_size() for rect in max_rects ])
            max_y_size = max([ rect.get_y_size() for rect in max_rects ])
            # If any of the sizes is not enough to cover the size then set the rect as non respecting
            if max_x_size < minimum or  max_y_size < minimum:
                continue
            respecting_rects.append(rect)
        # Make a new grid from the minimum rects which do respect the minimum size
        # Use self grid in case all self rects are respecting the minimum size
        # Return None if there are not respecting rects at this point
        if len(respecting_rects) == 0:
            return None
        respecting_region = self if len(self.rects) == len(respecting_rects) else Grid(respecting_rects)
        # Now, from the respecting regions, find groups of connected rects which respect the minimum size at the connections
        # Note that there may be more than one region (even when there are not regions which do not respect it)
        # Note that these regions may even overlap (see figure 01)
        pool = [ *respecting_region.max_rects ]
        minimum_regions = []
        while len(pool) > 0:
            next_group = [ pool[0] ]
            del pool[0]
            # Find connected rects to the current group of rects until there are no more
            while True:
                next_rect = None
                for max_rect in pool:
                    for group_rect in next_group:
                        # Get the overlap between the current maximum rectangle and the next rectangle in the group
                        overlap = max_rect[0].get_overlap_rect(group_rect[0])
                        if not overlap:
                            continue
                        # If the overlap is a point then rects are no connected
                        if isinstance(overlap, Point):
                            continue
                        # Find the overlap size
                        size = 0
                        if isinstance(overlap, Segment):
                            size = overlap.length
                        if isinstance(overlap, Rect):
                            size = overlap.get_diagonal_size()
                        # If the size is equal or greater than the minimum then we join this rect to the group
                        if size >= minimum:
                            next_rect = max_rect
                            break
                    # Exit the second for loop if we already found a new max_rect
                    if next_rect:
                        break
                # If we found a new rect then keep searching
                if next_rect:
                    next_group.append(next_rect)
                    pool.remove(next_rect)
                    continue
                # If there are not more connected rects then we are done
                # If the number of max rects in this group matches the number of max rects in the respecting region then we are done
                # It means the whole respecting region grid is to be returned (which may be self, so no change would be done)
                if len(next_group) == len(respecting_region.max_rects):
                    return respecting_region
                # Append the current group to the groups list. Save only the minimum rects as a grid
                minimum_rects = list(set(sum([ max_rect[1] for max_rect in next_group ], [])))
                minimum_region = Grid(minimum_rects) 
                minimum_regions.append(minimum_region)
                break
        # DANI: Esto no se ha provado
        print('DANI: Hay más de una minimum region, comprueba que todo esté bien')
        # At this point must be always more than one minimum region
        # Sort minimum regions by area (bigger goes first)
        minimum_regions.sort(key=lambda x: x.area, reverse=True)
        # Return the bigger minimum region
        return minimum_regions[0]

    # One by one for each *available rectangle, where available rectangles are the splitted rectangles
    # Get as many rectanges as possible which are connected horizontally to the current rectangle
    # Get as many rows of rectangles as possible which are connected vertically to all previous rectangles
    # Consider all previous rectangles as a single rectange
    # Repeat in the inverse order (first vertically, then horizontally)
    # Remove all previous rectangles from the *available rectangles list
    # Returned maximum rectangles contain both the overall maximum rectangle and the contained rectangles on it
    def find_maximum_rectangles (self) -> List[ Tuple [ Rect, List[Rect] ] ]:

        # Trace which rectangles have been already checked both horizontally and vertically to avoid repeating
        for rect in self.rects:
            rect.horizontal_check = False
            rect.vertical_check = False

        # Save all maximum rectangles in a list to be returned at the end
        maximum_rectangles = []

        # Group rectangles first horizontally and then vertically
        for rect in self.rects:
            # Skip already grouped rectangles horizontally
            if rect.horizontal_check == True:
                continue
            # Set the first row
            first_row = [ rect ]
            # Append all rectangles at right from current rectangle to the row
            rightest = rect
            while (True):
                rightest = self.get_right_rect(rightest)
                if not rightest:
                    break
                first_row.append(rightest)
            # Append all rectangles at left from current rectangle to the row
            leftest = rect
            while (True):
                leftest = self.get_left_rect(leftest)
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
                new_row = [ self.get_bottom_rect(rect) for rect in current_row ]
                if not all(new_row):
                    break
                group += new_row
                current_row = [ *new_row ]
            # Repeat the process but upperwards
            current_row = [ rect for rect in first_row ]
            while (True):
                new_row = [ self.get_upper_rect(rect) for rect in current_row ]
                if not all(new_row):
                    break
                group += new_row
                current_row = [ *new_row ]
            # Create a new rect which contains all group rects
            maximum_rect = merge_rectangles(group)
            # Add the new maximum rectnagle to the list if it is not there already
            if maximum_rect not in maximum_rectangles:
                maximum_rectangles.append((maximum_rect, group))

        # Now repeat the process in the inverse order (first vertically, then horizontally)
        for rect in self.rects:
            # Skip already grouped rectangles vertically
            if rect.vertical_check == True:
                continue
            # Set the first column
            first_column = [ rect ]
            # Append all rectangles at the bottom from current rectangle to the column
            bottomest = rect
            while (True):
                bottomest = self.get_bottom_rect(bottomest)
                if not bottomest:
                    break
                first_column.append(bottomest)
            # Append all rectangles upper from current rectangle to the column
            upperest = rect
            while (True):
                upperest = self.get_upper_rect(upperest)
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
                new_column = [ self.get_right_rect(rect) for rect in current_column ]
                if not all(new_column):
                    break
                group += new_column
                current_column = [ *new_column ]
            # Repeat the process but to the left
            current_column = [ rect for rect in first_column ]
            while (True):
                new_column = [ self.get_left_rect(rect) for rect in current_column ]
                if not all(new_column):
                    break
                group += new_column
                current_column = [ *new_column ]
            # Create a new rect which contains all group rects
            maximum_rect = merge_rectangles(group)
            # Add the new maximum rectnagle to the list if it is not there already
            if maximum_rect not in maximum_rectangles:
                maximum_rectangles.append((maximum_rect, group))

        return maximum_rectangles

    # One by one for each *available rectangle, where available rectangles are the splitted rectangles
    # Get as many rectanges as possible which are connected horizontally to the current rectangle
    # Consider all previous rectangles as a single rectange
    # Remove all previous rectangles from the *available rectangles list
    # Returned rows contain both the overall row rectangle and the contained rectangles in the row
    def find_rows(self) -> List[ Tuple [ Rect, List[Rect] ] ]:

        # Trace which rectangles have been already checked both horizontally and vertically to avoid repeating
        for rect in self.rects:
            rect.check = False

        # Save all maximum rectangles in a list to be returned at the end
        rows = []

        # Group rectangles first horizontally and then vertically
        for rect in self.rects:
            # Skip already grouped rectangles horizontally
            if rect.check == True:
                continue
            # Start the row
            row = [ rect ]
            # Append all rectangles at right from current rectangle to the row
            rightest = rect
            while (True):
                rightest = self.get_right_rect(rightest)
                if not rightest:
                    break
                row.append(rightest)
            # Append all rectangles at left from current rectangle to the row
            leftest = rect
            while (True):
                leftest = self.get_left_rect(leftest)
                if not leftest:
                    break
                row.append(leftest)
            # Set all rectangles in the row as checked
            for rect in row:
                rect.check = True
            # Merge all rectangles in row
            row_rectangle = merge_row_rectangles(row)
            # Save both the final row rectangle and the rectangles contained in the row
            rows.append((row_rectangle, row))

        return rows

    # One by one for each *available rectangle, where available rectangles are the splitted rectangles
    # Get as many rectanges as possible which are connected vertically to the current rectangle
    # Consider all previous rectangles as a single rectange
    # Remove all previous rectangles from the *available rectangles list
    # Returned columns contain both the overall column rectangle and the contained rectangles in the column
    def find_columns(self) -> List[ Tuple [ Rect, List[Rect] ] ]:

        # Trace which rectangles have been already checked both horizontally and vertically to avoid repeating
        for rect in self.rects:
            rect.check = False

        # Save all maximum rectangles in a list to be returned at the end
        columns = []

        # Group rectangles first horizontally and then vertically
        for rect in self.rects:
            # Skip already grouped rectangles horizontally
            if rect.check == True:
                continue
            # Start the column
            column = [ rect ]
            # Append all rectangles at the bottom from current rectangle to the column
            bottomest = rect
            while (True):
                bottomest = self.get_bottom_rect(bottomest)
                if not bottomest:
                    break
                column.append(bottomest)
            # Append all rectangles upper from current rectangle to the column
            upperest = rect
            while (True):
                upperest = self.get_upper_rect(upperest)
                if not upperest:
                    break
                column.append(upperest)
            # Set all rectangles in the column as checked
            for rect in column:
                rect.check = True
            # Merge all rectangles in column
            column_rectangle = merge_column_rectangles(column)
            columns.append((column_rectangle, column))

        return columns

    # Find groups of connected rectangles
    def find_connected_rect_groups (self):
        rects_to_group = self.rects
        rect_groups = []
        while len(rects_to_group) > 0:
            first_rect = rects_to_group[0]
            new_group = [ first_rect ]
            for current_rect in new_group:
                connected_rects = self.get_connected_rects(current_rect)
                new_group += [ rect for rect in connected_rects if rect not in new_group ]
            rects_to_group = [ rect for rect in rects_to_group if rect not in new_group ]
            rect_groups.append(new_group)
        return rect_groups

    # Given a segment, get all segment regions which overlap the grid
    def get_segment_overlap_segments (self, segment : Segment) -> List[Segment]:
        # Get the overlap of the segment with every rect in the grid
        overlap_segments = []
        for rect in self.rects:
            overlap_segment = rect.get_overlap_segment(segment)
            if overlap_segment:
                overlap_segments.append(overlap_segment)
        # Now merge all the overlap segmetns as much as possible
        overlap_segments = merge_inline_segments(overlap_segments)
        return overlap_segments
        
    # Given a list of segments, get all segment regions which overlap the grid
    def get_segments_overlap_segments (self, segments : List[Segment]) -> List[Segment]:
        overall_overlap_segments = []
        for segment in segments:
            overlap_segments = self.get_segment_overlap_segments(segment)
            overall_overlap_segments += overlap_segments
        return overall_overlap_segments

    # Given a polygon, for each segment in the polygon, get all segment regions which overlap the grid
    def get_polygon_overlap_segments (self, polygon : Polygon) -> List[Segment]:
        return self.get_segments_overlap_segments(polygon.segments)


# Auxiliar functions ---------------------------------------------------------------

# Return unique values of a list
def unique (values : list) -> list:
    return list(set(values))

# Set a special iteration system
# Return each value of the array and the next value
# By default, the final array value is skipped, since it has no next value
# However, if the 'retro' argument is True, the final array value is returned with the first array value as the next value
# By default, values are returned as follows: A with B, B with C, C with D ...
# However, if the 'loyals' argument is True, values are returned as follows: A with B, C with D, E with F ...
def pairwise (values : list, retro : bool = False, loyals : bool = False) -> Generator[tuple, None, None]:
    last = len(values) - 1
    step = 2 if loyals else 1
    for i in range(0, last, step):
        yield values[i], values[i+1]
    if retro:
        yield values[last], values[0]

# Set a special iteration system
# Return one value of the array and a new array with all other values for each value
def otherwise (values : list) -> Generator[tuple, None, None]:
    for v, value in enumerate(values):
        others = values[0:v] + values[v+1:]
        yield value, others

# x/y coordinates sorters
def sort_by_x (point : Point) -> number:
    return point.x
def sort_by_y (point : Point) -> number:
    return point.y
# Point sorter for points in a same line (otherwise it makes not sense)
def sort_points (points : List[Point]) -> List[Point]:
    return sorted( sorted( points, key=sort_by_x), key=sort_by_y)

# Set a function to merge multiple rects into a single big rect
# Do it by finding the most maximum x and y values
def merge_rectangles (rects : List[Rect]) -> Rect:
    min_x_min = min([ rect.x_min for rect in rects ])
    min_y_min = min([ rect.y_min for rect in rects ])
    max_x_max = max([ rect.x_max for rect in rects ])
    max_y_max = max([ rect.y_max for rect in rects ])
    return Rect(min_x_min, min_y_min, max_x_max, max_y_max)

# Set a function to merge multiple rects in a row into a single big rect
# Do it by finding the most maximum x and y values
# Note that all rects will match in both min and max y values
def merge_row_rectangles (rects : List[Rect]) -> Rect:
    min_x_min = min([ rect.x_min for rect in rects ])
    min_y_min = rects[0].y_min
    max_x_max = max([ rect.x_max for rect in rects ])
    max_y_max = rects[0].y_max
    return Rect(min_x_min, min_y_min, max_x_max, max_y_max)

# Set a function to merge multiple rects in a column into a single big rect
# Do it by finding the most maximum x and y values
# Note that all rects will match in both min and max x values
def merge_column_rectangles (rects : List[Rect]) -> Rect:
    min_x_min = rects[0].x_min
    min_y_min = min([ rect.y_min for rect in rects ])
    max_x_max = rects[0].x_max
    max_y_max = max([ rect.y_max for rect in rects ])
    return Rect(min_x_min, min_y_min, max_x_max, max_y_max)

# Given a group of segments which must be in the same line, merged them into bigger segments when connected
def merge_inline_segments (segments : List[Segment]) -> List[Segment]:
    remaining_segments = [ *segments ]
    merged_segments = []
    while len(remaining_segments) > 0:
        first_segment = remaining_segments[0]
        remaining_segments.remove(first_segment)
        merged_segment = first_segment
        while True:
            # Find a new segment which is connected to the merged segment
            connected_segment = next((segment for segment in remaining_segments if segment.is_connected_with(merged_segment)), None)
            # If there is no next segment then we are done with the curret merged segment
            if not connected_segment:
                merged_segments.append(merged_segment)
                break
            # Remove this segment from the segments list and add it to the polygon segments list
            remaining_segments.remove(connected_segment)
            merged_segment = merged_segment.combine_segment(connected_segment)
    return merged_segments

# Find groups of connected segments in a list of segments
# Yield any closed polygon when it is found
def connect_segments (segments : List[Segment]) -> Generator[Polygon, None, None]:
    remaining_segments = [ *segments ]
    while len(remaining_segments) > 0:
        first_segment = remaining_segments[0]
        remaining_segments.remove(first_segment)
        polygon_segments = [ first_segment ]
        while True:
            # Find a new segment which is connected to the last polygon segment
            # There must be always at least 1 segment connected
            last_segment = polygon_segments[-1]
            previous_segments = polygon_segments[0:-1]
            connected_segment = next(segment for segment in remaining_segments if segment.is_connected_with(last_segment))
            # Remove this segment from the segments list and add it to the polygon segments list
            remaining_segments.remove(connected_segment)
            polygon_segments.append(connected_segment)
            # Check if the new connected segment closes a polygon with any of the previous segments
            # Note that the closing segment may not be the first segment
            # This may happen in case we have 2 polygons with a same corner
            closing_segment = next((s for s, segment in enumerate(previous_segments) if segment.is_connected_with(connected_segment)), None)
            if closing_segment != None:
                polygon_segments = polygon_segments[closing_segment:]
                yield Polygon.non_canonical(polygon_segments)
                # Recover the discarded segments to the segments list in case the closing segment was not the first segment
                discarded_segments = polygon_segments[0:closing_segment]
                previous_segments += discarded_segments
                break

# Get the non overlap segments between a list of segments
# Segments must all be in the same line
def get_line_non_overlap_segments (segments : 'Segment') -> List['Segment']:
    # If all segments are not in the same line we cannot proceed
    if not all( current.same_line_as(nextone) for current, nextone in pairwise(segments) ):
        raise ValueError('Segments are not all in the same line')
    # Otherwise, order all segment points and check those regions between points where only 1 segment exists
    points = unique(sum([ list(segment.points) for segment in segments ], []))
    sorted_points = sort_points(points)
    # Track from which point only 1 segment is in
    first_point = None
    in_segments = { segment: False for segment in segments }
    non_overlap_segments = []
    for point in sorted_points:
        in_count = 0
        for segment in segments:
            if point in segment.points:
                in_segments[segment] = not in_segments[segment]
        if first_point:
            new_segment = Segment(first_point, point)
            non_overlap_segments.append(new_segment)
            first_point = None
        if sum(in_segments.values()) == 1:
            first_point = point
    # If segments do not overlap at any moment then return None
    return non_overlap_segments

# Get the non overlap segments between a list of segments
def get_non_overlap_segments (segments : 'Segment') -> List['Segment']:
    pool = [ *segments ]
    non_overlap_segments = []
    while len(pool) > 0:
        # Start from the first segment in the pool and find segments in the same line
        start_segment = pool[0]
        inline_segments = [ start_segment ] + [ other_segment for other_segment in pool[1:] if other_segment.same_line_as(start_segment) ]
        current_non_overlap_segments = get_line_non_overlap_segments(inline_segments)
        non_overlap_segments += current_non_overlap_segments
        # Remove current segments from the pool
        for segment in inline_segments:
            pool.remove(segment)
    return non_overlap_segments

# Given two catets, get the hipotenuse
def get_hipotenuse (segment_1 : Segment, segment_2 : Segment) -> Segment:
    if not segment_1.makes_corner_with(segment_2):
        raise ValueError('Input segment do not make a corner')
    unique_point_a = next(point for point in segment_1.points if point not in segment_2.points)
    unique_point_b = next(point for point in segment_2.points if point not in segment_1.points)
    return Segment(unique_point_a, unique_point_b)

# Get a random color for display tools
def generate_random_color ():
    r = round(random.random() * 255)
    g = round(random.random() * 255)
    b = round(random.random() * 255)
    return '#%02x%02x%02x' % (r,g,b)

# Create a cross with two segments to mark a specific point in the representation
# Return the segments
def mark_point (point : Point, color : str = 'black', size : number = 1):
    point_1a = point + Vector(-1, 1).normalized() * size
    point_1b = point + Vector(1, -1).normalized() * size
    segment_1 = Segment(point_1a, point_1b, color)
    point_2a = point + Vector(-1, -1).normalized() * size
    point_2b = point + Vector(1, 1).normalized() * size
    segment_2 = Segment(point_2a, point_2b, color)
    return segment_1, segment_2