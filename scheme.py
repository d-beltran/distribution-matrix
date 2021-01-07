
from scheme_display import setup_display, add_frame, plot_lines

from functools import reduce

from typing import List, Union, Optional

# CLASS DEFINITIONS ------------------------------------------------------------

# An x,y coordinate
class Point:

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
    
    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y)
        
# A segment defined by 2 coordinates (Points): 'a' and 'b'
class Line:

    def __init__(self, a : Point, b : Point):
        self.a = a
        self.b = b
        self.vector = b - a

    def __str__(self):
        return 'A: ' + str(self.a) + ', B: ' + str(self.b)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.a == other.a and self.b == other.b
        return False

    def __hash__(self):
        return hash((self.a, self.b))

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
        if other.x and other.y:
            in_x = other.x > self.pmin.x and other.x < self.pmax.x
            in_y = other.y > self.pmin.y and other.y < self.pmax.y
            return in_x and in_y
        return False

    # Return all rectangle points in a 'perimeter-friendly' order
    def get_points(self):
        point1 = self.pmin
        point2 = Point(self.pmin.x, self.pmax.y)
        point3 = self.pmax
        point4 = Point(self.pmax.x, self.pmin.y)
        return [ point1, point2, point3, point4 ]

    # Return all rectangle lines in a 'perimeter-friendly' order
    def get_lines(self):
        points = self.get_points()
        lines = [ Line(a,b) for a, b in pairwise(points, retro=True) ]
        return lines

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

    

# A perimeter defined by several lines
class Perimeter:

    def __init__(self, lines : list):
        self.lines = lines
        self.closed_check() 
        self.set_corners()
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

    # A point between 2 lines in the perimeter
    # The 'inside' variable defines if the corner is "pointing" to the inside of the perimeter
    # (e.g. rectangles have no inside corners while an 'L' has one inside corner)
    class Corner:

        def __init__(self, point : Point, lines : list, inside : bool):
            self.x = point.x
            self.y = point.y
            self.lines = lines
            self.inside = inside

        def __repr__(self):
            return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

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
    def closed_check(self):
        if not self.is_closed():
            raise NameError('The perimeter is not closed')

    # Set the perimeter corners
    def set_corners(self):
        # For each line, save the end point and if the direction of the next line is left
        # Count how many corners are there in one direction (left in this case)
        precorners = []
        left_count = 0
        for current, nextone in pairwise(self.lines, retro=True):
            point = current.b
            lefted_corner = next_left(current, nextone)
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
            corner = self.Corner(point, lines, inside)
            corners.append(corner)
        self.corners = corners

    # Get all perimeter points
    # DANI: Obsoleto
    def get_points(self):
        self.closed_check()
        return [ line.a for line in self.lines ]

    # Split the current perimeter in a list of rectangle perimeters
    # If the multisplit option is true it will return as many splitted rectangles as possible
    # Else, it will return the minimum splits to cover all the perimeter
    # (e.g. an 'L' will return 2 rectangles if multisplit is false but 3 rectangles if it is true)
    def split_in_rectangles(self):
        self.closed_check()

        # Get all outside corners
        outside_corners = [ corner for corner in self.corners if not corner.inside ]

        # Create a new rectangle for each outside corner
        rects = [ Rect.from_lines(corner.lines) for corner in outside_corners ]

        # Discard rectangles which contain any interior corner
        # i.e. this rectange cuts the perimeter at some point
        inside_corners = [ corner for corner in self.corners if corner.inside ]
        def is_real(rect):
            for corner in inside_corners:
                if corner in rect:
                    return False
            return True
        real_rects = list(filter(is_real, rects))

        # Join all rectangles
        joined_rects = join_rects(real_rects)

        # Save this data since it is expensive
        self.rects = joined_rects

        # And finally return rectangles
        return joined_rects


    # Calculate the area of the current perimeter
    def get_area(self):
        self.closed_check()
        # If the area is previously calculated then return it
        if self._area:
            return self._area
        # Get all rectangles which form the perimeter
        rects = self.split_in_rectangles()
        # Get the area of each rectangle
        areas = [ rect.get_area() for rect in rects ]
        self._area = sum(areas)
        return self._area

    # The area is treated appart since it may be an expensive calculation
    area = property(get_area, None, None, "The perimeter area")
        

# Set the main class which may contain several perimeters defined as 'rooms'
class Scheme:
    def __init__(self, limits : Perimeter, display : bool = False):
        self.limits = limits
        self.display = display
        if display:
            self.update_display()
            setup_display()
        print('Scheme ready')

    # A room is a smart perimeter with conservative area and size restrictions
    # The 'area' argument stablishes the expected room area. If no area is passed then the original perimeter area will be used
    class Room:
        def __init__(self, perimeter : Perimeter, area = None, min_size = None, max_size = None):
            self.perimeter = perimeter
            if area:
                self.area = area
            else:
                self.area = perimeter.getArea()
    
    # Get all lines in the scheme
    # i.e. the limits and lines in all room perimeters
    def get_all_lines (self):
        return self.limits.lines

    # Add a new frame in the display with the current lines of the scheme
    def update_display (self):
        if self.display:
            add_frame(self.get_all_lines())


# Auxiliar functions ---------------------------------------------------------------

# Set a special iteration system
# Return each value of the array and the next value
# By default, the final array value is skipped, since it has no next value
# However, if the 'retro' argument is True, the final array value is returned with the first array value as the next value
def pairwise (values : list, retro : bool = False):
    last = len(values) - 1
    for i in range(0, last):
        yield values[i], values[i+1]
    if retro:
        yield values[last], values[0]

# Operations with lines ----------------------------------------------------------

# Given two continue lines which inevitably form a corner, set if the second line goes "left" (true) or "right" (false) relative to the first line
# https://stackoverflow.com/questions/13221873/determining-if-one-2d-vector-is-to-the-right-or-left-of-another  º
def next_left (line1 : Line, line2 : Line) -> bool:
    a = line1.vector
    b = line2.vector
    dot = a.x * -b.y + a.y * b.x
    return dot < 0

# Get the intersection between two lines
# https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines
def get_intersection_point (line1 : Line, line2 : Line) -> Optional[Point]:
    xdiff = (line1.a.x - line1.b.x, line2.a.x - line2.b.x)
    ydiff = (line1.a.y - line1.b.y, line2.a.y - line2.b.y)

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       return None

    d = (det(line1.a, line1.b), det(line2.a, line2.b))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return Point(x, y)

# Check if two lines intersect
def intersect (line1 : Line, line2 : Line) -> bool:
    xdiff = (line1.a.x - line1.b.x, line2.a.x - line2.b.x)
    ydiff = (line1.a.y - line1.b.y, line2.a.y - line2.b.y)

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       return False
    return True

# Operations with rectangles ----------------------------------------------------------

# Given 2 rectangles, it returns the overlapping region, if exists, as a new rectangle
def get_rects_overlap (rect1 : Rect, rect2 : Rect) -> Optional[Rect]:
    # Find the overlap in the 'x' dimension
    # Get the maximum of the minimums
    x_minimum = max(rect1.pmin.x, rect2.pmin.x)
    # Get the minimum of the maximums
    x_maximum = min(rect1.pmax.x, rect2.pmax.x)
    # Check that the overlap range exists
    x_overlap = ( x_minimum, x_maximum )
    if x_minimum >= x_maximum:
        return None
    # Find the overlap in the 'y' dimension
    # Get the maximum of the minimums
    y_minimum = max(rect1.pmin.y, rect2.pmin.y)
    # Get the minimum of the maximums
    y_maximum = min(rect1.pmax.y, rect2.pmax.y)
    # Check that the overlap range exists
    y_overlap = ( y_minimum, y_maximum )
    if y_minimum >= y_maximum:
        return None
    # Build a rectangle with both dimensional overlaps
    pmin_overlap = Point(x_overlap[0], y_overlap[0])
    pmax_overlap = Point(x_overlap[1], y_overlap[1])
    overlap = Rect(pmin_overlap, pmax_overlap)
    return overlap

def split_rect (rect : Rect, x_splits : list = [], y_splits : list = []):
    # Get the rectangle minimum and maximum values
    pmin = rect.pmin
    pmax = rect.pmax
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


# Join 2 rectangles by returning a list of all rectangles which define the resulting polygon
# The overlapped region, if exists is transformed to a single rectangle
def join_2_rects (rect1 : Rect, rect2 : Rect) -> Optional[list]:
    # Find the overlap between these two rectangles
    overlap = get_rects_overlap(rect1, rect2)
    # If there is no overlap then just return both input rectangles
    if not overlap:
        return [rect1, rect2]
    # Split both input rectangles using the overlap maximum and minimum points as split points
    x_splits = [ overlap.pmin.x, overlap.pmax.x ]
    y_splits = [ overlap.pmin.y, overlap.pmax.y ]
    split1 = split_rect(rect1, x_splits, y_splits)
    split2 = split_rect(rect2, x_splits, y_splits)
    rects = list(split1) + list(split2)
    # Return only unique rectangles
    return set(rects)

# Join multiple rectangles by joining all of them by pairs
def join_rects (rects : list) -> Optional[list]:
    joined_rects = []
    for rect_1, rect_2 in pairwise(rects, retro=True):
        joined_rects += join_2_rects(rect_1, rect_2)
    # Return only unique rectangles
    return set(joined_rects)

# Substract the second rectangle form the first rectangle and return the rectangles which define the resulting polygon
def subtract_rects (rect1 : Rect, rect2 : Rect) -> Optional[list]:
    # Find the overlap between these two rectangles
    overlap = get_rects_overlap(rect1, rect2)
    # If there is no overlap then just return the first rectangle intact
    if not overlap:
        return [rect1]
    # Split the first input rectangle using the overlap maximum and minimum points as split points
    x_splits = [ overlap.pmin.x, overlap.pmax.x ]
    y_splits = [ overlap.pmin.y, overlap.pmax.y ]
    split = split_rect(rect1, x_splits, y_splits)
    # Return all splits but the overlapped rectangle
    rects = [ rect for rect in list(split) if rect != overlap ]
    return rects