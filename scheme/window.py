from math import sqrt, inf, tan, atan, degrees, radians
from typing import List, Tuple, Dict, Union, Optional

from utils.auxiliar import *
from vectorial_base import *

# A window is a segment in a boundary
# When boundaries are transformed to walls with tickness, windows become holes in the wall
# Windows belong to the root room (e.g. a building floor) and they are placed in its boundary (i.e. the facade)
# Windows are set before the solving process and they never move, so children rooms must adapt to them
class Window:
    def __init__ (self,
        # The point where the window is centered
        # It must be in the room boundary
        point : Optional[Point] = None,
        # Set how wide the window must be
        width : Optional[number] = None,
        # Set the minimum width of margins on each side of the window
        margin : Optional[number] = None,
        # Set the parent room
        room : Optional['Room'] = None,
        # Set a name for the window
        # This is a representation parameters and it has no effect in the logic
        name : Optional[str] = None,
    ):
        # Save input values as internal values
        # Width and margin are usually None at this point
        # They are usually set further from the window room 'window_args' value
        self.name = name
        self._width = width
        self._margin = margin
        self._point = point
        self._margined_width = None
        self._segment = None
        self._margined_segment = None
        # The room this window belongs to
        self.room = room
        # Store the direction from the window towards the "outside" of the room
        self._outside_direction = None

    def __repr__ (self):
        name = self.name if self.name else 'Unnamed'
        point = f'placed in {self._point}' if self._point else '(not placed)'
        width = f'with a width of {self._width}' if self._width else '(widthless)'
        margin = f'and with a margin of {self._margin}' if self._margin else '(marginless)'
        return f'<Window "{name}" {point} {width} {margin}>'

    # Get the width
    def get_width (self) -> number:
        # If we have a stored value already then return it
        if self._width != None:
            return self._width
        # Otherwise we must get it from the parent room args
        # If there is no parent room then we have nothing to do
        if not self.room:
            return None
        self._width = self.room.window_args['width']
        return self._width

    # Set the width (regular setter)
    def set_width (self, new_width : number):
        self._width = new_width

    # The window width
    width = property(get_width, set_width, None, "The window width")

    # Get the margin
    def get_margin (self) -> number:
        # If we have a stored value already then return it
        if self._margin != None:
            return self._margin
        # Otherwise we must get it from the parent room args
        # If there is no parent room then we have nothing to do
        if not self.room:
            return None
        self._margin = self.room.window_args['margin']
        return self._margin

    # Set the margin (regular setter)
    def set_margin (self, new_margin : number):
        self._margin = new_margin

    # The window margin
    margin = property(get_margin, set_margin, None, "The window margin")

    # Get the margined width
    def get_margined_width (self):
        if self._margined_width != None:
            return self._margined_width
        if self.width == None:
            raise ValueError('Window is missing width')
        if self.margin == None:
            raise ValueError('Window is missing margin')
        self._margined_width = self.width + self.margin * 2
        return self._margined_width

    # The window margined width
    margined_width = property(get_margined_width, None, None, "The window margined width")

    # Get the window point
    def get_point (self) -> Optional[Point]:
        return self._point

    # If the window point is set then reset its segment and margined segment
    def set_point (self, point : Optional[Point]):
        self._point = point
        self._segment = None
        self._margined_segment = None
        self._outside_direction = None
        if point:
            self.segment
            self.margined_segment

    # The window point
    point = property(get_point, set_point, None, "The window point")

    # Get the window segment
    def get_segment (self) -> Optional[Segment]:
        # Return internal value if it exists
        if self._segment:
            return self._segment
        self._segment = self.generate_segment(self.width)
        return self._segment

    # The window segment
    segment = property(get_segment, None, None, "The window segment")

    # Get the window margined segment
    def get_margined_segment (self) -> Optional[Segment]:
        # Return internal value if it exists
        if self._margined_segment:
            return self._margined_segment
        # If the window has no point then complain
        if not self.point:
            raise RuntimeError('Trying to get window margined segment when no point is defined')
        self._margined_segment = self.generate_segment(self.margined_width)
        return self._margined_segment

    # The window margined segment
    margined_segment = property(get_margined_segment, None, None, "The window margined segment")

    # Given a segment width, generate a new
    # The new segment will be centered in the window point
    # The new segment will be overlaped with the boundary segment where the window point is
    def generate_segment (self, width : number) -> Segment:
        # If width is 0 then the segment can not exist
        if width == 0:
            raise ValueError('Cannot generate a segment for a window of width 0')
        # If the window point is not assigned then we can not generate the segment
        if not self.point:
            return None
        # If we can not retrieve the boundary then we can not generate the segment
        room_boundary = self.get_room_boundary()
        if not room_boundary:
            return None
        # Otheriwse, generate the segment
        boundary_segment = next(( segment for segment in room_boundary.segments if self.point in segment ), None)
        if not boundary_segment:
            raise ValueError(f'The window point {self.point} is not over its room boundary ({self.room.name})')
        direction = boundary_segment.direction
        half_width = width / 2
        a = self.point - direction * half_width
        b = self.point + direction * half_width
        if a not in boundary_segment or b not in boundary_segment:
            raise ValueError(f'The window segment ({Segment(a,b)}) does not fit in its room boundary ({self.room.name})')
        return Segment(a,b)

    # Get the window outside direction
    # i.e. the direction from the window towards the outside side of its room boundary
    # This is the direction a window looks through
    def get_outside_direction (self) -> Optional[Vector]:
        # Return internal value if it exists
        if self._outside_direction:
            return self._outside_direction
        if not self.segment:
            return None
        room_boundary = self.get_room_boundary()
        if not room_boundary:
            return None
        self._outside_direction = -room_boundary.get_border_inside(self.segment)
        return self._outside_direction

    # The window outside direction
    # It crosses the window segment perpendicularly and it is a normalized vector
    outside_direction = property(get_outside_direction, None, None, "The window outside direction")

    # Get the boundary where the window is meant to be
    def get_room_boundary (self) -> Optional[Boundary]:
        room = self.room
        if not room:
            return None
        boundary = room.boundary
        if not boundary:
            return None
        return boundary

    # Make a copy of this window
    def copy(self) -> 'Window':
        return Window(
            point = self.point,
            width = self._width,
            margin = self._margin,
            room = self.room,
            name = self.name,
        )
