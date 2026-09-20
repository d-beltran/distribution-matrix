from math import sqrt, inf, tan, atan, degrees, radians
from typing import List, Tuple, Dict, Union, Optional

from utils.auxiliar import *
from vectorial_base import *

# A window is a segment in a boundary
# When boundaries are transformed to walls with tickness, windows become holes in the wall
# Windows are located in external walls
class Window:
    def __init__ (self,
        # A point may be passed. If no point is passed then it is assigned automatically
        # WARNING: This is only supported if the door room boundary is already set and the point matches on it
        point : Optional[Point] = None,
        # Set how wide the door must be
        width : Optional[number] = None,
        # Set the minimum width of margins on each side of the door
        margin : Optional[number] = None,
        # Set if the door is rigid
        # i.e. its point may not change as a result of the solving process
        rigid : bool = False,
        # Set the parent room
        room : Optional['Room'] = None,
    ):
        # Save input values as internal values
        # These values are usually None at this point
        # They are usually set further from the door room 'door_args' value
        self._width = width
        self._margin = margin
        self._point = point
        self._margined_width = None
        self._segment = None
        self._margined_segment = None
        self.rigid = rigid
        if self.rigid and not self.point:
            raise InputError('A point must be defined if the rigid flag is passed')
        self.reverse = reverse
        # The room this door belongs to
        self.room = room