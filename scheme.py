from typing import List, Union, Optional

from scheme_display import setup_display, add_frame, plot_lines

from vectorial_base import *
        
# A room is a smart perimeter that may contain other perimeters with conservative areas and size restrictions
# A start 'perimeter' may be passed. If no perimeters i passed it is assigned automatically according the room rules
# The 'forced_area' argument stablishes the expected final room area. If no area is passed then the original perimeter area will be used
# The 'min_size' and 'max_size' are the limits in both x and y axes
# The 'name' is only a representation parameter
class Room:
    def __init__(self,
        perimeter : Perimeter = None,
        forced_area : float = None,
        min_size : float = None,
        max_size : float = None,
        display : bool = False,
        name : str = 'Unnamed',
        children : list = [],
        ):
        
        self.perimeter = perimeter
        # Se the real area
        if perimeter:
            self.area = perimeter.area
        else:
            self.area = None
        # Set the expected final area
        if forced_area:
            self.forced_area = forced_area
        else:
            self.forced_area = self.area
        # Set the free area: area where there is no children rooms
        self.free_area = self.area
        # Set size limits
        self.min_size = min_size
        self.max_size = max_size
        # Set representation parameters
        self.display = display
        self.name = name
        # Set up the hierarchy of rooms
        self.parent = None
        self.children = []
        # Set up all children rooms
        self.add_children(children)
        # Update the representation after the setup
        self.update_display()

    # Add children rooms
    def add_children(self, rooms : list):
        if len(rooms) == 0:
            return
        # Check areas of all children rooms to do not sum up more than the parent area
        # In addition check if any of the children areas has perimeter and, if so, check the perimeter is inside the parent
        children_area = 0
        for room in rooms:
            children_area += room.forced_area
        if children_area > self.area:
            raise NameError('ERROR: Children require more area than the parent has')
        # Check all children are inside the perimeter, if they have a predefined perimeter
        for room in rooms:
            if room.perimeter and room.perimeter not in self.perimeter:
                raise NameError('ERROR: The child room "' + room.name + '" is out of the parent perimeter')
        # Set up each room by giving them a position and correct size to match the forced area
        for room in rooms:

            # Finally, update the room hierarchy
            room.parent = self
            self.children.append(room)

    # Get all lines from this room and all children room perimeters
    def get_lines_recuersive (self):
        lines = self.perimeter.lines if self.perimeter else []
        for room in self.children:
            lines += room.get_lines_recuersive()
        return lines

    # Add a new frame in the display with the current lines of this room and its children
    def update_display (self):
        if self.display:
            add_frame(self.get_lines_recuersive())