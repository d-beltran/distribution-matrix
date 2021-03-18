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
        # If the forced area does not cover the minimum size it makes no sense
        if min_size and forced_area < min_size**2:
            raise NameError('Input error: Forced area is not sufficient for the minimum size in room ' + name)
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
        # Set internal variables
        self._free_rects = None
        self._free_mrects = None

    # Get the available space inside de perimeter splitted in rects
    # i.e. space not filled by children rooms
    def get_free_rects (self):
        # If rects are previously calculated then return them
        if self._free_rects:
            return self._free_rects
        # Return none if there is not perimeter yet
        if not self.perimeter:
            return None
        # If there are no children then return the current perimeter rectangles
        if len(self.children) == 0:
            return self.perimeter.rects
        # Split in rectangles using the children as exclusion perimeters
        free_rects = self.perimeter.split_in_rectangles( exclusion_perimeters = [ child.perimeter for child in self.children if child.perimeter ] )
        self._free_rects = free_rects
        return free_rects

    # The area is treated appart since it may be an expensive calculation
    free_rects = property(get_free_rects, None, None, "The room free area splitted in rectangles")

    # Get the maximum joined rectangles from the free rectangles
    def get_free_mrects (self):
        # If rects are previously calculated then return them
        if self._free_mrects:
            return self._free_mrects
        # Return none if there is not perimeter yet
        if not self.perimeter:
            return None
        # If there are no children then return the current perimeter rectangles
        if len(self.children) == 0:
            return self.mrects
        # Split in rectangles using the children as exclusion perimeters
        free_mrects = self.perimeter.get_maximum_rectangles( splitted_rects = self.free_rects )
        self._free_mrects = free_mrects
        return free_mrects

    # The area is treated appart since it may be an expensive calculation
    free_mrects = property(get_free_mrects, None, None, "The maximum free are rectnagles")

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
            raise NameError('Input error: Children together require more area than the parent has')
        # Check all children are inside the perimeter, if they have a predefined perimeter
        for room in rooms:
            if room.perimeter and room.perimeter not in self.perimeter:
                raise NameError('Input error: The child room "' + room.name + '" is out of the parent perimeter')
        # Check all children minim sizes fit in the perimeter, if they have a predefined minimum size
        for room in rooms:
            size = room.min_size
            if size and not self.perimeter.fit(size, size):
                raise NameError('Input error: The child room "' + room.name + '" minimum size does not fit in the parent perimeter')    
        # Set up each room by giving them a position and correct size to match the forced area
        for room in rooms:
            # If the children has no perimeter it must be built
            if not room.perimeter:
                self.set_child_room_perimeter(room)
            # Finally, update the room hierarchy
            room.parent = self
            self.children.append(room)

    # Get all lines from this room and all children room perimeters
    def get_lines_recuersive (self, only_children : bool = False):
        lines = []
        if not only_children and self.perimeter:
            lines += self.perimeter.lines
        for room in self.children:
            lines += room.get_lines_recuersive()
        return lines

    # Check if a point is in the border of any children perimeter
    def in_children_border(self, point : Point):
        for room in self.children:
            for corner in room.perimeter.corners:
                if point == corner:
                    return True
            for line in room.perimeter.lines:
                if point in line:
                    return True
        return False

    # Set up a room perimeter
    def set_child_room_perimeter(self, room):
        # Find a suitable maximum free rectangle to deploy a starting base perimeter
        # The minimum base perimeter is a square which both sides as long as the room minimum size
        pass

    # Add a new frame in the display with the current lines of this room and its children
    def update_display (self):
        if self.display:
            add_frame(self.get_lines_recuersive())