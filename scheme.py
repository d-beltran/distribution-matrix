from typing import List, Union, Optional

from scheme_display import add_frame, plot_lines

from vectorial_base import *

import random  
        
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
        # Set internal variables
        self._free_rects = None
        self._free_mrects = None
        self._perimeter = None
        # Set representation parameters
        self.display = display
        self.name = name
        # Set up the hierarchy of rooms
        self.parent = None
        self.children = []
        # Set the perimeter
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
        if max_size:
            self.max_size = max_size
        elif self.forced_area and self.min_size:
            self.max_size = self.forced_area / self.min_size
        else:
            self.max_size = None
        # Set up all children rooms
        self.add_children(children)
        # Update the representation after the setup
        self.update_display()

    # Get the perimeter
    # Just return the internal perimeter value
    def get_perimeter (self):
        return self._perimeter

    # Set the perimeter
    # Reset the own rects and the parent rects also
    def set_perimeter (self, value):
        self.reset_rects()
        if self.parent:
            self.parent.reset_rects()
        self._perimeter = value

    # The area is treated appart since it may be an expensive calculation
    perimeter = property(get_perimeter, set_perimeter, None, "The room perimeter")

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
            return self.perimeter.mrects
        # Split in rectangles using the children as exclusion perimeters
        free_mrects = self.perimeter.get_maximum_rectangles( splitted_rects = self.free_rects )
        self._free_mrects = free_mrects
        return free_mrects

    # The area is treated appart since it may be an expensive calculation
    free_mrects = property(get_free_mrects, None, None, "The maximum free are rectnagles")

    # Reset all minimum and maximum free rects
    # This must be done each time the perimeter is modified since they are not valid anymore
    def reset_rects(self):
        self._free_rects = None
        self._free_mrects = None

    # Check if a rectangle fits somewhere in the perimeter
    # Iterate over all maximum rectangles and try to fit
    # If only the x size parameter is passed it is assumed to be both x and y size
    def fit (self, x_fit_size : float, y_fit_size : float = None):
        if not y_fit_size:
            y_fit_size = x_fit_size
        for rect in self.free_mrects:
            x_size, y_size = rect.get_size()
            if x_fit_size <= x_size and y_fit_size <= y_size:
                return True
        return False

    # Get all maximum free rectangles with the minimum specified x and y sizes
    # Iterate over all free maximum rectangles and try to fit
    # If only the x size parameter is passed it is assumed to be both x and y size
    def get_fit (self, x_fit_size : float, y_fit_size : float = None):
        fit_rects = []
        if not y_fit_size:
            y_fit_size = x_fit_size
        for rect in self.free_mrects:
            x_size, y_size = rect.get_size()
            if x_fit_size <= x_size and y_fit_size <= y_size:
                fit_rects.append(rect)
        return fit_rects

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

        # Sort children rooms by minimum size, with the biggest sizes first
        def sort_by_size(room):
            return room.min_size
        sorted_rooms = sorted( rooms, key=sort_by_size, reverse=True )

        # Set up each room by giving them a position and correct size to match the forced area
        for room in sorted_rooms:
            # Update the room hierarchy
            room.parent = self
            self.children.append(room)
            # If the children has no perimeter it must be built
            if not room.perimeter:
                self.set_child_room_perimeter(room)

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
        # The minimum base perimeter is a square with both sides as long as the room minimum size
        suitable_rects = self.get_fit(room.min_size)
        # Shuffle the suitable rects
        random.shuffle(suitable_rects)
        # Sort the suitable rects by minimum size, with the biggest sizes first
        def sort_by_size(rect):
            return min(rect.get_size())
        sorted_suitable_rects = sorted( suitable_rects, key=sort_by_size, reverse=True )
        # Stop here if there are no available rectangles
        # DANI: Esto no debería pasar nunca. Debería preveerse de antes.
        if len(suitable_rects) == 0:
            raise NameError('ERROR: The room ' + room.name + ' fits nowhere')
        # Try to set up the new room in all possible sites until we find one
        # Each 'site' means each corner in each suitable rectnagle
        # Check each site to allow other rooms to grow
        for rect in sorted_suitable_rects:
            for corner in rect.get_corners():
                minimum_rect = Rect.from_corner(corner, room.min_size, room.min_size)
                room.perimeter = Perimeter(minimum_rect.get_lines())
                self.update_display()
        # 

    # Set the initial room perimeter as the maximum possible rectangle
    # This is a shortcut to skip the difficult expansion protocol
    # It is useful to set a whole room at the begining, when there is plenty of free space
    # It is useful to set as much perimeter as possible at the start before we resolve space conflicts
    def set_maximum_initial_perimeter(self, corner : Point, directions : list, space : Rect):
        x_space, y_space = space.get_size()
        # If the room area is greater than the space then return the whole space as a permeter
        if space.area <= self.area:
            # If both space sizes are bigger than the maximum size we can not return the whole space
            # The size in one of both sides must be limited (the biggest size)
            if x_space > self.max_size and y_space > self.max_size:
                if x_space >= y_space:
                    # DANI: Me quedé aquí
                    pass

            return Perimeter(space.get_lines())
        # Else, fit the room in the space
        # DANI: Me quedé aquí
        
        #line_1 = Line(corner, corner - directions[0] * self.min_size)
        #line_2 = Line(corner, corner + directions[1] * self.min_size)
        #minimum_rect = Rect.from_lines([ line_1, line_2 ])
        #return Perimeter(minimum_rect.get_lines())

    # Go uppwards in the hyerarchy until you reach the room which has no parent
    def get_root_room(self):
        root = self
        while(root.parent):
            root = root.parent
        return root

    # Get all lines from this room and all children room perimeters
    def get_lines_recuersive (self, only_children : bool = False):
        lines = []
        if not only_children and self.perimeter:
            lines += self.perimeter.lines
        for room in self.children:
            lines += room.get_lines_recuersive()
        return lines

    # Add a new frame in the display with the current lines of this room and its children
    def update_display (self):
        # Find the root room
        root = self.get_root_room()
        if root.display:
            add_frame(root.get_lines_recuersive())