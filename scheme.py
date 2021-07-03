from typing import List, Union, Optional

from scheme_display import add_frame

from vectorial_base import *

from utils import InputError

import random
from math import sqrt
        
# A room is a smart perimeter that may contain other perimeters with conservative areas and size restrictions
# A start 'perimeter' may be passed. If no perimeters i passed it is assigned automatically according the room rules
# The 'forced_area' argument stablishes the expected final room area. If no area is passed then the original perimeter area will be used
# Forced area may be a number (absolute value) or a string (percent) with the format 'XX%'
# The 'min_size' and 'max_size' are the limits in both x and y axes
# The 'name' is only a representation parameter
class Room:
    def __init__ (self,
        perimeter : Optional[Perimeter] = None,
        forced_area : Optional[Union[number, str]] = None,
        min_size : Optional[number] = None,
        max_size : Optional[number] = None,
        display : bool = False,
        name : str = 'Unnamed',
        segments_color : str = 'black',
        fill_color : str = 'white',
        children : List['Room'] = [],
        ):
        # Set internal variables
        self._perimeter = None
        self._free_grid = None
        # Set representation parameters
        self.display = display
        self.name = name
        self.segments_color = segments_color
        self.fill_color = fill_color
        # Set up the hierarchy of rooms
        # Parent is never assigned from the instance itself, but it is assigned by the parent
        self.parent = None
        self.children = []
        # Set the perimeter
        # If the perimeter has been forced then update the display with the initial segments
        if perimeter:
            self.perimeter = perimeter
            self.update_display()
        else:
            self.perimeter = None
        # Set the expected final area
        self._forced_area_portion = None
        if forced_area:
            # If it is a number
            if is_number(forced_area):
                self.forced_area = forced_area
            # If it is a percent
            elif type(forced_area) == str and forced_area[-1] == '%':
                self.forced_area = None
                self._forced_area_portion = float(forced_area[0:-1]) / 100
            else:
                raise InputError('Forced area has a non-supported format')
        else:
            self.forced_area = self.area
        # If the forced area does not cover the minimum size it makes no sense
        if min_size and self.forced_area and self.forced_area < min_size**2:
            raise InputError('Forced area is not sufficient for the minimum size in room ' + name)
        # Set size limits
        if min_size or min_size == 0:
            self.min_size = min_size
        else:
            self.min_size = 0
        self._forced_max_size = max_size
        # Set up all children rooms
        self.add_children(children)

    # Get the perimeter
    # Just return the internal perimeter value
    def get_perimeter (self):
        return self._perimeter

    # Set the perimeter
    # Reset the own rects and the parent rects also
    def set_perimeter (self, perimeter):
        self._perimeter = perimeter
        self.reset_free_grid()
        if self.parent:
            self.parent.reset_free_grid()

    # The room perimeter
    perimeter = property(get_perimeter, set_perimeter, None, "The room perimeter")

    # Area inside the room perimeter (read only)
    def get_area(self):
        return self.perimeter.area
    area = property(get_area, None, None, "Area inside the room perimeter")

    # Get the available space inside the room perimeter as a rectangles grid
    # i.e. space not filled by children rooms
    def get_free_grid (self):
        # If rects are previously calculated then return them
        if self._free_grid:
            return self._free_grid
        # Return none if there is not perimeter yet
        if not self.perimeter:
            return None
        # If there are no children then return the current perimeter rectangles
        # If all children have no perimeter then return the current perimeter rectangles
        if len(self.children) == 0 or not any([ child.perimeter for child in self.children ]):
            free_grid = self.perimeter.grid
        # Otherwise, split in rectangles using the children as exclusion perimeters
        else:
            children_perimeters = [ child.perimeter for child in self.children if child.perimeter ]
            free_rects = self.perimeter.split_in_rectangles( exclusion_perimeters = children_perimeters )
            free_grid = Grid(free_rects)
        # Apply the current room colors to all rectangles
        for rect in free_grid.rects:
            rect.segments_color = self.segments_color
            rect.fill_color = self.fill_color
        self._free_grid = free_grid
        # ---------------------------------------------------------------------------------------------------
        # DANI: Muestra los rects la primera vez que se calculan
        colored_rects = [ rect.get_colored_rect(segments_color='red') for rect in free_grid.rects ]
        self.update_display(colored_rects)
        # ---------------------------------------------------------------------------------------------------
        return free_grid
    # Free space grid (read only)
    free_grid = property(get_free_grid, None, None, "The room free space grid")

    # Get rectangles in the free grid
    def get_free_rects (self):
        return self.free_grid.rects
    free_rects = property(get_free_rects, None, None, "The room free area splitted in rectangles")

    # Get the maximum rectangles from the free grid
    def get_free_max_rects (self):
        # ---------------------------------------------------------------------------------------------------
        # DANI: Muestra los max rects la primera vez que se calculan
        first_time = self.free_grid._max_rects == None
        if first_time:
            colored_rects = [ rect.get_colored_rect(segments_color='blue') for rect in self.free_grid.max_rects ]
            #self.update_display(colored_rects)
        # ---------------------------------------------------------------------------------------------------
        return self.free_grid.max_rects
    free_max_rects = property(get_free_max_rects, None, None, "The free maximum rectangles")

    # Free space area (read only)
    def get_free_area (self):
        return self.free_grid.area
    free_area = property(get_free_area, None, None, "Free space area (read only)")

    # Get the maximum size
    # If there is a forced maximum size then return it
    # Otherwise calculate it by dividing the forced area between the minimum size
    # i.e. return the maximum possible size in case the perimeter was the thinest rectangle
    def get_maximum_area (self):
        if self._forced_max_size:
            return self._forced_max_size
        elif self.forced_area and self.min_size and self.min_size != 0:
            return self.forced_area / self.min_size
        else:
            return None
    # Force the maximum size
    def set_maximum_size (self, value):
        self._forced_max_size = value
    max_size = property(get_maximum_area, set_maximum_size, None, "Maximum possible size")

    # Reset all minimum and maximum free rects
    # This must be done each time the perimeter is modified since they are not valid anymore
    def reset_free_grid (self):
        self._free_grid = None

    # Check if a rectangle fits somewhere in the perimeter
    # Iterate over all maximum rectangles and try to fit
    # If only the x size parameter is passed it is assumed to be both x and y size
    def fit (self, x_fit_size : number, y_fit_size : number = None):
        if not y_fit_size:
            y_fit_size = x_fit_size
        for rect in self.free_max_rects:
            x_size, y_size = rect.get_size()
            if x_fit_size <= x_size and y_fit_size <= y_size:
                return True
        return False

    # Get all maximum free rectangles with the minimum specified x and y sizes
    # Iterate over all free maximum rectangles and try to fit
    # If only the x size parameter is passed it is assumed to be both x and y size
    def get_fit (self, x_fit_size : number, y_fit_size : number = None):
        fit_rects = []
        if not y_fit_size:
            y_fit_size = x_fit_size
        for rect in self.free_max_rects:
            x_size, y_size = rect.get_size()
            if x_fit_size <= x_size and y_fit_size <= y_size:
                fit_rects.append(rect)
        return fit_rects

    # Add children rooms
    def add_children (self, rooms : list):
        if len(rooms) == 0:
            return
        # Check areas of all children rooms to do not sum up more than the parent area
        # In addition check if any of the children areas has perimeter and, if so, check the perimeter is inside the parent
        children_area = 0
        for room in rooms:
            # If the children has a percent forced area this is the time to calculate the absolute forced area
            if room._forced_area_portion:
                room.forced_area = self.forced_area * room._forced_area_portion
            children_area += room.forced_area
        if children_area > self.area:
            raise InputError('Children together require more area than the parent has')
        # Check all children are inside the perimeter, if they have a predefined perimeter
        for room in rooms:
            if room.perimeter and room.perimeter not in self.perimeter:
                raise InputError('The child room "' + room.name + '" is out of the parent perimeter')
        # Check all children minim sizes fit in the perimeter, if they have a predefined minimum size
        for room in rooms:
            size = room.min_size
            if size and not self.perimeter.fit(size, size):
                raise InputError('The child room "' + room.name + '" minimum size does not fit in the parent perimeter')

        # Sort children rooms by minimum size, with the biggest sizes first
        def sort_by_size (room):
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
    def in_children_border (self, point : Point):
        for room in self.children:
            for corner in room.perimeter.corners:
                if point == corner:
                    return True
            for segment in room.perimeter.segments:
                if point in segment:
                    return True
        return False

    # Set up a room perimeter
    def set_child_room_perimeter (self, room) -> bool:
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
            raise RuntimeError('ERROR: The room ' + room.name + ' fits nowhere')
        # Try to set up the new room in all possible sites until we find one
        # Each 'site' means each corner in each suitable rectangle
        # Check each site to allow other rooms to grow
        def set_base_perimeter():
            for rect in sorted_suitable_rects:
                for corner in rect.get_corners():
                    initial_perimeter = room.set_maximum_initial_perimeter(corner, rect)
                    if initial_perimeter:
                        return initial_perimeter
            return None
        base_perimeter = set_base_perimeter()
        # Set the child first perimeter, which automatically resets this room free grid
        room.perimeter = base_perimeter
        self.update_display()

        # Proceed with the expansion of this child room until it reaches its forced area
        if not room.expand():
            return False
        return True

    # Set the initial room perimeter as the maximum possible rectangle
    # This is a shortcut to skip the difficult expansion protocol
    # It is useful to set a whole room at the begining, when there is plenty of free space
    # It is useful to set as much perimeter as possible at the start before we resolve space conflicts
    def set_maximum_initial_perimeter (self, corner : Point, space : Rect) -> Perimeter:
        x_space, y_space = space.get_size()
        # If the room area is greater than the space then return the whole space as a permeter
        if space.area <= self.forced_area:
            # If both space sizes are bigger than the maximum size we can not return the whole space
            # The size in one of both sides must be limited (the biggest size is filled)
            if x_space > self.max_size and y_space > self.max_size:
                if x_space >= y_space:
                    maximum_rect = Rect.from_corner(corner, x_space, self.max_size)
                else:
                    maximum_rect = Rect.from_corner(corner, self.max_size, y_space)
                return Perimeter(maximum_rect.get_segments())
            return Perimeter(space.get_segments())
        # Else, fit the room in the space
        # First of all create the 3 rule rectangle from the space
        # i.e. calculate the relation of areas and apply it to the square root to both x and y sizes
        area_relation = self.forced_area / space.area # At this point the relation is always < 1
        new_x_size = x_space * sqrt(area_relation)
        new_y_size = y_space * sqrt(area_relation)
        # If any of the new sizes is shorter than the maximum size then the rectangle is valid
        if new_x_size < self.max_size or new_y_size < self.max_size:
            maximum_rect = Rect.from_corner(corner, new_x_size, new_y_size)
            return Perimeter(maximum_rect.get_segments())
        # If both new sizes are longer than the maximum size we must find another solution
        # The new rectangle will have the maximum size in one of its sides
        # Calculate the other side size
        second_size = self.forced_area / self.max_size
        # Now set which are the maximum and minimum sizes of the new rect
        new_sizes = [ self.max_size, second_size ]
        new_max_size = max(new_sizes)
        new_min_size = min(new_sizes)
        # Create the new rect fitting the biggest size in the biggest space and the opposite
        # For each size of the new rect use the calculated size only if it is not longer than the space
        if x_space >= y_space:
            new_x_size = min(new_max_size, new_x_size)
            new_y_size = min(new_min_size, new_y_size)
        else:
            new_x_size = min(new_min_size, new_x_size)
            new_y_size = min(new_max_size, new_y_size)
        maximum_rect = Rect.from_corner(corner, new_x_size, new_y_size)
        return Perimeter(maximum_rect.get_segments())

    # Calculate how much area we need to expand
    def get_required_area (self):
        return resolute(self.forced_area - self.area, base_resolution - 2)

    # Expand this room until it reaches the forced area
    def expand (self) -> bool:
        while self.get_required_area() > 0:
            if not self.expand_step():
                return False
        return True

    # Find the most suitable space for the current room to claim
    # If current room is expanded over another room then the other room must also expand to compensate
    # At the end, the extra space will be substracted from the parent free space
    def expand_step (self) -> bool:

        # Calculate how much area we need to expand
        required_area = self.get_required_area()
        print('Forced area: ' + str(self.forced_area))
        print('Current area: ' + str(self.area))
        print('Required area: ' + str(required_area))
        if required_area <= 0:
            return

        # Get the parent room
        parent_room = self.parent

        # ----------------------------------------------------------------------------------------------------
        # Try to expand in the easiest way: push an individual frontier segment
        # Only segment longer than the minimum size may be expanded this way
        # Start pushing free frontiers
        # - No matter if the space is divided as a result of the expansion
        # If there are no suitable free frontiers then proceed with brother frontiers
        # - Space must never be divided as a result of the expansion
        # ----------------------------------------------------------------------------------------------------

        def expand_frontier (frontier : Segment) -> bool:
            print("Expanding '" + self.name + "' -> " + str(frontier))
            # Set some parameters according to the room this frontier belongs to
            room = frontier.room
            grid = room.perimeter.grid if room else parent_room.free_grid
            margin_limit = room.min_size if room else parent_room.min_size
            # Find the maximum rectangles which are in contact with our frontier
            max_rects = grid.max_rects
            contact_max_rects = [ rect for rect in max_rects if frontier in rect ]
            # In case segment is vertical:
            # - The expansion uses maximum columns
            # - The expansion 'froward' is the x dimension
            # - The expansion 'sides' is the y dimension
            # In case segment is horizontal
            # - The expansion uses maximum rows
            # - The expansion 'froward' is the y dimension
            # - The expansion 'sides' is the x dimension
            # Here 'forward' and 'sides' may mean '0 -> x dimension' or '1 -> y dimension'
            # This is because rectangles size is given in a (x,y) tuple format
            # So size[0] = x and size[1] = 1
            if frontier.is_vertical():
                rects = grid.find_column_rectangles()
                forward = 0
                sides = 1
            elif frontier.is_horizontal():
                rects = grid.find_row_rectangles()
                forward = 1
                sides = 0
            else:
                raise RuntimeError('ERROR: diagonal segments are not supported for room expansion')

            if len(rects) == 0:
                raise RuntimeError('No pot ser: ' + str(grid))

            # One and only one of the rows/columns will always include the segment
            space = next((rect for rect in rects if frontier in rect), None)
            if space == None:
                # If this happens it may mean there is a problem with the grid
                raise RuntimeError('Frontier ' + str(frontier) + ' has no space in ' + str(rects))
            space_contact = next(segment for segment in space.segments if frontier in segment)
            space_forward_limit = space.get_size()[forward]

            # Find the direction of the expansion as a vector
            forward_direction = space_contact.get_middle_point() + space.get_middle_point()

            # Get the forward expansion limit according to maximum rectangles
            maximum_forward_limit = max([ rect.get_size()[forward] for rect in contact_max_rects ])
            margined_forward_limit = maximum_forward_limit - margin_limit

            # Now create a function to make a rectangle by pushing a segment (which may change)
            # The forward length of this rectangle will depend on the available area and limits
            # This segment may not be the original frontier, but a variation of it
            def push_segment (pushed_segment : Segment) -> bool:
                pushed_segment_size = pushed_segment.length
                push_length = required_area / pushed_segment_size
                # If the length exceeds the maximum limit then stay in the maximum limit
                if push_length > maximum_forward_limit:
                    push_length = maximum_forward_limit
                # If the length is between the maximum limit and the margin limit then stay at the margin
                elif push_length > margined_forward_limit:
                    push_length = margined_forward_limit
                # If at this point the length exceeds the space limit then stay in the space limit
                if push_length > space_forward_limit:
                    push_length = space_forward_limit
                # If the push length at this point is 0 or close to it then we can not push
                if push_length < 0.0001:
                    print('WARNING: The push length is too small')
                    return False
                # Create the new rect with the definitive length
                new_point = pushed_segment.a + forward_direction.normalized() * push_length
                #print('NEW POINT: ' + str(new_point))
                new_side = Segment(pushed_segment.a, new_point)
                new_rect = Rect.from_segments([pushed_segment, new_side])
                # Then add the new current rectangle by modifying the current perimeter
                # Substract the pushed segment from the perimeter and add the new created segments
                # First find the new created segments as a result of the expansion
                new_segments = [ seg for seg in new_rect.segments if seg != pushed_segment ]
                # In addition, it may happend that the new segments overlap with current segments
                # e.g. when the expanded frontier makes a corner with another segment from this room
                # In those cases we must substract the overlapped part from both new and current segments
                overlap_segments = []
                substract_segments = []
                for current_segment in self.perimeter.segments:
                    for new_segment in new_segments:
                        overlap = new_segment.get_overlap_segment(current_segment)
                        if overlap:
                            overlap_segments += [current_segment, new_segment]
                            current_substracted_segments = current_segment.substract_segments([overlap])
                            new_substracted_segments = new_segment.substract_segments([overlap])
                            substract_segments += current_substracted_segments + new_substracted_segments
                new_segments = [ seg for seg in new_segments if seg not in overlap_segments ]
                new_segments += substract_segments
                # We must find all current perimeter segments which will be modified by the push
                # Find the whole segment where the pushed segment comes from, which will always be modified
                deformed_segment = next( seg for seg in self.perimeter.segments if pushed_segment in seg )
                # Substract the pushed segment from the whole segment
                # WARNING: This step must be done after the subtsraction of overlapped segments
                new_segments += deformed_segment.substract_segments([ pushed_segment ])
                # Find the rest of segments in the perimeter that will stay the same
                modified_segments = [ deformed_segment, *overlap_segments ]
                remaining_segments = [ seg for seg in self.perimeter.segments if seg not in modified_segments ]
                # Join all previous segments to make the new perimeter
                new_perimeter = Perimeter.non_canonical([ *remaining_segments, *new_segments ])
                # In case the perimeter is extended over another room,
                # We must substract the claimed rect from the other room and make it expand to compensate
                if room:
                    # Make a backup of the current perimeter in case the further expansions fail an we have to go back
                    backup_perimeter = self.perimeter
                    self.perimeter = new_perimeter
                    # Substract the claimed rect from the other room
                    exclusion_perimeters = [ Perimeter(new_rect.segments) ]
                    truncated_grid = Grid( room.perimeter.split_in_rectangles( exclusion_perimeters ) )
                    truncated_perimeters = truncated_grid.find_perimeters()
                    # In case the invaded room has been splitted in 2 parts as a result of this expansion we go back
                    if len(truncated_perimeters) > 1:
                        print('WARNING: The room has been splitted -> Go back')
                        self.perimeter = backup_perimeter
                        return False
                    # Modify the invaded room perimeter
                    room.perimeter = truncated_perimeters[0]
                    self.update_display()
                    # Expand the invaded room as much as the invaded area
                    # In case the expansions fails go back
                    if not room.expand():
                        print('WARNING: The room cannot expand -> Go back')
                        self.perimeter = backup_perimeter
                        return False
                else:
                    self.perimeter = new_perimeter
                    self.update_display()
                return True
                 

            # ----------------------------------------------------------------------------------------------------
            # First of all check if frontier points are connected to the space limits
            # In this case, we do not have to bother about margins
            # i.e. minimum length between the claimed space and space borders
            # Otherwise, we must check that borders are respected
            # ----------------------------------------------------------------------------------------------------

            # Sort points using the space contact 'a' point as reference
            points = [ frontier.a, frontier.b ]
            def by_distance(point):
                return space_contact.a.get_distance_to(point)
            sorted_points = sorted(points, key=by_distance)
            frontier_a = sorted_points[0]
            frontier_b = sorted_points[1]

            # Define margins at both sides of the current frontiers
            margin_a = Segment(space_contact.a, frontier_a) if space_contact.a != frontier_a else None
            margin_b = Segment(space_contact.b, frontier_b) if space_contact.b != frontier_b else None

            # If a margin exists and it is not as long as required we have a problem
            problem_a = margin_a and margin_a.length < margin_limit
            problem_b = margin_b and margin_b.length < margin_limit

            # If there is no problem we can just push the frontier
            if not problem_a and not problem_b:
                return push_segment(frontier)

            # ----------------------------------------------------------------------------------------------------
            # If a margin is not respected then we have 2 options (no option will always be possible):
            # - Cut the frontier to be expanded in order to respect the margin
            #   * If the cutted frontier is shorter than the minimum size we can not expand
            # - Claim also all the space between the claimed space and the space border
            #   * If the expanded space perpendicular space is not equal or longer to the minimum we cannot expand sideways
            #   * If the extra claimed space exceeds the required expand area we cannot claim it
            # ----------------------------------------------------------------------------------------------------

            print('PROBLEMA')
            return False

        # ----------------------------------------------------------------------------------------------------
        # Split the room perimeter in segments according to what each region is connected to
        # Regions connected to the parent free space are desired to expand
        # Regions connected to other rooms may be expanded if the colliding room is able to expand also
        # Regions connected to the parent limits will never be expanded
        # ----------------------------------------------------------------------------------------------------

        # The parent limits are not allowed for expansion
        parent_frontiers = self.get_frontiers(parent_room)
        # Other rooms inside the same parent may be displaced if there is no free space available
        brother_rooms = [ room for room in parent_room.children if room is not self ]
        brother_frontiers = []
        for room in brother_rooms:
            frontiers = self.get_frontiers(room)
            if frontiers:
                brother_frontiers += frontiers
        # The prefered limits to expand are those connected to free space inside the current parent
        # First of all get all already found frontier segments
        all_frontiers = [*parent_frontiers, *brother_frontiers]
        # Now susbstract all frontiers to each room perimeter segments to find free frontiers
        free_frontiers = []
        for segment in self.perimeter.segments:
            free_frontiers += segment.substract_segments(all_frontiers)
        # Set the room all free segment belong to as None
        for segment in free_frontiers:
            segment.room = None

        # Conserve only frontiers which reach the minimum size length
        # As an exception, frontiers connected to at least 1 inside corner may be expanded with no size restrictions
        inside_corners = [ corner for corner in self.perimeter.corners if corner.inside == True ]
        def is_suitable (frontier : Segment) -> bool:
            return frontier.length >= self.min_size or frontier.a in inside_corners or frontier.b in inside_corners
        suitable_free_frontiers = [ seg for seg in free_frontiers if is_suitable(seg) ]
        suitable_brother_frontiers = [ seg for seg in brother_frontiers if is_suitable(seg) ]

        #----------------------------------------------------------------------------------------------
        # Set now the strategy to sort the frontiers (i.e. to choose which frontiers will be tried first)
        # This is a critical step, since many frontiers may be expanded but only a few may be useful to expand
        # There are many possible situtations to take in count
        # A wrong strategy could lead to an endless loop or dead end
        # A good strategy could make the algorithm faster
        #----------------------------------------------------------------------------------------------

        # Random (no strategy)
        random.shuffle(suitable_free_frontiers)
        random.shuffle(suitable_brother_frontiers)

        # Try to expand one by one all frontiers until one of them is expanded sucessfully
        # Once an expansion has been done we must recalculate frontiers if we want to keep expanding
        # Priorize free frontiers over bother frontiers
        available_frontiers = [ *suitable_free_frontiers, *suitable_brother_frontiers ]
        for frontier in available_frontiers:
            if expand_frontier(frontier):
                return True
        return False


    # Get all overlapped segments between the current room and others
    def get_frontiers (self, other) -> list:
        self_segments = self.perimeter.segments
        other_segments = other.perimeter.segments
        overlap_segments = []
        for self_segment in self_segments:
            for other_segment in other_segments:
                overlap_segment = self_segment.get_overlap_segment(other_segment)
                if overlap_segment:
                    # Save the room this segment belongs to inside the segment object
                    overlap_segment.room = other
                    overlap_segments.append(overlap_segment)
        return overlap_segments

    # Go uppwards in the hyerarchy until you reach the room which has no parent
    def get_root_room (self):
        root = self
        while(root.parent):
            root = root.parent
        return root

    # Get all segments from this room and all children room perimeters
    def get_rooms_recuersive (self, only_children : bool = False):
        rooms = []
        if not only_children and self.perimeter:
            rooms.append(self)
        for room in self.children:
            rooms.append(room)
        return rooms

    # Add a new frame in the display with the current segments of this room and its children
    # Also an 'extra' argument may be passed with extra segments to be represented
    def update_display (self, extra : list = []):
        # Find the root room
        root = self.get_root_room()
        if root.display:
            elements_to_display = [ *root.get_rooms_recuersive(), *extra ]
            add_frame(elements_to_display)