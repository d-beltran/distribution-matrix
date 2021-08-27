from typing import List, Union, Optional

from scheme_display import add_frame

from vectorial_base import *

import random
from math import sqrt, inf

# Set the seed and print it
seed = None
#seed = 168166 # burbuja imposible perimetro 4
#seed = 414509 # infinite loop perimetro 4
#seed = 402486 # infinite loop perimetro 1
if not seed:
    seed = round(random.random() * 999999)
print('Seed ' + str(seed))
random.seed(seed)
        
# A room is a smart boundary that may contain other boundary with conservative areas and size restrictions
# A start boundary may be passed. If no boundary is passed it is assigned automatically according the room rules
# The 'forced_area' argument stablishes the expected final room area. If no area is passed then the original boundary area will be used
# Forced area may be a number (absolute value) or a string (percent) with the format 'XX%'
# The 'min_size' and 'max_size' are the limits in both x and y axes
# The name and color parameters are only representation parameters and they have no effect in the logic
class Room:
    def __init__ (self,
        boundary : Optional[Boundary] = None,
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
        self._boundary = None
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
        # Set the boundary
        # If the boundary has been forced then update the display with the initial segments
        if boundary:
            self.boundary = boundary
        else:
            self.boundary = None
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

    def __str__(self):
        return '<Room "' + str(self.name) + '">'

    def __repr__(self):
        return '<Room "' + str(self.name) + '">'

    # Get the boundary
    # Just return the internal boundary value
    def get_boundary (self):
        return self._boundary

    # Set the boundary
    # Reset the own rects and the parent rects also
    def set_boundary (self, boundary):
        self._boundary = boundary
        self.reset_free_grid()
        if self.parent:
            self.parent.reset_free_grid()
        self.update_display()

    # The room boundary
    boundary = property(get_boundary, set_boundary, None, "The room boundary")

    # Area inside the room boundary (read only)
    def get_area(self):
        boundary = self.boundary
        if not boundary:
            return 0
        return boundary.area
    area = property(get_area, None, None, "Area inside the room boundary")

    # Get the available space inside the room boundary as a rectangles grid
    # i.e. space not filled by children rooms
    def get_free_grid (self):
        # If rects are previously calculated then return them
        if self._free_grid:
            return self._free_grid
        # Return none if there is not boundary yet
        if not self.boundary:
            return None
        # If there are no children then return the current boundary grid
        # If all children have no boundary then return the current boundary grid
        if len(self.children) == 0 or not any([ child.boundary for child in self.children ]):
            free_grid = self.boundary.grid
        # Otherwise, find out the free space grid inside the room boundary
        # Add the children room boundary exterior polygons to self room boundary interior polygons
        else:
            children_polygons = [ child.boundary.exterior_polygon for child in self.children if child.boundary ]
            free_boundary = Boundary(self.boundary.exterior_polygon, self.boundary.interior_polygons + children_polygons)
            free_grid = free_boundary.grid
        # Apply the current room colors to all rectangles
        for rect in free_grid.rects:
            rect.segments_color = self.segments_color
            rect.fill_color = self.fill_color
        self._free_grid = free_grid
        # ---------------------------------------------------------------------------------------------------
        # DANI: Muestra los rects la primera vez que se calculan
        colored_rects = [ rect.get_colored_rect(segments_color='red') for rect in free_grid.rects ]
        #self.update_display(colored_rects)
        # ---------------------------------------------------------------------------------------------------
        return free_grid
    # Free space grid (read only)
    free_grid = property(get_free_grid, None, None, "The room free space grid")

    # Reset all minimum and maximum free rects
    # This must be done each time the boundary is modified since they are not valid anymore
    def reset_free_grid (self):
        self._free_grid = None

    # Get rectangles in the free grid
    def get_free_rects (self):
        return self.free_grid.rects
    free_rects = property(get_free_rects, None, None, "The room free area splitted in rectangles")

    # Get the maximum rectangles from the free grid
    def get_free_max_rects (self):
       return self.free_grid.max_rects
    free_max_rects = property(get_free_max_rects, None, None, "The free maximum rectangles")

    # Free space area (read only)
    def get_free_area (self):
        return self.free_grid.area
    free_area = property(get_free_area, None, None, "Free space area (read only)")

    # Get the maximum size
    # If there is a forced maximum size then return it
    # Otherwise calculate it by dividing the forced area between the minimum size
    # i.e. return the maximum possible size in case the boundary was the thinest rectangle
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

    # Check if a room fits in this room according to its minimum size
    # Check free space by default and all space if the argument 'force' is passed
    def does_room_fit (self, room : 'Room', force : bool = False) -> bool:
        size = room.min_size
        # If there is no minimum size it will always fit
        if not size:
            return True
        # Generate fitting rects. If there is at least one then the room fits
        grid = self.boundary.grid if force else self.free_grid
        fitting_rects = grid.get_fitting_space(size, size)
        if next(fitting_rects, None):
            return True
        return False

    # Search all maximum rects where the specified room minimum rectangle would fit
    # Search in free space by default and all space if the argument 'force' is passed
    def get_room_fitting_rects (self, room : 'Room', force : bool = False) -> List[Rect]:
        size = room.min_size
        # Get all fitting rects
        grid = self.boundary.grid if force else self.free_grid
        fitting_rects = grid.get_fitting_space(size, size)
        return list(fitting_rects)

    # Add children rooms
    def add_children (self, rooms : list):
        if len(rooms) == 0:
            return
        # Check areas of all children rooms to do not sum up more than the parent area
        # In addition check if any of the children room has boundary and, if so, check the boundary is inside the parent
        children_area = 0
        for room in rooms:
            # If the children has a percent forced area this is the time to calculate the absolute forced area
            if room._forced_area_portion:
                room.forced_area = self.forced_area * room._forced_area_portion
            children_area += room.forced_area
        if children_area > self.area:
            raise InputError('Children together require more area than the parent has')
        # Check all children are inside the parent boundary, if they have a predefined boundary
        for room in rooms:
            if room.boundary and room.boundary not in self.boundary:
                raise InputError('The child room "' + room.name + '" is out of the parent boundary')
        # Check all children minim sizes fit in the parent boundary, if they have a predefined minimum size
        for room in rooms:
            if not self.does_room_fit(room, force=True):
                raise InputError('The child room "' + room.name + '" minimum size does not fit in the parent boundary')

        # Sort children rooms by minimum size, with the biggest sizes first
        def sort_by_size (room):
            return room.min_size
        sorted_rooms = sorted( rooms, key=sort_by_size, reverse=True )

        # Set up each room by giving them a position and correct size to match the forced area
        for room in sorted_rooms:
            # Update the room hierarchy
            room.parent = self
            self.children.append(room)
            # If the children has no boundary it must be built
            if not room.boundary:
                if not self.set_child_room_boundary(room):
                    raise RuntimeError('Child ' + room.name + ' has failed to be set')

    # Set up a child room boundary
    def set_child_room_boundary (self, room) -> bool:
        # Find a suitable maximum free rectangle to deploy a starting base boundary
        # The minimum base boundary is a square with both sides as long as the room minimum size
        suitable_rects = self.get_room_fitting_rects(room)
        # If there are no suitable rects it means this child fits nowhere in the free space
        # This may happen with the last childs, when previous childs have take almost all space
        # In this case we take as availbale space all the room space and then we invade overlapped children
        forced = False
        if len(suitable_rects) == 0:
            print('WARNING: The room ' + room.name + ' fits nowhere in the free space')
            suitable_rects = self.get_room_fitting_rects(room, force=True)
            forced = True
            if len(suitable_rects) == 0:
                # DANI: Esto no debería pasar nunca. Debería preveerse de antes
                raise RuntimeError('The room ' + room.name + ' fits nowhere')
        # Shuffle the suitable rects
        random.shuffle(suitable_rects)
        # Sort the suitable rects by minimum size, with the biggest sizes first
        def sort_by_size(rect):
            return min(rect.get_size())
        sorted_suitable_rects = sorted( suitable_rects, key=sort_by_size, reverse=True )
            
        # Try to set up the new room in all possible sites until we find one
        # Each 'site' means each corner in each suitable rectangle
        # Check each site to allow other rooms to grow
        sites = [ (corner, rect) for rect in sorted_suitable_rects for corner in rect.get_corners() ]
        for corner, rect in sites:
            # In case we forced the base boundary we must check which children were invaded
            # In addition, the base boundary must be as small as possible
            if forced:
                initial_boundary = room.set_minimum_initial_boundary(corner, rect)
            # Otherwise we set freely the maximum possible boundary
            else:
                initial_boundary = room.set_maximum_initial_boundary(corner, rect)
            # There must be always an initial boundary at this point
            if not initial_boundary:
                raise RuntimeError('Something went wrong with initial boundary')
            # Set the child first boundary, which automatically will reset self room free grid
            # In case it was forced, we must check that the overlapped children rooms are fine with the invasion
            if forced:
                if not self.invade_children(room, initial_boundary):
                    room.boundary = None
                    continue
            # Otherwise just claim the boundary
            else:
                room.boundary = initial_boundary
            # If the child room has no initial boundary at this point there must be something wrong
            if not room.boundary:
                raise RuntimeError('Failed to set an initial boundary for room "' + room.name + '"')
            # Proceed with the expansion of this child room until it reaches its forced area
            if not room.expand():
                room.boundary = None
                continue
            break
        # If at this point we still have no boundary it means there is no available place to set the perimeter
        if room.boundary == None:
            return False
        return True

    # Set the initial room boundary as the maximum possible rectangle
    # This is a shortcut to skip the difficult expansion protocol as much as possible
    # It is useful to set a whole room at the begining, when there is plenty of free space
    def set_maximum_initial_boundary (self, corner : Point, space : Rect) -> Boundary:
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
                return Boundary(Polygon.from_rect(maximum_rect))
            return Boundary(Polygon.from_rect(space))
        # Else, fit the room in the space
        # First of all create the 3 rule rectangle from the space
        # i.e. calculate the relation of areas and apply it to the square root to both x and y sizes
        area_relation = self.forced_area / space.area # At this point the relation is always < 1
        new_x_size = x_space * sqrt(area_relation)
        new_y_size = y_space * sqrt(area_relation)
        # If any of the new sizes is shorter than the maximum size then the rectangle is valid
        if new_x_size < self.max_size or new_y_size < self.max_size:
            maximum_rect = Rect.from_corner(corner, new_x_size, new_y_size)
            return Boundary(Polygon.from_rect(maximum_rect))
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
        return Boundary(Polygon.from_rect(maximum_rect))

    # Set the initial room boundary as the minimum possible rectangle
    # This is used when the initial boundary must be forced over other children rooms
    def set_minimum_initial_boundary (self, corner : Point, space : Rect) -> Boundary:
        minimum_rect = Rect.from_corner(corner, self.min_size, self.min_size)
        return Boundary(Polygon.from_rect(minimum_rect))

    # Calculate how much area we need to expand
    def get_required_area (self):
        return resolute(self.forced_area - self.area, -2)

    # Expand this room until it reaches the forced area
    # In case it is not able to expand at some point recover the original situation
    def expand (self) -> bool:
        # Make a boundary bakcup of this room and all its brothers
        affected_rooms = self.parent.children
        boundaries_backup = [ room.boundary for room in affected_rooms ]
        while self.get_required_area() > 0:
            if not self.expand_step():
                for i, room in enumerate(affected_rooms):
                    room.boundary = boundaries_backup[i]
                return False
        return True

    # Find the most suitable space for the current room to claim
    # If current room is expanded over another room then the other room must also expand to compensate
    # At the end, the extra space will be substracted from the parent free space
    def expand_step (self) -> bool:

        print('Expanding ' + self.name)
        # Calculate how much area we need to expand
        required_area = self.get_required_area()
        if required_area <= 0:
            raise ValueError('Something went wrong with expansion')

        # Check the required area is big enought to be meaningfull according to the resolution
        # i.e. check if expanding a segment with the minimum length would make it move more than the minimum resolution
        if required_area / self.min_size < minimum_resolution:
            print('WARNING: The required area is not enought to expand according to resolution: ' + str(required_area))
            return False

        # Get the parent room
        parent_room = self.parent
        # Get the exterior polygon of the room boundary
        exterior_polygon = self.boundary.exterior_polygon
        # Get the inside corners
        inside_corners = [ corner for corner in exterior_polygon.corners if corner.inside == True ]

        # ----------------------------------------------------------------------------------------------------
        # Try to expand in the easiest way: push an individual frontier segment
        # Only segment longer than the minimum size may be expanded this way
        # Start pushing free frontiers
        # - No matter if the space is divided as a result of the expansion
        # If there are no suitable free frontiers then proceed with brother frontiers
        # - Space must never be divided as a result of the expansion
        # ----------------------------------------------------------------------------------------------------

        def expand_frontier (frontier : Segment) -> bool:
            print("Push -> " + str(frontier))
            # Set some parameters according to the rooms this frontier belongs to
            rooms = frontier.rooms
            # Get the grid of the rooms we are about to invade as a reference to calculate how much we can expand
            # When the invaded room is self room it means we are expanding over free space
            first_room = rooms[0]
            grid = first_room.free_grid if first_room == parent_room else first_room.boundary.grid
            # In case this frontier has more than 1 rooms it means it is a combined frontier
            # Create a 'ficticious' room with all implicated rooms as a reference
            for next_room in rooms[1:]:
                next_grid = next_room.free_grid if next_room == parent_room else next_room.boundary.grid
                grid = grid.get_merge_grid(next_grid)
            margin_limit = max([ room.min_size for room in rooms ])
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
                raise ValueError('ERROR: diagonal segments are not supported for room expansion')

            if len(rects) == 0:
                raise RuntimeError('No pot ser: ' + str(grid))

            # One and only one of the rows/columns will always include the segment
            space = next((rect for rect in rects if frontier in rect), None)
            if space == None:
                # If this happens it may mean there is a problem with the grid
                raise RuntimeError('Frontier ' + str(frontier) + ' has no space in ' + str(rects))
            space_contact = next(segment for segment in space.segments if frontier in segment)
            space_forward_limit = space.get_size()[forward]
            margined_space_forward_limit = space_forward_limit - margin_limit

            # Find the direction of the expansion as a vector
            forward_direction = space.get_side_direction_to_center(space_contact)

            # Get the forward expansion limit according to maximum rectangles
            maximum_forward_limit = max([ rect.get_size()[forward] for rect in contact_max_rects ])
            margined_maximum_forward_limit = maximum_forward_limit - margin_limit

            # Now create a function to make a rectangle by pushing a segment (which may change)
            # The forward length of this rectangle will depend on the available area and limits
            # This segment may not be the original frontier, but a variation of it
            def push_segment (pushed_segment : Segment, greedy : bool = True) -> bool:
                push_length = required_area / pushed_segment.length
                # In case this is an insider segment which is not wide enought to be pushed alone,
                # Find out how much we can push this segment
                # i.e. find the connected frontier/s and get the maximum length of these segments
                if pushed_segment.length < self.min_size:
                    push_limit = None
                    other_segments = [ segment for segment in exterior_polygon.segments if segment != pushed_segment ]
                    for point in pushed_segment.points:
                        if point in inside_corners:
                            insider = next(segment for segment in other_segments if point in segment.points)
                            if not push_limit or push_limit < insider.length:
                                push_limit = insider.length
                    # If the push length exceeds the insider limit then stay in the limit
                    if push_length > push_limit:
                        push_length = push_limit
                # First, try to expand using the maximum available space
                # This is risky since it may split the invaded room in two parts or make regions which do not respect the minimum size
                if greedy:
                    # If the length exceeds the maximum limit then stay in the maximum limit
                    if push_length > maximum_forward_limit:
                        push_length = maximum_forward_limit
                    # If the length is between the maximum limit and the margined maximum limit then stay at the margin
                    elif push_length > margined_maximum_forward_limit:
                        push_length = margined_maximum_forward_limit
                # If the greedy try fails use the moderate try
                # Expand only using the closest row/column rectangle
                # This expansion is smaller but safe
                else:
                    # # If at this point the length exceeds the space limit then stay in the space limit
                    if push_length > space_forward_limit:
                        push_length = space_forward_limit
                    # If the length is between the space limit and the margined space limit then stay at the margin
                    elif push_length > margined_space_forward_limit:
                        push_length = margined_space_forward_limit
                    # In case the length is bigger than the margined maximum we stay at the margin
                    if push_length > margined_maximum_forward_limit:
                        push_length = margined_maximum_forward_limit
                # If the push length at this point is 0 or close to it then we can not push
                # WARNING: This usually happens because of limits, not because the area was not big enought
                # For this reason, trying to reduce the segment and push again will have no effect almost always
                if push_length < minimum_resolution:
                    print('WARNING: The push length is too small')
                    return False
                # Create the new rect with the definitive length
                new_point = pushed_segment.a + forward_direction.normalized() * push_length
                new_side = Segment(pushed_segment.a, new_point)
                new_rect = Rect.from_segments([pushed_segment, new_side])
                # Then add the new current rectangle by modifying the current boundary
                # Add the new rectangle segment regions which do not overlap with current boundary segments
                # Substract the new rectangle segment regions which overlap with current boundary segments
                # e.g. the pushed segment will always be an overlapped region
                new_segments = new_rect.segments
                current_segments = self.boundary.exterior_polygon.segments
                added_segments = [ seg for segment in new_segments for seg in segment.substract_segments(current_segments) ]
                remaining_segments = [ seg for segment in current_segments for seg in segment.substract_segments(new_segments) ]
                # Join all previous segments to make the new polygon
                new_exterior_polygon = Polygon.non_canonical([ *added_segments, *remaining_segments ])
                new_boundary = Boundary(new_exterior_polygon, self.boundary.interior_polygons)
                # In case the boundary is extended over another room,
                # We must substract the claimed rect from the other room and make it expand to compensate
                # Make a backup of the current boundary in case the further expansions fail an we have to go back
                backup_boundary = self.boundary
                self.boundary = new_boundary
                invaded_region = Boundary(Polygon.from_rect(new_rect))
                # Substract the claimed rect from other rooms
                # Make a backup of all other room current boundaries
                backup_room_boundaries = [ room.boundary for room in rooms ]
                for room in rooms:
                    if room == parent_room:
                        continue
                    current_room_invaded_regions = room.boundary.get_overlap_boundaries(invaded_region)
                    if not room.invade(current_room_invaded_regions):
                        self.boundary = backup_boundary
                        for i, modified_room in enumerate(rooms):
                            modified_room.boundary = backup_room_boundaries[i]
                        # If the invasion fails then try it again with the non greedy push
                        # If we are already in the non greedy push then return false
                        if greedy:
                            return push_segment(pushed_segment, False)
                        return False
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
            point_a = sorted_points[0]
            point_b = sorted_points[1]

            # Define margins at both sides of the current frontiers
            margin_a = Segment(space_contact.a, point_a) if space_contact.a != point_a else None
            margin_b = Segment(space_contact.b, point_b) if space_contact.b != point_b else None

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

            # Try to cut the frontier
            reduced_frontier = frontier
            for i, problem in enumerate([ problem_a, problem_b ]):
                if not problem:
                    continue
                margin = [ margin_a, margin_b ][i]
                reduction = margin_limit - margin.length
                problem_point = sorted_points[i]
                other_point = next(point for point in reduced_frontier.points if point != problem_point)
                direction = (problem_point + other_point).normalized()
                new_point = problem_point + direction * reduction
                # In case the new pont has passed the other point we stop here
                if new_point not in reduced_frontier:
                    print('WARNING: The frontier has been fully consumed')
                    return False
                if new_point.get_distance_to(other_point) < self.min_size:
                    print('WARNING: The reduced frontier is not wide enought')
                    return False
                reduced_frontier = Segment(new_point, other_point)

            return push_segment(reduced_frontier)

        # ----------------------------------------------------------------------------------------------------
        # Split the room boundary exterior polygon in segments according to what each region is connected to
        # Regions connected to the parent free space are desired to expand
        # Regions connected to other rooms may be expanded if the colliding room is able to expand also
        # Regions connected to the parent limits will never be expanded
        # ----------------------------------------------------------------------------------------------------

        free_frontiers, brother_frontiers, parent_frontiers = self.get_frontiers()
        expandable_frontiers = free_frontiers + brother_frontiers

        #----------------------------------------------------------------------------------------------
        # Set now the strategy to sort the frontiers (i.e. to choose which frontiers will be tried first)
        # This is a critical step, since many frontiers may be expanded but only a few may be useful to expand
        # There are many possible situtations to take in count
        # A wrong strategy could lead to an endless loop or dead end
        # A good strategy could make the algorithm faster
        #----------------------------------------------------------------------------------------------

        # Make a function to check which frontiers are suitable for expansion alone
        # It means frontiers whose length reaches the self room minimum size
        # As an exception, short frontiers connected to at least 1 inside corner may be expanded
        # However these frontiers will have a push length limit
        def is_suitable (frontier : Segment) -> bool:
            return frontier.length >= self.min_size or frontier.a in inside_corners or frontier.b in inside_corners

        # Short frontiers could be also expanded together with other connected and aligned frontiers
        # These compound frontiers must be taken in count also although they are harder to expand
        # Each frontier may form 0, 1 or 2 combined frontiers
        def get_combined_frontiers (frontier : Segment) -> List[Segment]:
            combined_frontiers = []
            for point in [ frontier.a, frontier.b ]:
                other_frontiers = [ seg for seg in expandable_frontiers if seg != frontier ]
                combined_frontier = frontier
                implicated_rooms = frontier.rooms
                next_point = point
                # Join connected frontiers until we have enough length or the next frontier is not aligned
                while combined_frontier.length < self.min_size:
                    # Find the next connected frontier
                    # There may be no next connected frontier if next would be a parent frontier
                    connected_frontier = next((seg for seg in other_frontiers if seg.has_point(next_point)), None)
                    # If it is not aligned we can not combine it so we stop here
                    if not connected_frontier or not connected_frontier.same_line_as(combined_frontier):
                        break
                    # Otherwise remove the connected frontier from the others list and combine both frontiers
                    other_frontiers = [ seg for seg in other_frontiers if seg != connected_frontier ]
                    combined_frontier = combined_frontier.combine_segment(connected_frontier)
                    implicated_rooms += connected_frontier.rooms
                    next_point = next(point for point in connected_frontier.points if point != next_point)
                # If the combined frontier reaches the minimum length add it to the list to be returned
                if combined_frontier.length >= self.min_size:
                    # Add all combined frontier rooms to the combined frontier
                    combined_frontier.rooms = unique(implicated_rooms)
                    combined_frontiers.append(combined_frontier)
            return combined_frontiers

        # Start trying to expand free frontiers first
        # Check if frontiers are available for expansion alone
        # In case they are not, skip them and let them for the end
        for frontiers_groups in [ free_frontiers, brother_frontiers ]:
            # Sort frontiers according to their type
            # First of all and always we shuffle them at random
            random.shuffle(frontiers_groups)
            # Then, in case they are brother frontiers, we sort them with some extra logic
            if frontiers_groups == brother_frontiers:
                # Find for each brother room the number of colliding rooms we must jump to find free space
                # Then use this value to set the "score" of each brother room frontiers
                colliding_rooms = unique([ frontier.rooms[0] for frontier in frontiers_groups ])
                for colliding_room in colliding_rooms:
                    previous_rooms = [ self ]
                    current_rooms = [ colliding_room ]
                    counter = 1
                    searching_free = True
                    # Get frontiers from all current rooms
                    # If any of them has free frontiers we are done
                    # Otherwise, get the rooms from all brother rooms and repeat the whole process
                    while True:
                        current_frontiers = []
                        for current_room in current_rooms:
                            free, brother, parent = current_room.get_frontiers()
                            if len(free) > 0:
                                searching_free = False
                                break
                            current_frontiers += brother
                        if not searching_free:
                            print('ROOM ' + colliding_room.name + ' -> SCORE ' + str(counter))
                            break
                        previous_rooms = previous_rooms + current_rooms
                        current_rooms = unique([ frontier.rooms[0] for frontier in current_frontiers if frontier.rooms[0] not in previous_rooms ])
                        if len(current_rooms) == 0:
                            counter = None
                            print('ROOM ' + colliding_room.name + ' -> DEAD END')
                            break
                        counter += 1
                    # Set the scores for all frontiers
                    colliding_room_frontiers = [ frontier for frontier in frontiers_groups if frontier.rooms[0] == colliding_room ]
                    if counter:
                        for frontier in colliding_room_frontiers:
                            frontier.score = counter
                    else:
                        for frontier in colliding_room_frontiers:
                            frontier.score = inf
                # Now sort frontiers using previous scores
                def by_score (frontier):
                    return frontier.score
                frontiers_groups = sorted(frontiers_groups, key=by_score)
            hard_frontiers = []
            for frontier in frontiers_groups:
                # If it cannot be expanded alone we skip it by now
                if not is_suitable(frontier):
                    hard_frontiers.append(frontier)
                    continue
                # Otherwise try to expand it
                if expand_frontier(frontier):
                    return True
            # If none of the suitable frontiers was expanded successfully then try with the combined ones
            for frontier in hard_frontiers:
                combined_frontiers = get_combined_frontiers(frontier)
                for combined_frontier in combined_frontiers:
                    if expand_frontier(combined_frontier):
                        return True
        print('WARNING: There are no more frontiers available')
        return False

    # Invade this room by substracting part of its space
    # Then this room must expand to recover the lost area
    def invade (self, regions : List[Boundary]) -> bool:
        # Calculate the boundary of this room after substracting the invaded regions
        truncated_grid = self.boundary.grid
        for region in regions:
            truncated_grid = truncated_grid.get_substract_grid(region.grid)
        # Check the truncated grid to still respecting the minimum size
        if not truncated_grid.check_minimum(self.min_size):
            print('WARNING: The room is not respecting the minimum size')
            return False
        truncated_boundaries = truncated_grid.find_boundaries()
        # In case the room has been splitted in 2 parts as a result of the invasion we go back
        if len(truncated_boundaries) > 1:
            print('WARNING: The room has been splitted -> Go back')
            return False
        # DANI: Es posible que se coma toda la habitación??
        if len(truncated_boundaries) == 0:
            print('WARNING: The invaded room has been fully consumed -> Go back')
            return False
        # Modify the room boundary but save a backup in case we have to go back further
        backup_boundary = self.boundary
        self.boundary = truncated_boundaries[0]
        # Expand the invaded room as much as the invaded area
        # In case the expansions fails go back
        if not self.expand():
            print('WARNING: The room cannot expand -> Go back')
            self.boundary = backup_boundary
            return False
        return True

    # Invade children in this room by substracting part of their space
    # Then all children must expand to recover the lost area
    # The invasor is the children room which will claim this new space
    def invade_children (self, invasor : 'Room', region : Boundary) -> bool:
        # Make a backup of all children, which includes the invasor
        # All children boundaries will be recovered if only 1 child fails to get invaded
        children_backup = [ child.boundary for child in self.children ]
        # Claim the invaded region for the invasor room
        invasor.merge_boundary(region)
        # Now find which children are overlaped by the invade region and then invade them
        for child in self.children:
            # Skip the invasor room
            if child == invasor:
                continue
            # Get the overlap regions between the invaded region and this child
            overlap_boundaries = child.boundary.get_overlap_boundaries(region)
            if len(overlap_boundaries) == 0:
                continue
            # In case something went wrong with the invasion recover all children backups and stop
            if not child.invade(overlap_boundaries):
                for c, child in enumerate(self.children):
                    child.boundary = children_backup[c]
                return False
        return True

    # Fuse a boundary to self room boundary
    # Check that boundary can be joined to current boundary as a single room
    # i.e. both boundaries must be colliding and the colliding region must be only one and as wide as the minimum size or more
    def merge_boundary (self, boundary : Boundary):
        if self.boundary == None:
            self.boundary = boundary
        else:
            self.boundary = self.boundary.merge_boundary(boundary, self.min_size)

    # Get all overlapped segments between the current room and others
    def get_frontiers_with (self, other : 'Room') -> list:
        self_segments = self.boundary.exterior_polygon.segments
        other_segments = other.boundary.exterior_polygon.segments
        overlap_segments = []
        for self_segment in self_segments:
            for other_segment in other_segments:
                overlap_segment = self_segment.get_overlap_segment(other_segment)
                if overlap_segment:
                    # Save the room this segment belongs to inside the segment object
                    overlap_segment.rooms = [other]
                    overlap_segments.append(overlap_segment)
        return overlap_segments

    # Get all self frontiers segments
    # They come inside a tuple separated by free frontiers, brother frontiers and parent frontiers respectively
    def get_frontiers (self) -> tuple:
        parent_room = self.parent
        exterior_polygon = self.boundary.exterior_polygon
        # The parent limits are not allowed for expansion
        parent_frontiers = self.get_frontiers_with(parent_room)
        # Other rooms inside the same parent may be displaced if there is no free space available
        brother_rooms = [ room for room in parent_room.children if room is not self ]
        brother_frontiers = []
        for room in brother_rooms:
            # This function assign the frontier.rooms already
            frontiers = self.get_frontiers_with(room)
            if frontiers:
                brother_frontiers += frontiers
        # The prefered limits to expand are those connected to free space inside the current parent
        # First of all get all already found frontier segments
        conflict_frontiers = [*parent_frontiers, *brother_frontiers]
        # Now susbstract all previous frontiers from the rest of the external polygon to find free frontiers
        free_frontiers = []
        for segment in exterior_polygon.segments:
            free_frontiers += segment.substract_segments(conflict_frontiers)
        # Set the room all free segment belong to as None
        for frontier in free_frontiers:
            frontier.rooms = [parent_room]
        return free_frontiers, brother_frontiers, parent_frontiers

    # Go uppwards in the hyerarchy until you reach the room which has no parent
    def get_root_room (self):
        root = self
        while(root.parent):
            root = root.parent
        return root

    # Get this room and all children rooms recursively
    # Get only rooms with boundary
    def get_rooms_recuersive (self, only_children : bool = False):
        rooms = []
        if not only_children and self.boundary:
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

# Exception for when user input is wrong
class InputError(Exception):
    pass