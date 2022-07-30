from typing import List, Tuple, Dict, Union, Optional

from scheme_display import add_frame

from vectorial_base import *

import random
from math import sqrt, inf

# Set the seed and print it
seed = None
#seed = 975359
if not seed:
    seed = round(random.random() * 999999)
print('Seed ' + str(seed))
random.seed(seed)

# Set overall solving options

# Set if the solving process must be displayed
display_solving_process = False
# Sets how much greater is the effective minimum size in the steps previous to corridor and wall tickness build
# The most preventive margin we have, the most respect for the minimum size but the less flexibility when solving the puzzle
# This means a strict protocol may have not solution
# 0 - Use the minimum size as it is
#     This is prone to not respect minimum size at several regions
#     Moderate risk of large room space loose during the corridor build
# 1 - Use the minimum size + half the corridor size
#     This is prone to not respect minimum size at a few regions
#     Lesser risk of large room space loose during the corridor build
# 2 - Use the minimum size + the corridor size
#     This should respect minimum size in most cases
#     Negligible risk of large room space loose during the corridor build
# 3 - Use the minimum size + twice the corridor size
#     This should respect minimum size even in the worst possible scenario
preventive_min_size_protocol = 0

# A door is a segment in a polygon
# When boundaries are transformed to walls with tickness, doors become holes in the wall
class Door:
    def __init__ (self,
        # A point may be passed. If no point is passed then it is assigned automatically
        # WARNING: This is only supported if the door room boundary is already set and the point matches on it
        point : Optional[Point] = None,
        # Set how wide the door must be
        width : Optional[number] = None,
        # Set the minimum width of margins on each side of the door
        margin : Optional[number] = None,
        # Set the direction the door is open thorugh
        direction : Optional[Vector] = None,
        # Set the pivot point for the door to open
        pivot : Optional[Point] = None
    ):  
        self.width = width
        self.margin = margin
        self._point = point
        self._margined_width = None
        self._segment = None
        self._margined_segment = None
        self._direction = direction
        self._pivot = pivot
        # The room this door belongs to
        self.room = None

    def __repr__ (self):
        point = str(self.point) if self.point else 'No point'
        width = str(self.width) if self.width else 'No width'
        margin = str(self.margin) if self.margin else 'No margin'
        return '<Door ' + point + ' ' + width + '(' + margin + ')' + '>'

    # Get the margined width
    def get_margined_width (self):
        if self._margined_width != None:
            return self._margined_width
        if self.width == None:
            raise ValueError('Door is missing width')
        if self.margin == None:
            raise ValueError('Door is missing margin')
        self._margined_width = self.width + self.margin * 2
        return self._margined_width

    # The door margined width
    margined_width = property(get_margined_width, None, None, "The door margined width")

    # Get the door point
    def get_point (self) -> Optional[Point]:
        return self._point

    # If the door point is set then reset its segment and margined segment
    def set_point (self, point : Optional[Point]):
        self._point = point
        self._segment = None
        self._margined_segment = None
        self._direction = None
        self._pivot = None
        if point:
            self.segment
            self.margined_segment
            self.direction
            self.pivot

    # The door segment
    point = property(get_point, set_point, None, "The door point")

    # Get the door segment
    def get_segment (self) -> Optional[Segment]:
        # Return internal value if it exists
        if self._segment:
            return self._segment
        self._segment = self.generate_segment(self.width)
        return self._segment

    # The door segment
    segment = property(get_segment, None, None, "The door segment")

    # Get the door margined segment
    def get_margined_segment (self) -> Optional[Segment]:
        # Return internal value if it exists
        if self._margined_segment:
            return self._margined_segment
        self._margined_segment = self.generate_segment(self.margined_width)
        return self._margined_segment

    # The door margined segment
    margined_segment = property(get_margined_segment, None, None, "The door margined segment")

    # Given a segment width, generate a new
    # The new segment will be centered in the door point
    # The new segment will be overlaped with the polygon segment where the door point is
    def generate_segment (self, width : number) -> Segment:
        # If the door point or polygon are not assigned we can not generate the segment
        room_polygon = self.get_room_polygon()
        if not room_polygon or not self.point:
            return None
        # Otheriwse, generate the margined segment
        polygon_segment = next(( segment for segment in room_polygon.segments if self.point in segment ), None)
        if not polygon_segment:
            raise ValueError('The door point ' + str(self.point) + ' is not over its polygon (' + self.room.name + ')')
        direction = polygon_segment.direction
        half_width = width / 2
        a = self.point - direction * half_width
        b = self.point + direction * half_width
        if a not in polygon_segment or b not in polygon_segment:
            raise ValueError('The door segment (' + str(Segment(a,b)) + ') does not fit in its room polygon (' + self.room.name + ')')
        return Segment(a,b)

    # Get the door direction
    # i.e. the direction the door is open thorugh
    # By default the direction points to the door polygon inside side
    def get_direction (self) -> Optional[Vector]:
        # Return internal value if it exists
        if self._direction:
            return self._direction
        if not self.segment:
            return None
        room_polygon = self.get_room_polygon()
        if not room_polygon:
            return None
        direction = room_polygon.get_border_inside(self.segment)
        self._direction = direction
        return direction
    # The direction crosses the door segment perpendicularly
    # It is a normalized vector
    direction = property(get_direction, None, None, "The door direction")

    # Get the door pivot
    # i.e. the point where the door would rotate around to get open / closed
    # By default the pivot is the door segment point which is closer to an outside corner
    def get_pivot (self) -> Optional[Point]:
        # Return internal value if it exists
        if self._pivot:
            return self._pivot
        if not self.segment:
            return None
        room_polygon = self.get_room_polygon()
        if not room_polygon:
            return None
        door_points = self.segment.points
        polygon_segment = next( segment for segment in room_polygon.segments if self.point in segment )
        polygon_segment_points = polygon_segment.points
        outside_corners = [ corner for corner in room_polygon.corners if corner in polygon_segment_points and not corner.inside ]
        # If both segment corners are inside corners then just set the first segment point as the pivot
        if len(outside_corners) == 0:
            self._pivot = door_points[0]
            return self._pivot
        # Get minimum distance to a closer corner for both door segment points and then select the shortest distance point as the pivot
        minimum_distances = []
        for point in door_points:
            distances = [ point.get_distance_to(segment_point) for segment_point in outside_corners ]
            minimum_distance = min(distances)
            minimum_distances.append(minimum_distance)
        absolute_minimum_distance = min(minimum_distances)
        minimum_distance_door_point = door_points[ minimum_distances.index(absolute_minimum_distance) ]
        self._pivot = minimum_distance_door_point
        return self._pivot
    # The door pivot
    pivot = property(get_pivot, None, None, "The door pivot")

    # Get the polygon where the door is meant to be
    def get_room_polygon (self) -> Optional[Polygon]:
        room = self.room
        if not room:
            return None
        boundary = room.boundary
        if not boundary:
            return None
        return boundary.exterior_polygon

    # Generate a new segment which represents the door open
    # Note that this function is used for display pourposes only
    def get_open_door (self) -> 'Segment':
        segment = self.segment
        pivot = self.pivot
        if not segment or not pivot:
            return None
        return Segment(pivot, pivot + self.direction * segment.length)

    # Get all suitable segments to place the point for this door
    # Get also suitable points (i.e. segment where the door fits in 1 exact point, which is something prone to happen)
    # There are 2 conditions for a segment to be suitable:
    # - The door and its margins fit in the segment
    # - The door does not overlap with other doors (margins may overlap)
    # Note that there is no problem if the door is in contact with more than one parent/children room
    # Corridors are meant to fix these situations
    def find_suitable_regions (self, available_segments : Optional[List[Segment]] = None) -> Tuple[ List['Segment'], List['Point'] ]:
        # Set the minimum length a segment must have in order to fit the door and its margins
        minimum_segment_length = self.margined_width
        # If not available segments are passed then we use the room exterior polygon segments after substracting other doors
        if available_segments == None:
            # The door must have a polygon
            room_polygon = self.get_room_polygon()
            if not room_polygon:
                raise ValueError('Cannot find a suitable region for a door without polygon')
            # Get segments in its polygon which are wide enought for the door
            candidate_segments = [ segment for segment in room_polygon.segments if segment.length >= minimum_segment_length ]
            # Set the doors which must be substracted from the suitable segments (i.e. other doors already set)
            already_set_door_segments = [ door.segment for door in self.room.doors if door != self and door.point ]
            # Substract the neighbour doors from the available segments
            fit_segments = []
            for segment in candidate_segments:
                # Substract already set doors (without their margins) from the current segment
                free_segments = segment.substract_segments(already_set_door_segments)
                new_fit_segments = [ segment for segment in free_segments if segment.length >= minimum_segment_length ]
                fit_segments += new_fit_segments
            # If there is not wide enought room segments we stop here
            # This may happen if the minimum size of the room is not enought to fit the door width including its margins
            if len(fit_segments) == 0:
                raise ValueError('There is not a segment wide enought to fit the door')
            available_segments = fit_segments
        # Get from each suitable segment the region where the door may fit (i.e. the margined suitable segment)
        # Wide segments will provide a suitable segment
        # Exact segments (i.e. same length that the margined door) will provide a suitable point
        suitable_segments = []
        suitable_points = []
        for segment in available_segments:
            if segment.length < minimum_segment_length:
                continue
            # Cut the margins of all available segments for the reminaing segments to be available to store the door point (center)
            elif segment.length == minimum_segment_length:
                suitable_point = segment.get_middle_point()
                suitable_points.append(suitable_point)
            else:
                suitable_segment = segment.get_margined_segment(minimum_segment_length / 2)
                suitable_segments.append(suitable_segment)
        return suitable_segments, suitable_points

        
# A room is a smart boundary that may contain other boundaries with conservative areas and size restrictions
class Room:
    def __init__ (self,
        # A start boundary may be passed. If no boundary is passed it is assigned automatically
        boundary : Optional[Boundary] = None,
        # The 'forced_area' argument stablishes the expected final room area. If no area is passed then the original boundary area will be used
        # Forced area may be a number (absolute value) or a string (percent) with the format 'XX%'
        forced_area : Optional[Union[number, str]] = None,
        # Minimum size in both x and y dimensions
        min_size : Optional[number] = None,
        # Maximum size in at least one dimension
        max_size : Optional[number] = None,
        # Set if the room boundary may me modified
        # DANI: No está implementado todavía
        rigid : bool = False,
        # Corridor size (the minimum size by default)
        corridor_size : Optional[number] = None,
        # The doors to enter the room and start corridors
        doors : Optional[List[Door]] = None,
        # The default inputs for all doors whose inputs are not specified
        # These options are inherited by children rooms whose options are not specified
        door_args : Optional[Dict] = None,
        # The name and color parameters are only representation parameters and they have no effect in the logic
        display : bool = False,
        name : str = 'Unnamed',
        segments_color : str = 'black',
        fill_color : str = 'white',
        # Set other rooms inside this room
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
        self._children = None
        self.children = children
        # Set the boundary
        # If the boundary has been forced then update the display with the initial segments
        if boundary:
            self.boundary = boundary
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
                raise InputError('Forced area has a non-supported format in room ' + name)
        else:
            self.forced_area = self.area
        # If the forced area does not cover the minimum size it makes no sense
        if min_size and self.forced_area and self.forced_area < min_size**2:
            raise InputError('Forced area is not sufficient for the minimum size in room ' + name)
        # Set if the boundary is rigid
        self.rigid = rigid
        # Set size limits
        if min_size or min_size == 0:
            self.min_size = min_size
        else:
            self.min_size = 0
        # Save the input max size
        self.input_max_size = max_size
        # Set the maximium size as infinite by now
        # This may be changed further
        self.max_size = inf
        # Check if the parent has a stablished boundary
        # All the logic changes drastically whether the parent has or not a boundary
        self._parent_has_boundary = False
        # Parent free limit is set by the parent while setting the child boundary
        # Parent free limit is the maximum min size of all parent children but this child
        self.parent_free_limit = None
        # Save input door args
        self.door_args = door_args
        # Save input doors
        # If there is no input doors set a single defualt door
        self.doors = doors if doors else [ Door() ]
        # Set each door room
        for door in self.doors:
            door.room = self
        # Set the room corridor size
        self.corridor_size = corridor_size
        if not self.corridor_size:
            self.corridor_size = self.min_size
        if not self.corridor_size:
            self.corridor_size = self.get_min_min_size()
        # In order to respect the minimum size we are going to use a greater minimum size at the first steps
        # This is done to have margin for performing invasive steps (e.g. corridor and wall thickness build)
        if preventive_min_size_protocol == 0:
            self.preventive_min_size = self.min_size
        elif preventive_min_size_protocol == 1:
            self.preventive_min_size = self.min_size + self.corridor_size / 2
        elif preventive_min_size_protocol == 2:
            self.preventive_min_size = self.min_size + self.corridor_size
        elif preventive_min_size_protocol == 3:
            self.preventive_min_size = self.min_size + self.corridor_size * 2
        else:
            raise SystemExit('Preventive min size protocol ' + str(preventive_min_size_protocol) + ' not defined')
        # Children handling:
        if self.children:
            # Now set some parameters in children
            for child in self.children:
                child._post_init(self)
            # Check areas of all children rooms to do not sum up more than the parent area
            # In addition check if any of the children room has boundary and, if so, check the boundary is inside the parent
            children_area = sum([ child.forced_area for child in self.children ])
            if self.area and children_area > self.area + minimum_resolution:
                raise InputError('Children together require more area than the parent has')
            # In case children do not cover the whole parent area set a new dummy room to cover this free area first
            if children_area < self.forced_area - minimum_resolution:
                remaining_area = self.forced_area - children_area
                dummy_room_name = self.name + ' (free)'
                dummy_room = Room(forced_area=remaining_area, min_size=self.min_size, name=dummy_room_name, fill_color=self.fill_color)
                self.children = [ dummy_room ] + self.children

    # Set a function which sets other initial values once the parent has been set
    # This function is called by the parent once it has been initiated
    def _post_init (self, parent):
        # Check if the parent has a stablished boundary
        self._parent_has_boundary = self.parent.boundary != None
        # Set the absolute forced area in case the area was set as a percent of the parent area
        if self._forced_area_portion:
            if not parent.forced_area:
                raise InputError('Room "' + self.name + '" has a portion area but the parent has not area')
            self.forced_area = parent.forced_area * self._forced_area_portion
        # Calculate the coherent max size, which is the forced area divided by the minimum size
        # i.e. the maximum possible size in case the boundary was the thinest rectangle
        if self.forced_area:
            self.max_size = self.forced_area / self.preventive_min_size
        # In case there was a forced maximum size and it is smaller than the coherent max size, apply it
        if self.input_max_size and self.input_max_size < self.max_size:
            self.max_size = self.input_max_size
        # If the preventive min size is bigger than the maximum size we have to adapt to it
        if self.preventive_min_size > self.max_size:
            self.preventive_min_size = self.max_size
        # In case the parent has a boundary, check this room may fit on it
        if self._parent_has_boundary:
            # Check this room is inside the parent boundary, if they have a predefined boundary
            if self.boundary and self.boundary not in parent.boundary:
                raise InputError('The child room "' + self.name + '" is out of the parent boundary')
            # Check if this rooms minimum size fits in the parent boundary
            if not self.boundary and not parent.does_room_fit(self, force=True):
                raise InputError('The child room "' + self.name + '" minimum size does not fit in the parent boundary')

    def __str__(self):
        return '<Room "' + str(self.name) + '">'

    def __repr__(self):
        return '<Room "' + str(self.name) + '">'

    # Get the children rooms
    def get_children (self):
        return self._children

    # Set the children rooms
    # Update hierarchy
    def set_children (self, children):
        for child in children:
            child.parent = self
        self._children = children

    # The children rooms
    children = property(get_children, set_children, None, "The children rooms")

    # Get the boundary
    # Just return the internal boundary value
    def get_boundary (self):
        return self._boundary

    # Set the boundary
    # Reset the own rects and the parent rects also
    def set_boundary (self, boundary, skip_update_display : bool = False):
        self._boundary = boundary
        self.reset_free_grid()
        if self.parent:
            self.parent.reset_free_grid()
        if not skip_update_display:
            self.update_display()

    # The room boundary
    boundary = property(get_boundary, set_boundary, None, "The room boundary")

    # Get the grid
    # Just return the internal boundary grid value
    def get_grid (self):
        return self._boundary.grid

    # The room boundary
    grid = property(get_grid, None, None, "The room grid")

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
            free_grid = self.grid
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

    # Free space area (read only)
    def get_free_area (self):
        return self.free_grid.area
    free_area = property(get_free_area, None, None, "Free space area (read only)")

    # Set the children boundaries according to the room configuration
    # This function triggers the logic to solve room distributions
    # If the recursive flag is passed then set each child's children boundaries and so on recursively
    # All children rooms must have their boundary fully set before solving the next generation of children
    def set_children_boundaries (self, recursive : bool = False):
        rooms = self.children
        # Sort children rooms by minimum size, with the biggest sizes first
        def sort_by_size (room):
            return room.min_size
        sorted_rooms = sorted( rooms, key=sort_by_size, reverse=True )
        # Set up each room by giving them a position and correct size to match the forced area
        for room in sorted_rooms:
            # Configure the child room to respect the parent free min size limit according to its brothers
            parent_free_limit = max([ other.min_size for other in rooms if other != room ])
            room.parent_free_limit = parent_free_limit
            # If the children has no boundary it must be built
            if not room.boundary:
                if not self.set_child_room_boundary(room):
                    raise SystemExit('Child ' + room.name + ' has failed to be set')

        # Now that all children bondaries are set we must set the doors and the corridor

        # Setup the doors
        # Set the default door arguments in case they were not specified
        if not self.door_args:
            # If we are not the root, inherit parent default door arguments
            if self.parent:
                self.door_args = self.parent.door_args
            # Otherwise, we have to guess the most suitable arguments
            # Make the margined width of all doors equal to the minimum room minimum size
            # Make the width of all doors the 80% of the margined width
            margined_width = self.get_min_min_size()
            width = margined_width * 0.8
            margin = margined_width * 0.1
            self.door_args = {
                'width': width,
                'margin': margin
            }
        # Set each missing door args
        for door in self.doors:
            for arg, value in self.door_args.items():
                current_value = getattr(door, arg)
                if not current_value:
                    setattr(door, arg, value)
        # Set each children missing door args
        for room in rooms:
            for door in room.doors:
                for arg, value in self.door_args.items():
                    current_value = getattr(door, arg)
                    if not current_value:
                        setattr(door, arg, value)

        # Set the corridor
        if len(rooms) > 0:
            self.set_corridor()

        # Set children boundaries recurisvely if the recursive flag was passed
        if recursive:
            for room in rooms:
                room.set_children_boundaries(recursive=True)

    # Check if a room fits in this room according to its minimum size
    # Check free space by default and all space if the argument 'force' is passed
    def does_room_fit (self, room : 'Room', force : bool = False) -> bool:
        size = room.min_size
        # If there is no minimum size it will always fit
        if not size:
            return True
        # Generate fitting rects. If there is at least one then the room fits
        grid = self.grid if force else self.free_grid
        fitting_rects = grid.get_fitting_space(size, size)
        if next(fitting_rects, None):
            return True
        return False

    # Search all maximum rects where the specified room minimum rectangle would fit
    # Search in free space by default and all space if the argument 'force' is passed
    def get_room_fitting_rects (self, room : 'Room', force : bool = False) -> List[Rect]:
        size = room.min_size
        # Get all fitting rects
        grid = self.grid if force else self.free_grid
        fitting_rects = grid.get_fitting_space(size, size)
        return list(fitting_rects)

    # Set up a child room boundary
    def set_child_room_boundary (self, room) -> bool:
        # If the parent (self) has a boundary
        if self.boundary:
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
            # Sort the suitable rects by minimum size
            def sort_by_size(rect):
                return min(rect.get_size())
            sorted_suitable_rects = sorted( suitable_rects, key=sort_by_size )

            # Make a backup in case we have to force since other rooms will be modified
            if forced:
                backup = self.make_children_backup()
                
            # Try to set up the new room in all possible sites until we find one
            # Each 'site' means each corner in each suitable rectangle
            # Check each site to allow other rooms to grow
            sites = [ (corner, rect) for rect in sorted_suitable_rects for corner in rect.get_corners() ]
            previous_initial_boundary = None
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
                # At this point, check if the boundary is equal to a previous tried (and failed) boundary
                # This may happend with the 4 corners of the same rect
                if previous_initial_boundary == room.boundary:
                    continue
                previous_initial_boundary = room.boundary
                # Proceed with the expansion of this child room until it reaches its forced area
                if room.fit_to_required_area():
                    return True
                # If the expansion failed then clean the boundary (i.e. the initial boundary)
                room.boundary = None
                if forced:
                    self.restore_children_backup(backup)
            # If at this point we still have no boundary it means there is no available place to set the perimeter
            return False
        # If the parent (self) has not boundary
        else:
            # There are no boundary restrictions
            # All children will be set in one single step, as the maximum initial boundaries
            # To do so, we must find a suitable corner and space (rect) for the boundary to be set
            # This space will be as long as the double of the room maximum size, for it to be able to expand without problems
            huge_size = room.max_size * 2
            # First of all, find all childrens with already set boundaries
            # They are the only reference to place the new room
            brother_boundaries = [ child.boundary for child in self.children if child.boundary ]
            # If this is the first room then there is no reference at all, so we set the first corner as the (0,0) point
            # Then the "available space" will be a fake rect which is huge enough for this room to set its whole boundary
            if len(brother_boundaries) == 0:
                space = Rect(0, 0, huge_size, huge_size)
                # This is the "real" corner, which include its segments
                initial_point = Point(0,0)
                corner = space.get_corner(initial_point)
                room.boundary = room.set_maximum_initial_boundary(corner, space)
                return True
            # If there are other childs with boundaries already then we must find the best location for the new room
            # It has to be next to its brother rooms, in a point which minimizes the corridor length further
            # If there is one brother only then it makes not sense to find a corridor
            if len(brother_boundaries) == 1:
                space = Rect(-huge_size, 0, 0, huge_size)
                reference_point = Point(0,0)
                corner = space.get_corner(reference_point)
                room.boundary = room.set_maximum_initial_boundary(corner, space)
                return True
            # Otherwise, calculate the "pre-corridor" and find the closest point in the current parent "pre-exterior-polygon"
            # First, get the current corridor points
            corridor_segments = self.set_corridor(0)
            corridor_points = list(set(sum([ segment.points for segment in corridor_segments ], ())))
            # Now get all points in the exterior polygon which may be closer to the corridor
            # Note that these point may not be a corner in the exterior polygon, but a point in the middle of any segment
            # This should be easy to find since the current parent is made of rectangular rooms together
            boundary_segments = sum([ boundary.exterior_polygon.segments for boundary in brother_boundaries ], [])
            exterior_polygon_segments = get_non_overlap_segments(boundary_segments)
            exterior_polygon_segment_points = list(set(sum([ segment.points for segment in exterior_polygon_segments ], ())))
            # Now calculate the actual exterior polygon from the segments we found before
            # This polygon will be useful to know the direction towards we must expand the new child further
            exterior_polygon = Polygon.non_canonical(exterior_polygon_segments)
            # Find the point in the current exterior polygon which is close to any point in the current corridor
            # Note that this point may not be a corner in the exterior polygon
            minimum_distance_point = None
            # If we find a point in the polygon which is already in the corrdior then we are done
            common_point = next(( point for point in exterior_polygon_segment_points if point in corridor_points ), None)
            if common_point:
                minimum_distance_point = common_point
            # Otherwise iterate all point combinations to find the closest ones
            minimum_distance = inf
            for exterior_polygon_segment_point in exterior_polygon_segment_points:
                for corridor_point in corridor_points:
                    distance = exterior_polygon_segment_point.get_distance_to(corridor_point)
                    if distance < minimum_distance:
                        minimum_distance = distance
                        minimum_distance_point = exterior_polygon_segment_point
            # Now find which directions we must expand
            # In case the point is a corner (it will always be an inside corner) we will have both directions already
            corner = exterior_polygon.get_corner(minimum_distance_point)
            if corner:
                direction_a = -corner.segments[0].direction
                direction_b = corner.segments[1].direction
                space_segment_a = Segment(corner + direction_a * huge_size, corner)
                space_segment_b = Segment(corner, corner + direction_b * huge_size)
                corner = Corner(corner.x, corner.y, space_segment_a, space_segment_b)
                space = Rect.from_corner(corner)
                room.boundary = room.set_maximum_initial_boundary(corner, space)
                return True
            # In case the point is in a segment we have only one direction to expand and we can choose the other
            polygon_segment = exterior_polygon.get_border_element(minimum_distance_point)
            direction_a = -exterior_polygon.get_border_inside(polygon_segment)
            direction_b = direction_a.rotate(90)
            space_segment_a = Segment(minimum_distance_point + direction_a * huge_size, minimum_distance_point)
            space_segment_b = Segment(minimum_distance_point, minimum_distance_point + direction_b * huge_size)
            corner = Corner(minimum_distance_point.x, minimum_distance_point.y, space_segment_a, space_segment_b)
            space = Rect.from_corner(corner)
            room.boundary = room.set_maximum_initial_boundary(corner, space)
            return True

    # Build a provisional exterior perimeter from the child boundaries
    # This function is meant to be used only when the parent room has not boundary
    def get_provisional_exterior_polygon (self) -> 'Polygon':
        child_boundary_segments = []
        for child in self.children:
            if not child.boundary:
                continue
            child_boundary_segments += child.boundary.exterior_polygon.segments
        exterior_polygon_segments = get_non_overlap_segments(child_boundary_segments)
        exterior_polygon = Polygon.non_canonical(exterior_polygon_segments)
        return exterior_polygon

    # Build a provisional exterior grid
    # This function is meant to be used only when the parent room has not boundary
    # Create a fake square which contains all the children and then add a margin as long as the distance argument
    # This is used by children to expand outwards
    def generate_exterior_free_grid (self, margin_size : number) -> 'Grid':
        exterior_polygon = self.get_provisional_exterior_polygon()
        exterior_box = exterior_polygon.get_box()
        expanded_box = exterior_box.expand_margins(margin_size)
        exterior_free_boundary = Boundary(Polygon.from_rect(expanded_box), [exterior_polygon])
        return exterior_free_boundary.grid

    # Set the initial room boundary as the maximum possible rectangle
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
        # Try to create the most "squared" possible rectangle
        # Calculate how long would be the side of a perfect square
        square_side_length = sqrt(self.forced_area)
        # Set how long will be the short (restricted) side of the rectangle
        # Then calculate the other side length
        # The limit may come from the square side limit, the sapce limit or the own room limit
        first_side_length = min(square_side_length, x_space, y_space, self.max_size)
        second_side_length = self.forced_area / first_side_length
        # Create the new rect fitting the biggest size in the biggest space and the opposite
        # For each size of the new rect use the calculated size only if it is not longer than the space
        if x_space >= y_space:
            new_x_size = second_side_length
            new_y_size = first_side_length
        else:
            new_x_size = first_side_length
            new_y_size = second_side_length
        maximum_rect = Rect.from_corner(corner, new_x_size, new_y_size)
        return Boundary(Polygon.from_rect(maximum_rect))

    # Set the initial room boundary as the minimum possible rectangle
    # This is used when the initial boundary must be forced over other children rooms
    def set_minimum_initial_boundary (self, corner : Point, space : Rect) -> Boundary:
        minimum_rect = Rect.from_corner(corner, self.preventive_min_size, self.preventive_min_size)
        return Boundary(Polygon.from_rect(minimum_rect))

    # Using all boundaries, calculate which segments make the shortest path to connect parent doors and all children rooms
    # Missing children doors are also set during this process and they are placed according to make the corridor as short as possible
    # Finally, the area around the selected corridor segments is claimed to build the actual corridor
    # The protocol allows to select how far you want the logic to go (claiming the area is not reversible)
    # 0 - Find the corridor segments and return the corridor segments
    # 1 - Place the children doors (thus possibly expanding a bit the corridor segments) and return the corridor segments
    # 2 - Claim the corridor area, return nothing
    # Note that this function may be used several times to find a "pre-corridor" before definitely setting it
    def set_corridor (self, protocol : int = 2) -> Optional[List['Segment']]:

        # Set the corridor nodes
        # This is a preprocessing step which is required to find the shortest corridor and other similar calculations
        
        # Use both self boundary and all children room boundary segments
        available_segments = []
        door_points = []
        # Split self boundaries at the doors
        if self.boundary:
            door_points = [ door.point for door in self.doors ]
            for self_segment in self.boundary.segments:
                available_segments += self_segment.split_at_points(door_points)
        # Add all children segments
        for child in self.children:
            if not child.boundary:
                continue
            available_segments += child.boundary.segments
        # Split and merge available segments according to overlaps
        splitted_segments = []
        for available_segment, other_segments in otherwise(available_segments):
            # For each segment, find all cutting point (overlaps are implicit) with the rest of segments
            other_points = []
            for other_segment in other_segments:
                other_points += other_segment.points
            splitted_segments += available_segment.split_at_points(other_points)
        # Remove duplicates
        # At this point there should be no overlaps between splitted segments  
        splitted_segments = list(set(splitted_segments))
        # Now define "nodes"
        # Nodes are points between splitted segments (with no duplicates)
        # Each node may have from 2 up to 4 connected segments
        # Each node may have from 1 up to 4 contact rooms
        # First find all connected segments as nodes are defined
        # Create a dictionary with the points as keys
        nodes = {}
        for segment in splitted_segments:
            points = segment.points
            for point in points:
                current_node = nodes.get(point, None)
                if current_node:
                    current_node['connected_segments'].append(segment)
                else:
                    nodes[point] = {'connected_segments': [segment]}
        # Find which rooms are in contact to each node and if nodes are in the exterior boundary
        # Find also which nodes are doors
        exterior_polygon = self.boundary.exterior_polygon if self.boundary else self.get_provisional_exterior_polygon()
        for node_point, node_data in nodes.items():
            rooms = []
            for child in self.children:
                if not child.boundary:
                    continue
                if node_point in child.boundary.exterior_polygon:
                    rooms.append(child)
            if self.boundary and node_point in exterior_polygon:
                rooms.append(self)
            node_data['rooms'] = rooms
            node_data['is_door'] = node_point in door_points
        # Now find "non-redundant" nodes and the "paths" between them
        # Redundant nodes are those whose contact rooms are already included in all connected nodes
        # Knwoing this is useful when we are expanding our corridor since a reundant node will never solve the puzzle
        # When we expand thorugh redundant nodes we can claim several nodes until we find a non-redundant node
        # Note that nodes with 2 connected segments will always be redundant while others will be non-redundant
        # IMPORTANT: Doors are also non-redundant nodes
        for node_data in nodes.values():
            node_data['is_redundant'] = len(node_data['connected_segments']) == 2 and not node_data['is_door']
        # Then a path between non-redundant nodes is a list of segments, which are connected by redundant nodes
        for node_point, node_data in nodes.items():
            if node_data['is_redundant']:
                continue
            # Find all current node paths to other non-redundant nodes
            paths = []
            # Save the non-redundant node each path is leading to
            path_nodes = []
            for starting_segment in node_data['connected_segments']:
                last_segment = starting_segment
                last_point = next(point for point in last_segment.points if point != node_point)
                last_node = nodes[last_point]
                path = [ last_segment ]
                # Keep moving forward until we find a non-redundant node
                while last_node['is_redundant']:
                    # If the node is redundant then there will be always only 2 connected segments and one of them will be the last
                    last_segment = next( segment for segment in last_node['connected_segments'] if segment != last_segment )
                    last_point = next(point for point in last_segment.points if point != last_point)
                    last_node = nodes[last_point]
                    path.append(last_segment)
                    if last_point == node_point:
                        print(self.name)
                        raise ValueError('This should never happen. Do we have one room only?')
                # Once we have the reached a non-redundant node save the accumulated path and the node itself
                # In case the point is already in the list it means we have two paths for the same node
                # In this case, get the shortest path
                if last_point in path_nodes:
                    index = path_nodes.index(last_point)
                    previous_path = paths[index]
                    if get_path_length(path) < get_path_length(previous_path):
                        paths[index] = path
                else:
                    path_nodes.append(last_point)
                    paths.append(path)
            # Add paths and path nodes to the node data
            node_data['paths'] = paths
            node_data['path_nodes'] = path_nodes

        # ------------------------------------------------------------------------------------------------------------------------------

        # Now find all possible path combinations until we cover all doors and rooms
        # Set the rooms which must be reached by the corridor
        required_children_rooms = [ child for child in self.children if child.boundary and len(child.doors) > 0 ]
        required_rooms = ([ self ] if self.boundary and len(self.doors) > 0 else []) + required_children_rooms
        # If there are less than 2 required rooms then it makes not sense finding a corridor
        if len(required_rooms) < 2:
            raise ValueError('Trying to find a corridor between less than 2 rooms')
        # Set the doors which must be reached by the corridor
        # Doors may have a point (i.e. they are already set) or not
        # Doors which have a point must be reached by the corridor
        # Doors which do not have a point must be set after the corridor and then the corrdior must be expanded to cover them
        doors = sum([ room.doors for room in required_rooms ], [])
        already_set_doors = []
        unset_doors = []
        for door in doors:
            if door.point:
                already_set_doors.append(door)
            else:
                unset_doors.append(door)
        current_corridor = None
        current_corridor_length = None
        current_corridor_nodes = None
        # Trak which combinations of path nodes we have tried allready
        # Combinations of path nodes are equivalent to combinations of paths, but easier to compare
        # This way we do not analyze the same corridor multiple times
        already_covered_path_nodes = []
        # Set a function to generate corridors by recuersively joining node paths
        def get_following_paths (
            current_path : list,
            current_path_nodes : list,
            available_paths : list,
            available_path_nodes : list,
            current_rooms : set
        ):
            nonlocal current_corridor
            nonlocal current_corridor_length
            nonlocal current_corridor_nodes
            nonlocal already_covered_path_nodes
            # Check if the current path nodes have been covered already and stop here if so
            current_path_nodes_set = set(current_path_nodes)
            if any( current_path_nodes_set == set(path_nodes) for path_nodes in already_covered_path_nodes ):
                return
            # Add current path nodes to the covered list in order to avoid repeating this path further
            already_covered_path_nodes.append(current_path_nodes)
            # Try to expand the current corridor using all available paths
            for i, next_path in enumerate(available_paths):
                # The following node will be the other next path's node
                following_node = available_path_nodes[i]
                # Get the corresponding node data
                following_node_data = nodes[following_node]
                # Get the remaining available paths/nodes after substracting the current next path
                following_available_paths = available_paths[0:i] + available_paths[i+1:]
                following_available_path_nodes = available_path_nodes[0:i] + available_path_nodes[i+1:]
                # Get the following node paths which are not already included in the current path and its nodes
                following_node_paths = following_node_data['paths']
                following_node_path_nodes = following_node_data['path_nodes']
                # Add the following node paths/nodes to the remaning available paths/nodes
                # Then we get the available paths/nodes for the following path
                for j, path in enumerate(following_node_paths):
                    node = following_node_path_nodes[j]
                    # If the node is already in the current path nodes list then we skip it
                    if node in current_path_nodes:
                        continue
                    # In case the point is already in the list of available nodes it means we have two paths for the same node
                    # In this case, get the shortest path
                    if node in following_available_path_nodes:
                        index = following_available_path_nodes.index(node)
                        previous_path = following_available_paths[index]
                        if get_path_length(path) < get_path_length(previous_path):
                            following_available_paths[index] = path
                    # Otherwise, add the current new available path and node to the lists
                    else:
                        following_available_paths.append(path)
                        following_available_path_nodes.append(node)
                # Get the new following path after adding the last path while getting the next node
                following_path = current_path + next_path
                following_path_nodes = current_path_nodes + [ following_node ]
                # Get the following path covered rooms
                following_rooms = current_rooms.union(set(following_node_data['rooms']))
                # If following path includes all rooms then it is a candidate to be the corridor
                contains_all_rooms = all(room in following_rooms for room in required_rooms)
                contains_all_doors = all(door.point in following_path_nodes for door in already_set_doors)
                if contains_all_rooms and contains_all_doors:
                    # Check if this path is shorter than the current corridor
                    # The shorter path will remain as the current corridor
                    following_path_length = get_path_length(following_path)
                    if not current_corridor or following_path_length < current_corridor_length:
                        current_corridor = following_path
                        current_corridor_length = following_path_length
                        current_corridor_nodes = following_path_nodes
                    continue
                # If the follwoing path does not cover all rooms yet then keep expanding it
                get_following_paths(
                    following_path,
                    following_path_nodes,
                    following_available_paths,
                    following_available_path_nodes,
                    following_rooms
                )
        # Call the function starting to solve the corridor at the first set door in case we have doors
        # The door is set when it has at least a point
        first_set_door = self.doors and next((door for door in self.doors if door.point), None)
        if first_set_door:
            start_point = self.doors[0].point
            start_node = nodes[start_point]
            start_rooms = set(start_node['rooms'])
            start_path = []
            start_path_points = [ start_point ]
            start_available_paths = start_node['paths']
            start_available_path_nodes = start_node['path_nodes']
            get_following_paths(
                start_path,
                start_path_points,
                start_available_paths,
                start_available_path_nodes,
                start_rooms
            )
        # Otherwise, there is no node which we know for sure it will be part from the corridor
        # We must solve the corridor several times starting from diferent nodes
        # We start on each non-redundant node from the room with less non-redundant nodes
        else:
            all_room_ocurrences = sum([ node_data['rooms'] for node_data in nodes.values() ], [])
            rooms = list(set(all_room_ocurrences))
            node_room_counts = { room: all_room_ocurrences.count(room) for room in rooms  }
            room_with_less_nodes = min(node_room_counts, key=node_room_counts.get)
            nodes_from_room_with_less_nodes = [
                node_point for node_point, node_data in nodes.items() if node_data['is_redundant'] == False and room_with_less_nodes in node_data['rooms']
            ]
            for node in nodes_from_room_with_less_nodes:
                start_point = node
                start_node = nodes[start_point]
                start_rooms = set(start_node['rooms'])
                start_path = []
                start_path_points = [ start_point ]
                start_available_paths = start_node['paths']
                start_available_path_nodes = start_node['path_nodes']
                get_following_paths(
                    start_path,
                    start_path_points,
                    start_available_paths,
                    start_available_path_nodes,
                    start_rooms
                )

        # Display the current corridor
        elements_to_display = current_corridor
        for segment in elements_to_display:
            segment.color = 'red'
        self.update_display(extra=elements_to_display)

        # ------------------------------------------------------------------------------------------------------------------------------

        if protocol < 1:
            return current_corridor

        # Now we have to set doors which are not set already in children rooms and expand the current corridor to cover these doors if needed

        # Set a function to find the closest point to the corridor in a list of suitable segments to place a door
        # Then find the path to this point and add it to the current corridor
        # This is used to place a door which can not be placed in the current corridor because it does not fit anywhere
        def expand_corridor_to_place_door (door : 'Door', suitable_segments : List[Segment], suitable_points : List[Point]):
            nonlocal current_corridor
            nonlocal current_corridor_nodes
            nonlocal nodes
            # Get the corridor non redundant nodes in contact with the door room
            room = door.room
            room_nodes = [ node for node in current_corridor_nodes if not nodes[node]['is_redundant'] and room in nodes[node]['rooms'] ]
            # The closest point will be always one of the ends of one of the suitable segments
            suitable_points += list(set(sum([ list(segment.points) for segment in suitable_segments ], [])))
            # In order to find the shortest path, first, get the room exterior perimeter
            perimetral_segments = room.boundary.exterior_polygon.segments
            # Now split those segments by both nodes and suitable points
            split_points = room_nodes + suitable_points
            splitted_segments = sum([ list(segment.split_at_points(split_points)) for segment in perimetral_segments ], [])
            # Then build paths using these splitted segments and keep the shortest path
            # Start at each node and keep joining connected segments until a suitable point is found
            current_path = None
            current_path_length = None
            for node in room_nodes:
                # Note that 2 paths may be built from each node
                start_segment_generator = ( segment for segment in splitted_segments if segment.has_point(node) )
                start_segments = [ next(start_segment_generator), next(start_segment_generator) ]
                for start_segment in start_segments:
                    path = [ start_segment ]
                    current_segment = start_segment
                    current_point = start_segment.get_other_point(node)
                    while current_point not in split_points:
                        current_segment = next((segment for segment in splitted_segments if current_point in segment and segment != current_segment ))
                        path.append(current_segment)
                        current_point = current_segment.get_other_point(current_point)
                    # If another node is found while joining segments then skip this path
                    if current_point in room_nodes:
                        continue
                    # If this path is shorter than the current path then keep this path as the current
                    path_length = get_path_length(path)
                    if not current_path or path_length < current_path_length:
                        current_path_length = path_length
                        # Edit the last segment in the path in order to cover the whole door with the corridor
                        previous_point = current_segment.get_other_point(current_point)
                        current_direction = (previous_point + current_point).normalized()
                        further_point = current_point + current_direction * door.margined_width / 2
                        last_segment = Segment(previous_point, further_point)
                        path[-1] = last_segment
                        # Update the current path and the door point
                        current_path = path
                        door.point = current_point
            # Now add the current path to the whole corridor
            current_corridor += current_path
            # Add the node and a fake node data in order to make the corridor further expansable at this new point
            # WARNING: The fake node data coontains way less data than other nodes and it is not fully functional
            current_corridor_nodes.append(further_point)
            nodes[further_point] = {'is_redundant': False, 'rooms': [room]}

        # For each door whose point is not yet set, set its point now
        for door in unset_doors:
            # Get the corridor non redundant nodes in contact with the door room
            room = door.room
            room_nodes = [ node for node in current_corridor_nodes if not nodes[node]['is_redundant'] and room in nodes[node]['rooms'] ]
            # Check if any node is suitable to fit the door
            # Note that there is no problem if the door is in contact with more than one parent/children room
            # Corridors are meant to fix these situations
            suitable_segments, suitable_points = door.find_suitable_regions()
            final_point = None
            for node in room_nodes:
                if node in suitable_points:
                    final_point = node
                    break
                suitable_segment = next((segment for segment in suitable_segments if node in segment), None)
                if suitable_segment:
                    final_point = node
                    break
            if final_point:
                door.point = final_point
                continue
            # If not, we must check if there is any path between those nodes and try to fit the door in this path
            for segment in current_corridor:
                suitable_point = next((point for point in suitable_points if point in segment), None)
                if suitable_point:
                    final_point = suitable_point
                    break
                for suitable_segment in suitable_segments:
                    overlap = segment.get_overlap_segment(suitable_segment)
                    if overlap:
                        final_point = overlap.get_random_point()
                        break
                if final_point:
                    break
            if final_point:
                door.point = final_point
                continue
            # If not, we must find a point as close as posible to the corridor nodes
            # We must also find the path from the closest node to this closest point
            expand_corridor_to_place_door(door, suitable_segments, suitable_points)

        # Display the final corridor
        elements_to_display = current_corridor
        for segment in elements_to_display:
            segment.color = 'red'
        self.update_display(extra=elements_to_display)

        # ------------------------------------------------------------------------------------------------------------------------------

        if protocol < 2:
            return current_corridor

        # Build the corridor by claiming area around the corridor path
        # To do so we must set a boundary
        # With the current implementation there should never be cyclic (closed) corridors
        # For this reason the boundary should only have an exterior polygon
        # However if this happens in the future there should be no problem since boundaries have interior polygons as well

        # Set the corridor size (width)
        corridor_size = self.corridor_size
        half_corridor_size = corridor_size / 2

        # Split corridor segments at the parent exterior perimter inside corners
        # This is to avoid segments which are partially in the parent to get the full offset* in the whole segment
        # * Offset means the wall displacement when claiming the corridor area, it is explained below
        cut_points = [ corner for corner in exterior_polygon.corners if corner.inside ]
        for segment in current_corridor:
            cutted_segments = list(segment.split_at_points(cut_points))
            if len(cutted_segments) > 1:
                current_corridor.remove(segment)
                current_corridor += cutted_segments
        # Set the children doors
        # These doors are susceptible of beeing moved while the corridor is build
        children_doors = sum([ child.doors for child in required_children_rooms ], [])

        # Set a function to generate the corridor boundary
        # This process is wrapped in a function because we may have to change the corridor and redo the boundary further
        # e.g. a door can not be relocated in the boundary so it must be relocated now and the corridor will change
        def generate_corridor_boundary ():
            # Generate data for each point between segments (similar to nodes) by recording the connected segments
            # Generate data for each segment by generating 2 boundary lines in the corridor
            point_connected_segments = {}
            segment_lines = {}
            for segment in current_corridor:
                # Get the points connected segments
                points = segment.points
                for point in points:
                    connected_segments = point_connected_segments.get(point, None)
                    if connected_segments:
                        connected_segments.append(segment)
                    else:
                        point_connected_segments[point] = [segment]
                # Get the children doors inside the current segment
                segment_doors = [ door for door in children_doors if door.point in segment and door.segment.same_line_as(segment) ]
                # Get the segment inside direction, which is the direction it will be offseted in case the segment is too close to the parent boundary
                # get also how much it will be offseted (the distance). It may range from 0 to half the corridor size.
                inside_direction = None
                inside_distance = 0
                # If the segment is overlapping one fo the segments in the parent boundary then we have to offset its expansion
                # One of the boundary lines will remain in the same segment line and the other will be translated the whole corridor size intead of a half
                overlap_segment = next(exterior_polygon.get_segment_overlap_segments(segment), None)
                if overlap_segment:
                    inside_direction = exterior_polygon.get_border_inside(segment)
                    inside_distance = half_corridor_size
                # If there is an inside corner from the parent boundary close enought to the segment then we have to offset it as well
                for corner in exterior_polygon.get_inside_corners():
                    line_closer_point = segment.line.get_closer_point(corner)
                    # If the corner is not in the "perpendicular range" then there is no conflict
                    if line_closer_point not in segment:
                        continue
                    # If the corner is far enough then there is no conflict
                    corner_distance = line_closer_point.get_distance_to(corner)
                    if corner_distance > half_corridor_size:
                        continue
                    # We have a conflict here
                    # Find the direction we must push the corridor to solve the conflict
                    if corner_distance == 0:
                        affecting_segment = next(( seg for seg in corner.segments if segment.is_paralel_to(seg) ), None)
                        if not affecting_segment:
                            raise ValueError('This should not happen. Does the exteropr polygon has diagonal segments?')
                        corner_inside_direction = exterior_polygon.get_border_inside(affecting_segment)
                    else:
                        corner_inside_direction = (corner + line_closer_point).normalized()
                    # If there was a previous inside direction and it is not matching the current one it must be the opposite
                    # This means we must offset corridors in both perpendicular directions at the same time, which is impossible
                    # This may happen in very elaborate situations, and if it happens it may have a solution, but it is too much work
                    # I rather not to support this situation
                    if inside_direction and inside_direction != corner_inside_direction:
                        raise ValueError('Can not push corridor in both directions at the same time (see segment ' + str(segment) + ')')
                    inside_direction = corner_inside_direction
                    # Calculate how much we must push the line inside to avoid the conflict corner
                    # The highest inside distance stays
                    corner_inside_distance = half_corridor_size - corner_distance
                    if inside_distance < corner_inside_distance:
                        inside_distance = corner_inside_distance

                # Set the segment boundary lines
                # Each segment will have 2 lines: one on each side
                # We call these sides as clockwise and counter-clockwise sides
                segment_direction = segment.direction
                clockwise_direction = segment_direction.rotate(90)
                counterclockwise_direction = segment_direction.rotate(-90)
                # Use segment 'a' point as a reference to set each line point, but it could be the segment 'b' point as well
                if inside_direction:
                    # Push the line points
                    inside_offset = half_corridor_size + inside_distance
                    outside_offset = half_corridor_size - inside_distance
                    if inside_direction == clockwise_direction:
                        clockwise_point = segment.a + clockwise_direction * inside_offset
                        counterclockwise_point = segment.a + counterclockwise_direction * outside_offset
                    elif inside_direction == counterclockwise_direction:
                        clockwise_point = segment.a + clockwise_direction * outside_offset
                        counterclockwise_point = segment.a + counterclockwise_direction * inside_offset
                    else:
                        raise ValueError('Impossible inside direction when solving corridor boundary')
                    # Push the doors
                    for door in segment_doors:
                        if door.direction == inside_direction:
                            setattr(door, 'offset', inside_offset)
                        elif door.direction == -inside_direction:
                            setattr(door, 'offset', outside_offset)
                        else:
                            raise ValueError('Impossible door direction when solving corridor boundary')
                else:
                    clockwise_point = segment.a + clockwise_direction * half_corridor_size
                    counterclockwise_point = segment.a + counterclockwise_direction * half_corridor_size
                    for door in segment_doors:
                        setattr(door, 'offset', half_corridor_size)
                clockwise_line = Line(clockwise_point, segment_direction)
                counterclockwise_line = Line(counterclockwise_point, segment_direction)
                # Save both lines in a data dictionary
                data = { 'clockwise': clockwise_line, 'counterclockwise': counterclockwise_line }
                # Save the data dictionary inside another dictionary where the key is the segment itself
                segment_lines[segment] = data

            # Join all segment lines together in a dictionary where lines are the keys
            # Each line value will be a list of intersections which will be set empty at this moment
            # At the end of the next step each line must have exactly 2 intersections
            # Note that duplicated lines will remain as a single key
            line_intersections = {}
            for lines in segment_lines.values():
                clockwise_line = lines['clockwise']
                counterclockwise_line = lines['counterclockwise']
                for line in [ clockwise_line, counterclockwise_line ]:
                    line_intersections[line] = []

            # Now for each segment in the corridor get the segments of the boundary
            # First, segments must be built by finding the intersection point between segment lines
            boundary_segments = []
            for point, connected_segments in point_connected_segments.items():
                # In case we have only 1 connected segment it means this is a death end of the corridor
                # In this case we generate a new segment perpendicular to the only segment and which crosses the point itself
                # This segment will be generated from a line which intersects both of the boundary lines in the only segment
                if len(connected_segments) == 1:
                    segment = connected_segments[0]
                    # Using the current point as the reference point of view
                    is_segment_pointing_outside = point == segment.a
                    lines = segment_lines[segment]
                    clockwise_line = lines['clockwise'] if is_segment_pointing_outside else lines['counterclockwise']
                    counterclockwise_line = lines['counterclockwise'] if is_segment_pointing_outside else lines['clockwise']
                    perpendicular_vector = clockwise_line.vector.rotate(90)
                    perpendicular_line = Line(point, perpendicular_vector)
                    clockwise_line_intersection = clockwise_line.get_intersection_point(perpendicular_line)
                    counterclockwise_line_intersection = counterclockwise_line.get_intersection_point(perpendicular_line)
                    line_intersections[clockwise_line].append((clockwise_line_intersection, segment))
                    line_intersections[counterclockwise_line].append((counterclockwise_line_intersection, segment))
                    new_segment = Segment(clockwise_line_intersection, counterclockwise_line_intersection)
                    boundary_segments.append(new_segment)
                    continue
                # Using the current point as the reference point of view:
                # Sort segments according to their order around the point
                reference_vector = Vector(0,1) # This could be any vector
                def get_reference_angle (segment : 'Segment') -> number:
                    pointing_outside_vector = segment.vector if point == segment.a else -segment.vector
                    return reference_vector.get_angle_with(pointing_outside_vector)
                sorted_connected_segments = sorted(connected_segments, key=get_reference_angle)
                # print('SORTED ' + str(point))
                # print([ segment.vector for segment in sorted_connected_segments ])

                # For each pair of segments, there is a pair of boundary lines (one from each segment) which must intersect
                # As an exception, if segments are paralel, we must check if lines are the same line
                # In this case there the intersection point will be the middle point (bot segments will be merged further)
                # Otherwise, we will have to add a perpendicular segment to intercept both lines to close the corridor at some point
                for current, nextone in pairwise(sorted_connected_segments, retro=True):
                    # Using the current point as the reference point of view:
                    # Get the intersection point between the line in the clockwise side of the current segment and
                    #   the line in the counterclockwise side of the next one
                    is_current_pointing_outside = point == current.a
                    current_side = 'clockwise' if is_current_pointing_outside else 'counterclockwise'
                    current_lines = segment_lines[current]
                    current_clockwise_line = current_lines[current_side]
                    is_nextone_pointing_outside = point == nextone.a
                    nextone_side = 'counterclockwise' if is_nextone_pointing_outside else 'clockwise'
                    nextone_lines = segment_lines[nextone]
                    nextone_counterclockwise_line = nextone_lines[nextone_side]
                    intersection = current_clockwise_line.get_intersection_point(nextone_counterclockwise_line)
                    # If there is an intersection then save this point as an intersection for both lines
                    if intersection:
                        line_intersections[current_clockwise_line].append((intersection, current))
                        line_intersections[nextone_counterclockwise_line].append((intersection, nextone))
                    # If there is no intersection it means lines are paralel
                    else:
                        # If they are the same line then there is no intersection at this point
                        if current_clockwise_line == nextone_counterclockwise_line:
                            continue
                        # Otherwise, generate a new paralel segment which cuts both lines thus closing the corridor
                        perpendicular_vector = current_clockwise_line.vector.rotate(90)
                        # Find which line is closer to the point and find if this line 
                        current_line_distance = current_clockwise_line.get_distance_to(point)
                        nextone_line_distance = nextone_counterclockwise_line.get_distance_to(point)
                        current_is_closer = current_line_distance < nextone_line_distance
                        closer_line = current_clockwise_line if current_is_closer else nextone_counterclockwise_line
                        # If the line is intersecting with the point itself then the new segment must start at the point itself
                        # (Note that we are in a corner of the parent exterior boundary)
                        if point in closer_line:
                            new_line = Line(point, perpendicular_vector)
                        # Otherwise, the new line must be pushed in one direction in order to make space for the corridor
                        # The direction of the push must be through where is the segment whom the closer line comes from
                        else:
                            if current_is_closer:
                                offset_direction = current.direction
                                if not is_current_pointing_outside:
                                    offset_direction = -offset_direction
                            else:
                                offset_direction = nextone.direction
                                if not is_nextone_pointing_outside:
                                    offset_direction = -offset_direction
                            offset = offset_direction * corridor_size
                            offset_point = point + offset
                            new_line = Line(offset_point, perpendicular_vector)
                        # Find the intersection points with each line and create a new segment from both interactions
                        current_line_intersection = current_clockwise_line.get_intersection_point(new_line)
                        nextone_line_intersection = nextone_counterclockwise_line.get_intersection_point(new_line)
                        line_intersections[current_clockwise_line].append((current_line_intersection, current))
                        line_intersections[nextone_counterclockwise_line].append((nextone_line_intersection, nextone))
                        new_segment = Segment(current_line_intersection, nextone_line_intersection)
                        boundary_segments.append(new_segment)

            # Build segments out of all found intersection points
            for line, intersections in line_intersections.items():
                # In the majority of the ocasions we will have 2 intersections per line
                if len(intersections) == 2:
                    intersection_points = [ intersection[0] for intersection in intersections ]
                    # It may happen in a fe ocassions that a line has the same intersection twice
                    # This happens when a full offset corridor eats a perpendicular small corridor with the exact length of the corridor size
                    if intersection_points[0] == intersection_points[1]:
                        continue
                    new_segment = Segment(intersection_points[0], intersection_points[1])
                    boundary_segments.append(new_segment)
                    continue
                # It may happen in a few ocassions that a line has more than 2 intersections (e.g. 4)
                # This happens when corridor boundaries overlap in a 'U' shaped corridor
                original_segments = list(set([ intersection[1] for intersection in intersections ]))
                new_segments = []
                for original_segment in original_segments:
                    intersection_points = [ intersection[0] for intersection in intersections if intersection[1] == original_segment ]
                    # At this point there should be always 2 intersections only
                    if len(intersection_points) != 2:
                        # Display the corridor
                        elements_to_display = current_corridor
                        for segment in elements_to_display:
                            segment.color = 'red'
                        self.update_display(extra=elements_to_display)
                        print('Line ' + str(line) + ' from the original segment ' + str(original_segment))
                        raise ValueError('Different number of intersections than 2: ' + str(intersection_points))
                    if intersection_points[0] == intersection_points[1]:
                        continue
                    new_segment = Segment(intersection_points[0], intersection_points[1])
                    new_segments.append(new_segment)
                # In case the new segments overlap, remove the overlapping region
                # Note that this step is performed in order to make the corridor boundary clean but,
                #   if not done, the segment would be removed anyway since the polygon is set through the non canonical pathway
                final_segments = get_line_non_overlap_segments(new_segments)
                boundary_segments += final_segments

            # Display the corridor boundary
            elements_to_display = boundary_segments
            for segment in elements_to_display:
                segment.color = 'red'
            self.update_display(extra=elements_to_display)

            # Generate a boundary from the previous segments
            # DANI: Con la implementación actual debería haber siempre un único polígono
            # DANI: Es posible que esto cambie en el futuro
            polygons = list(connect_segments(boundary_segments))
            if len(polygons) > 1:
                raise ValueError('There should be only 1 polygon at this point')
            corridor_boundary = Boundary(polygons[0])
            return corridor_boundary

        # Run the function above to generate the corridor
        corridor_boundary = generate_corridor_boundary()
        corridor_polygon = corridor_boundary.exterior_polygon

        # Check doors to be able to relocate into this new boundary before proceeding to claim it
        for door in children_doors:
            # Check the door to fit and, in case it does not fit, exapand the corridor (this usually does not happend)
            # Repeat this process until the door fits (if it happens, it happens usually 1 time to be solved)
            while True:
                # Get those corridor segments which are not in the parent exterior perimeter
                # WARNING: This is important since these segments in the corridor may overlap the current door room polygon
                # WARNING: However, they will not exist once the corridor has been set (it is hard to imagine if you don't see it)
                avaliable_corridor_segments = corridor_polygon.get_polygon_non_overlap_segments(exterior_polygon)
                common_segments = door.room.grid.get_segments_overlap_segments(avaliable_corridor_segments)
                suitable_segments, suitable_points = door.find_suitable_regions(common_segments)
                if len(suitable_segments) > 0 or len(suitable_points) > 0:
                    break
                # In case there is not available space to relocate the door in the boundary we must relocate the door now
                # Then we will expand the corridor to cover the door and remake the boundary
                print('WARNING: A door from ' + door.room.name + ' could not be relocated to the corridor boundary so it must be expanded')
                # Get the segments of the door room after substracting the corridor from them
                uncommon_segments = door.get_room_polygon().get_non_overlap_segments(current_corridor)
                suitable_segments, suitable_points = door.find_suitable_regions(uncommon_segments)
                if len(suitable_segments) == 0 and len(suitable_points) == 0:
                    raise ValueError('There is no place to relocate a door from ' + door.room.name)
                # Expand the corridor to place the door somewhere in these segments
                expand_corridor_to_place_door(door, suitable_segments, suitable_points)
                # Remake the boundary now that the corridor has been expaned
                corridor_boundary = generate_corridor_boundary()
                # WARNING: We must get the corridor polygon again
                corridor_polygon = corridor_boundary.exterior_polygon
            # Save the already found suitable segments and points in case we need them further
            setattr(door, 'suitable_segments', suitable_segments)
            setattr(door, 'suitable_points', suitable_points)


        # Substract this region from the rest of rooms and claim it back for the parent
        # Note that, at this point, other rooms will not expand to compensate the lost at this time
        for child in self.children:
            if not child.truncate([corridor_boundary], force=True, skip_update_display=True):
                raise ValueError('The space required by the corridor cannot be claimed from ' + child.name)
        # Set the corridor as a new independen room which will occupy the new freed space
        # This is useful to prevent this space to be claimed back if we further expand the rooms
        corridor_room = Room(boundary=corridor_boundary, name=self.name + ' corridor', fill_color=self.fill_color)
        corridor_room.parent = self
        self.children.append(corridor_room)

        # Show the corridor area
        self.update_display()

        # Now, if the parent has no boundary, we expand child rooms to compensate for their area loss
        if not self.boundary:
            for child in self.children:
                if not child.fit_to_required_area():
                    raise ValueError(child.name + ' failed to fit to required area after corridor area truncation')
        
        # Finally relocate doors to the new corridor boundary
        for door in doors:
            offset = getattr(door, 'offset', None)
            if not offset:
                continue
            # Try to simply push the door as much as the corrdior was expanded
            # If it does not work then we relocate the door in any other available segment/point
            suitable_segments = getattr(door, 'suitable_segments', None)
            suitable_points = getattr(door, 'suitable_points', None)
            push = door.direction * offset
            new_point = door.point + push
            if new_point in suitable_points or any(new_point in segment for segment in suitable_segments):
                door.point = new_point
                continue
            # If it is not in the corridor boundary then we must relocate it
            # There should always be space for the door a this point since we have already checked this
            suitable_segment_points = sum([ list(segment.points) for segment in suitable_segments ], [])
            available_points = suitable_points + suitable_segment_points
            door.point = random.choice(available_points)

        # Show the corrdior area and the relocated doors
        self.update_display()

    # Get the minimum of all minimum sizes
    # Get to the the root and the check all children min sizes recuersively in order to get the minimum
    # This function is meant to be used only once by the root, so its value is not stored
    # WARNING: 0 values are removed
    def get_min_min_size (self) -> number:
        root = self.get_root_room()
        all_rooms = root.get_rooms_recuersive()
        all_min_sizes = [ room.min_size for room in all_rooms if room.min_size > 0 ]
        if len(all_min_sizes) == 0:
            return 0
        return min(all_min_sizes)

    # Calculate how much area we need to expand
    def get_required_area (self):
        return resolute(self.forced_area - self.area)

    # Expand or contract this room until it reaches the forced area
    # In case it is not able to fit at some point recover the original situation
    # Restricted segments are segments which must remain as are
    def fit_to_required_area (self, restricted_segments : list = []) -> bool:
        # Make a boundary backup of this room and all its brothers
        backup = self.parent.make_children_backup()
        # Calculate how much area we need to expand
        required_area = self.get_required_area()
        # Keep expanding or contracting until the current room reaches the desired area
        # Check the required area is big enought to be meaningfull according to the resolution
        # i.e. check if expanding a segment with the minimum length would make it move more than the minimum resolution
        # If not, then we have finished the fitting
        # The resolution is multiplied by the minimum size since area error will always be bigger than length error
        # DANI: Esto último es nuevo, hay que ver que tal
        # DANI: En realidad este problema no tiene solución fácil
        # DANI: Es mejor con la multiplicación que sin ella. Le da más flexibilidad a la resolución del puzle acelerando así el proceso
        # DANI: Es muy peligrosos que queden espacios libres sin reclamar (cuando no tenga que haberlos)
        # DANI: Pero una resolución pequeña no hará que no queden espacios libres, simplemente hará que esos espacios sean muy pequeños
        while abs(required_area / self.max_size) >= minimum_resolution * self.max_size:
            # If the required are is positive it means we must expand our room
            if required_area > 0:
                # Set a function to supervise if each expansion step is succesful or not
                def expand_step () -> bool:
                    # Get the most suitable frontier to expand and try to expand it
                    # If the expansions fails, try with the next one
                    for frontier, loan_permission in self.get_best_frontiers(restricted_segments):
                        if self.expand_frontier(frontier, required_area, loan_permission):
                            return True
                    return False
                # Expand
                if expand_step():
                    # Refresh how much area we need to expand and go for the next step
                    required_area = self.get_required_area()
                    continue
            # If the required are is negative it means we need to contract our room
            if required_area < 0:
                # Set a function to supervise if each expansion step is succesful or not
                def contract_step () -> bool:
                    # Get the most suitable frontier to contract and try to contract it
                    # If the contraction fails, try with the next one
                    # DANI: De momento uso la misma lógica que la de la expansión porque no va mal
                    # DANI: i.e. evitar contraer fronteras del padre y priorizar fronteras libres es bueno
                    for frontier, loan_permission in self.get_best_frontiers(restricted_segments, contraction=True):
                        if self.contract_frontier(frontier, -required_area):
                            return True
                    return False
                # Contract
                if contract_step():
                    # Refresh how much area we need to expand and go for the next step
                    required_area = self.get_required_area()
                    continue
            # In case the fitting failed restore the backup and exit
            self.parent.restore_children_backup(backup)
            return False
        return True

    # Yield all room frontiers in the most suitable order:
    # - Free frontiers before borther frontiers
    # - Single frontiers before combined frontiers
    # Expansion (default):
    # - In case of a borther frontier, the one which makes shorter the path to free space
    # Contraction:
    # - Parent frontiers are also suitable, but they are the last try
    def get_best_frontiers (self,
        restricted_segments : List[Segment] = [],
        contraction : bool = False,
        ) -> Generator[ Tuple[ Segment, bool ], None, None ]:

        # Get the exterior polygon of the room boundary
        exterior_polygon = self.boundary.exterior_polygon
        # Get the inside corners
        inside_corners = exterior_polygon.get_inside_corners()

        # ----------------------------------------------------------------------------------------------------
        # Split the room boundary exterior polygon in segments according to what each region is connected to
        # Regions connected to the parent free space are desired to expand
        # Regions connected to other rooms may be expanded if the colliding room is able to expand also
        # Regions connected to the parent limits will never be expanded
        # ----------------------------------------------------------------------------------------------------

        free_frontiers, brother_frontiers, parent_frontiers = self.get_frontiers()

        # Filter out frontiers which are inside restricted segments
        def is_restricted (segment : Segment) -> bool:
            for restricted_segment in restricted_segments:
                if segment.get_overlap_segment(restricted_segment):
                    return True
            return False

        free_frontiers = [ frontier for frontier in free_frontiers if not is_restricted(frontier) ]
        brother_frontiers = [ frontier for frontier in brother_frontiers if not is_restricted(frontier) ]
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
            return frontier.length >= self.preventive_min_size or frontier.a in inside_corners or frontier.b in inside_corners

        # Short frontiers could be also expanded together with other connected and aligned frontiers
        # These compound frontiers must be taken in count also although they are harder to expand
        # Each frontier may form 0, 1 or 2 combined frontiers
        def get_combined_frontiers (frontier : Segment) -> Generator[Segment, None, None]:
            for point in frontier.points:
                other_frontiers = [ seg for seg in expandable_frontiers if seg != frontier ]
                combined_frontier = frontier
                implicated_rooms = frontier.rooms
                next_point = point
                # Join connected frontiers until we have enough length or the next frontier is not aligned
                while combined_frontier.length < self.preventive_min_size:
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
                    next_point = next(p for p in connected_frontier.points if p != next_point)
                # If the combined frontier reaches the minimum length add it to the list to be returned
                if combined_frontier.length >= self.preventive_min_size:
                    # Get unique implicated rooms
                    implicated_rooms = unique(implicated_rooms)
                    # Create a fit combined frontier which takes the minimum possible part from the last added frontier
                    # This way the expansion takes as much possible from the main frontier room
                    extended_direction = (point + next_point).normalized()
                    oposite_point = next(p for p in frontier.points if p != point)
                    fit_combined_frontier = Segment(oposite_point, oposite_point + extended_direction * self.preventive_min_size)
                    # Add all combined frontier rooms to the fit combined frontier
                    fit_combined_frontier.rooms = implicated_rooms
                    yield fit_combined_frontier
                    # Add all combined frontier rooms to the combined frontier
                    combined_frontier.rooms = implicated_rooms
                    # Now yield the whole combined frontier
                    yield combined_frontier

        # Check each frontier to be available for expansion alone
        # In case it is, yield it
        # In case it is not, skip it and let it for the end
        def priorize_single_frontiers (frontiers_group : list) -> Generator[Segment, None, None]:
            hard_frontiers = []
            for frontier in frontiers_group:
                # If it cannot be expanded alone we skip it by now
                if not is_suitable(frontier):
                    hard_frontiers.append(frontier)
                    continue
                # Otherwise try to expand it
                yield frontier
            # If none of the suitable frontiers was expanded successfully then try with the combined ones
            for frontier in hard_frontiers:
                combined_frontiers = get_combined_frontiers(frontier)
                for combined_frontier in get_combined_frontiers(frontier):
                    yield combined_frontier

        # Ask for all segments first without any loaned push
        # If all of them fail then retry with the loaned push allowed
        def first_normal_then_loaned (frontiers : Generator[Segment, None, None]) -> Generator[ Tuple[ Segment, bool ], None, None ]:
            already_tried_frontiers = []
            for frontier in frontiers:
                yield frontier, False
                already_tried_frontiers.append(frontier)
            # Pushed loans are not allowed when contracting
            if contraction:
                return
            for frontier in already_tried_frontiers:
                yield frontier, True

        # Find for each brother room the number of colliding rooms we must jump to find free space
        # Then use this value to set the "score" of each brother room frontiers and sort them
        def sort_by_shortest_path (frontiers_group : list) -> list:
            colliding_rooms = unique([ frontier.rooms[0] for frontier in frontiers_group ])
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
                        free_frontiers, brother_frontiers, parent_frontiers = current_room.get_frontiers()
                        if len(free_frontiers) > 0:
                            searching_free = False
                            break
                        current_frontiers += brother_frontiers
                    if not searching_free:
                        #print('ROOM ' + colliding_room.name + ' -> SCORE ' + str(counter))
                        break
                    previous_rooms = previous_rooms + current_rooms
                    current_rooms = unique([ frontier.rooms[0] for frontier in current_frontiers if frontier.rooms[0] not in previous_rooms ])
                    if len(current_rooms) == 0:
                        counter = None
                        #print('ROOM ' + colliding_room.name + ' -> DEAD END')
                        break
                    counter += 1
                # Set the scores for all frontiers
                colliding_room_frontiers = [ frontier for frontier in frontiers_group if frontier.rooms[0] == colliding_room ]
                if counter:
                    for frontier in colliding_room_frontiers:
                        frontier.score = counter
                else:
                    for frontier in colliding_room_frontiers:
                        frontier.score = inf
            # Now sort frontiers using previous scores
            def by_score (frontier):
                return frontier.score
            return sorted(frontiers_group, key=by_score)

        # Start trying to expand free frontiers first
        # First of all and always we shuffle them at random
        random.shuffle(free_frontiers)
        # Sort free frontiers priorizing single frontiers since they are easier to expand
        # Try all free frontiers with loaned push not allowed
        # If all of them fail then try with loaned push allowed
        sorted_free_frontiers = first_normal_then_loaned(priorize_single_frontiers(free_frontiers))
        # Yield frontiers and loan permissions
        for frontier, loan_permission in sorted_free_frontiers:
            yield frontier, loan_permission
        # Then we try with the brother frontiers
        random.shuffle(brother_frontiers)
        # In case we are expanding, sort them according to the shortest route to the free space
        if not contraction:
            brother_frontiers = sort_by_shortest_path(brother_frontiers)
        # Sort brother frontiers priorizing single frontiers since they are easier to expand
        # Try all brother frontiers with loaned push not allowed
        # If all of them fail then try with loaned push allowed
        sorted_brother_frontiers = first_normal_then_loaned(priorize_single_frontiers(brother_frontiers))
        # Yield frontiers and loan permissions
        for frontier, loan_permission in sorted_brother_frontiers:
            yield frontier, loan_permission
        # In case we are contracting, yield also the parent frontiers
        # Note that there is no priority sort for this situation
        if contraction:
            for frontier in parent_frontiers:
                yield frontier, False
        print('WARNING: There are no more frontiers available')

    # Try to expand a specific room frontier
    # Note that the frontier must contain the room it belongs to
    # Set the required (maximum) area it can expand
    # Return True if the expansion was succesful or False if there was no expansion
    def expand_frontier (self, frontier : Segment, required_area : number, allowed_loan_push : bool = False) -> bool:
        # Set the push segment protocol according to the loan permission
        push_protocol = 3 if allowed_loan_push else 1
        # Get the parent room
        parent_room = self.parent
        # Get the exterior polygon of the room boundary
        exterior_polygon = self.boundary.exterior_polygon
        # Get the inside corners
        inside_corners = [ corner for corner in exterior_polygon.corners if corner.inside == True ]
        # Set some parameters according to the rooms this frontier belongs to
        rooms = frontier.rooms

        # Get the grid of the rooms we are about to invade as a reference to calculate how much we can expand
        # When the invaded room is self room it means we are expanding over free space
        # In case this frontier has more than 1 rooms it means it is a combined frontier
        # Create a 'ficticious' room with all implicated rooms as a reference
        grid = None
        for next_room in rooms:
            next_grid = next_room.free_grid if next_room == parent_room else next_room.grid
            # If we have not grid at this point it must mean we are expanding over a parent without grid
            # In this case we must generate a fake grid which allows the room as far as it may need but without conflicting with its brothers
            if not next_grid:
                if next_room != parent_room:
                    raise ValueError('The parent room is the only room which may lack a grid at this point')
                next_grid = parent_room.generate_exterior_free_grid(margin_size=self.max_size*2)
            # Add the current grid to the previous accumulated grid. if any
            if not grid:
                grid = next_grid
                continue
            grid = grid.get_merge_grid(next_grid)
        # Find the maximum rectangles which are in contact with our frontier
        max_rects = grid.max_rects
        contact_max_rects = [ rect for rect in max_rects if frontier in rect ]
        # In case segment is vertical:
        # - The expansion uses maximum columns
        # - The expansion 'forward' is the x dimension
        # - The expansion 'sides' is the y dimension
        # In case segment is horizontal
        # - The expansion uses maximum rows
        # - The expansion 'forward' is the y dimension
        # - The expansion 'sides' is the x dimension
        # Here 'forward' and 'sides' may mean '0 -> x dimension' or '1 -> y dimension'
        # This is because rectangles size is given in a (x,y) tuple format
        # So size[0] = x and size[1] = 1
        if frontier.is_vertical():
            rects = [ column[0] for column in grid.columns ]
            forward = 0
            sides = 1
        elif frontier.is_horizontal():
            rects = [ row[0] for row in grid.rows ]
            forward = 1
            sides = 0
        else:
            raise ValueError('Diagonal segments are not supported for room expansion')

        if len(rects) == 0:
            raise RuntimeError('Empty grid: ' + str(grid))

        # One and only one of the rows/columns will always include the segment
        space = next((rect for rect in rects if frontier in rect), None)
        if space == None:
            # If this happens it may mean there is a problem with the grid
            raise RuntimeError('Frontier ' + str(frontier) + ' has no space in ' + str(rects))
        space_contact = next(segment for segment in space.segments if frontier in segment)
        space_forward_limit = space.get_size()[forward]

        # Find the direction of the expansion as a vector
        forward_direction = space.get_side_direction_to_center(space_contact)

        # Get the forward expansion limit according to maximum rectangles
        maximum_forward_limit = max([ rect.get_size()[forward] for rect in contact_max_rects ])

        # Set the margin according to the invaded rooms
        room_limits = [ self.parent_free_limit if room == parent_room else room.min_size for room in rooms ]
        margin_limit = max(room_limits)

        # Use the margin limit to set the margined forward limits
        margined_space_forward_limit = space_forward_limit - margin_limit
        margined_maximum_forward_limit = maximum_forward_limit - margin_limit

        # Now create a function to make a rectangle by pushing a segment (which may change)
        # The forward length of this rectangle will depend on the available area and limits
        # This segment may not be the original frontier, but a variation of it
        # The push may happen using 3 different protocols
        # - Protocol 1 (greedy): It tries to expand as much as it can, respecting the required area
        # - Protocol 2 (moderate): It tries to expand the minimum possible according to maximum rects
        # - Protocol 3 (loaned): It tries to expand as much as it can even taking more area than needed
        def push_segment (pushed_segment : Segment, protocol : int = 1) -> bool:
            push_length = required_area / pushed_segment.length
            # If the push length at this point is 0 or close to it then we can not push
            # Try to reduce the segment
            if push_length < minimum_resolution:
                print('WARNING: The push length is too small for segment ' + str(pushed_segment))
                new_point_a = pushed_segment.a
                new_point_b = pushed_segment.a + pushed_segment.direction * self.preventive_min_size
                pushed_segment = Segment(new_point_a, new_point_b)
                print('WARNING: segment has been reduced to ' + str(pushed_segment))
                push_length = required_area / self.preventive_min_size
            # In case this is an insider segment which is not wide enought to be pushed alone,
            # Find out how much we can push this segment
            # i.e. find the connected frontier/s and get the maximum length of these segments
            corner_push_limit = None
            if pushed_segment.length < self.preventive_min_size - minimum_resolution:
                other_segments = [ segment for segment in exterior_polygon.segments if pushed_segment not in segment]
                for point in pushed_segment.points:
                    if point in inside_corners:
                        insider = next(segment for segment in other_segments if point in segment)
                        if not corner_push_limit or corner_push_limit < insider.length:
                            corner_push_limit = insider.length
                # If the push length exceeds the insider limit then stay in the limit
                if push_length > corner_push_limit:
                    push_length = corner_push_limit
            # First, try to expand using the maximum available space
            # This is risky since it may split the invaded room in two parts or make regions which do not respect the minimum size
            if protocol == 1:
                # If the length exceeds the maximum limit then stay in the maximum limit
                if push_length > maximum_forward_limit:
                    push_length = maximum_forward_limit
                # If the push length is equal to the maximum forward limit
                elif equal(push_length, maximum_forward_limit):
                    pass
                # If the length is between the maximum limit and the margined maximum limit then stay at the margin
                elif push_length > margined_maximum_forward_limit:
                    push_length = margined_maximum_forward_limit
            # If the greedy try fails use the moderate try
            # Expand only using the closest row/column rectangle
            # This expansion is smaller but safe
            elif protocol == 2:
                # # If at this point the length exceeds the space limit then stay in the space limit
                if push_length > space_forward_limit:
                    push_length = space_forward_limit
                # If the length is between the space limit and the margined space limit then stay at the margin
                elif push_length > margined_space_forward_limit:
                    push_length = margined_space_forward_limit
                # In case the length is bigger than the margined maximum we stay at the margin
                if push_length > margined_maximum_forward_limit:
                    push_length = margined_maximum_forward_limit
            # If the moderate try fails too then use the loaned try
            # Expand as much as possible not tanking in count the required area
            # This expansion may fix a situation where a corner is not claimed because the required area is not enought
            # Note that this room will need to contract further in order to finally get the required area
            elif protocol == 3:
                limits = [ maximum_forward_limit ]
                if corner_push_limit != None:
                    limits.append(corner_push_limit)
                push_length = min(limits)
            # If the push length at this point is 0 or close to it then we can not push
            # This is because of forward limits, not because the area was not big enought
            if push_length < minimum_resolution:
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
            # In case we made a loaned push now we may have more area than we need (usually)
            # We need to return the extra area now
            # Otherwise, if we invade other rooms it may happen that there is not free space enought for them to expand after
            # WARNING: Note that this is an exceptional situation and there are a few extra things to take in count
            # WARNING: At this point there is interior polygons overlap in the parent boundary
            # WARNING: At this point the room is not coordianted with other rooms or the parent free grid
            #          The recently pushed segments must remain as they are
            #          They do not exist for other rooms so me must exlcude them during the fitting to avoid inconsistency
            if protocol == 3 and self.get_required_area() < 0:
                #print('LOANED PUSH -> NEW REQUIRED AREA: ' + str(self.get_required_area()))
                if not self.fit_to_required_area(restricted_segments=new_segments):
                    self.boundary = backup_boundary
                    for i, modified_room in enumerate(rooms):
                        modified_room.boundary = backup_room_boundaries[i]
                    return False
            for room in rooms:
                # If we claimed parent (free) space there is no need to check anything
                if room == parent_room:
                    continue
                # Now invade other rooms
                current_room_invaded_regions = room.boundary.get_overlap_boundaries(invaded_region)
                if not room.invade(current_room_invaded_regions):
                    self.boundary = backup_boundary
                    for i, modified_room in enumerate(rooms):
                        modified_room.boundary = backup_room_boundaries[i]
                    # If the invasion fails then try it again with the next protocol
                    # If we are in the last protocol then return False
                    if protocol == 1:
                        return push_segment(pushed_segment, 2)
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
            return push_segment(frontier, push_protocol)

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
            if new_point.get_distance_to(other_point) < self.preventive_min_size:
                print('WARNING: The reduced frontier is not wide enought')
                return False
            reduced_frontier = Segment(new_point, other_point)

        return push_segment(reduced_frontier, push_protocol)

    # Try to contract a specific room frontier
    # Note that the frontier must contain the room it belongs to
    # Set the required (maximum) area it can contract
    # Return True if the contraction was succesful or False if there was no contraction
    def contract_frontier (self, frontier : Segment, required_area : number) -> bool:
        # Get the parent room
        parent_room = self.parent
        # Get the exterior polygon of the room boundary
        exterior_polygon = self.boundary.exterior_polygon
        # Get the inside corners
        inside_corners = [ corner for corner in exterior_polygon.corners if corner.inside == True ]
        # Find the maximum rectangles which are in contact with our frontier
        grid = self.grid
        max_rects = grid.max_rects
        contact_max_rects = [ rect for rect in max_rects if frontier in rect ]
        # In case segment is vertical:
        # - The expansion uses maximum columns
        # - The expansion 'forward' is the x dimension
        # - The expansion 'sides' is the y dimension
        # In case segment is horizontal
        # - The expansion uses maximum rows
        # - The expansion 'forward' is the y dimension
        # - The expansion 'sides' is the x dimension
        # Here 'forward' and 'sides' may mean '0 -> x dimension' or '1 -> y dimension'
        # This is because rectangles size is given in a (x,y) tuple format
        # So size[0] = x and size[1] = 1
        if frontier.is_vertical():
            rects = [ column[0] for column in grid.columns ]
            forward = 0
            sides = 1
        elif frontier.is_horizontal():
            rects = [ row[0] for row in grid.rows ]
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
        margin_limit = self.preventive_min_size
        margined_space_forward_limit = space_forward_limit - margin_limit

        # Find the direction of the expansion as a vector
        forward_direction = space.get_side_direction_to_center(space_contact)

        # Get the forward expansion limit according to maximum rectangles
        maximum_forward_limit = max([ rect.get_size()[forward] for rect in contact_max_rects ])
        margined_maximum_forward_limit = maximum_forward_limit - margin_limit

        # Now create a function to make a rectangle by pulling a segment (which may change)
        # The forward length of this rectangle will depend on the available area and limits
        # This segment may not be the original frontier, but a variation of it
        # The pull may happen using 2 different protocols, which are tried in the following order
        # - Protocol 1 (greedy): It tries to pull as much as it can, respecting the required area
        # - Protocol 2 (moderate): It tries to pull the minimum possible according to maximum rects
        def pull_segment (pushed_segment : Segment, protocol : int = 1) -> bool:
            push_length = required_area / pushed_segment.length
            # First, try to pull using the maximum available space
            # This is risky since it may split the invaded room in two parts or make regions which do not respect the minimum size
            if protocol == 1:
                # If the length exceeds the maximum limit then stay in the maximum limit
                if push_length > maximum_forward_limit:
                    push_length = maximum_forward_limit
                # If the length is between the maximum limit and the margined maximum limit then stay at the margin
                elif push_length > margined_maximum_forward_limit:
                    push_length = margined_maximum_forward_limit
            # If the greedy try fails use the moderate try
            # Pull only using the closest row/column rectangle
            # This expansion is smaller but safe
            elif protocol == 2:
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
            # WARNING: This usually happens because of forward limits, not because the area was not big enought
            # For this reason, trying to reduce the segment and push again will have no effect almost always
            if push_length < minimum_resolution:
                print('WARNING: The pull length is too small: ' + str(push_length))
                return False
            # Create the new rect with the definitive length
            new_point = pushed_segment.a + forward_direction.normalized() * push_length
            new_side = Segment(pushed_segment.a, new_point)
            new_rect = Rect.from_segments([pushed_segment, new_side])
            removed_boundary = Boundary(Polygon.from_rect(new_rect))
            # We must substract the new rect from this room and check everything is fine after
            if not self.truncate([removed_boundary]):
                # If the contraction fails then try it again with the next protocol
                # If we are in the last protocol then return False
                if protocol != 2:
                    return pull_segment(pushed_segment, 2)
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
            return pull_segment(frontier)

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
            if new_point.get_distance_to(other_point) < self.preventive_min_size:
                print('WARNING: The reduced frontier is not wide enought')
                return False
            reduced_frontier = Segment(new_point, other_point)

        return pull_segment(reduced_frontier)

    # Remove part of the room space and check everything is fine after
    # Use the force argument to skip all checkings and simply truncate the boundary
    # Use the skip_update_display to avoid this truncate to generate a display frame
    # * This is useful when truncating several rooms at the same time, so we do not have to recalculate the parent free grid every time
    def truncate (self, regions : List[Boundary], force : bool = False, skip_update_display : bool = False) -> bool:
        # Calculate the boundary of this room after substracting the invaded regions
        truncated_grid = self.grid
        for region in regions:
            truncated_grid = truncated_grid.get_substract_grid(region.grid)
        # Check the truncated grid to still respecting the minimum size
        if not force and not truncated_grid.check_minimum(self.preventive_min_size):
            print('WARNING: The room is not respecting the minimum size -> Go back')
            return False
        truncated_boundaries = truncated_grid.find_boundaries()
        # In case the room has been splitted in 2 parts as a result of the invasion we go back
        if not force and len(truncated_boundaries) > 1:
            print('WARNING: The room has been splitted -> Go back')
            return False
        # DANI: Es posible que se coma toda la habitación??
        if not force and len(truncated_boundaries) == 0:
            print('WARNING: The invaded room has been fully consumed -> Go back')
            return False
        # Modify the room boundary
        self.set_boundary(truncated_boundaries[0], skip_update_display=skip_update_display)
        return True

    # Invade this room by substracting part of its space
    # Then this room must expand to recover the lost area
    def invade (self, regions : List[Boundary], force : bool = False, skip_update_display : bool = False) -> bool:
        # Truncate the room boundary but save a backup in case we have to go back further
        backup_boundary = self.boundary
        if not self.truncate(regions, force, skip_update_display):
            print('WARNING: The invaded region can not be truncated -> Go back')
            self.boundary = backup_boundary
            return False
        # Expand the invaded room as much as the invaded area
        # In case the expansions fails go back
        if not self.fit_to_required_area():
            print('WARNING: The invaded room can not expand -> Go back')
            self.boundary = backup_boundary
            return False
        return True

    # Invade children in this room by substracting part of their space
    # Then all children must expand to recover the lost area
    # The invasor is the children room which will claim this new space
    def invade_children (self, invasor : 'Room', region : Boundary) -> bool:
        # Make a backup of all children, which includes the invasor
        # All children boundaries will be recovered if only 1 child fails to get invaded
        backup = self.make_children_backup()
        # Claim the invaded region for the invasor room
        invasor.merge_boundary(region)
        # Now find which children are overlaped by the invade region and then truncate them
        children = [ child for child in self.children if child != invasor ]
        for child in children:
            # Get the overlap regions between the invaded region and this child
            overlap_boundaries = child.boundary.get_overlap_boundaries(region)
            if len(overlap_boundaries) == 0:
                continue
            # In case something went wrong with any child truncation recover all children backups and stop
            if not child.truncate(overlap_boundaries):
                self.restore_children_backup(backup)
                return False
        # Finally expand all truncated children
        for child in children:
            # In case something went wrong with any child expansion recover all children backups and stop
            if not child.fit_to_required_area():
                self.restore_children_backup(backup)
                return False
        return True

    # Fuse a boundary to self room boundary
    # Check that boundary can be joined to current boundary as a single room
    # i.e. both boundaries must be colliding and the colliding region must be only one and as wide as the minimum size or more
    def merge_boundary (self, boundary : Boundary):
        if self.boundary == None:
            self.boundary = boundary
        else:
            self.boundary = self.boundary.merge_boundary(boundary, self.preventive_min_size)

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
    # i.e. "free of conflict" frontiers, "conflict" frontiers and "forbidden" frontiers respectively
    def get_frontiers (self) -> tuple:
        parent_room = self.parent
        # In case we have a parent boundary, frontiers with this boundary are totally forbidden to expansion
        if parent_room.boundary:
            # The parent limits are not allowed for expansion
            parent_frontiers = self.get_frontiers_with(parent_room)
            # Other rooms inside the same parent may be displaced if there is no free space available
            brother_rooms = [ room for room in parent_room.children if room is not self and room.boundary ]
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
            for segment in self.boundary.exterior_polygon.segments:
                free_frontiers += segment.substract_segments(conflict_frontiers)
            # Set the room all free segment belong to as None
            for frontier in free_frontiers:
                frontier.rooms = [parent_room]
            return free_frontiers, brother_frontiers, parent_frontiers
        # If there is no paren boundary then any frontier which has no brothers conflict is suitable for expansion
        else:
            # Other rooms inside the same parent may be displaced if there is no free space available
            brother_rooms = [ room for room in parent_room.children if room is not self and room.boundary ]
            brother_frontiers = []
            for room in brother_rooms:
                # This function assign the frontier.rooms already
                frontiers = self.get_frontiers_with(room)
                if frontiers:
                    brother_frontiers += frontiers
            # Now susbstract all previous frontiers from the rest of the external polygon to find free frontiers
            free_frontiers = []
            for segment in self.boundary.exterior_polygon.segments:
                free_frontiers += segment.substract_segments(brother_frontiers)
            # Set the room all free segment belong to as None
            for frontier in free_frontiers:
                frontier.rooms = [parent_room]
            return free_frontiers, brother_frontiers, []

    # Go uppwards in the hyerarchy until you reach the room which has no parent
    def get_root_room (self) -> 'Room':
        root = self
        while root.parent:
            root = root.parent
        return root

    # Get this room and all children rooms recursively
    def get_rooms_recuersive (self) -> List['Room']:
        rooms = [self]
        for room in self.children:
            rooms += room.get_rooms_recuersive()
        return rooms

    # Make a backup of current children boundaries
    def make_children_backup (self) -> List[Boundary]:
        return [ child.boundary for child in self.children ]

    # Restore children boundaries using a backup
    def restore_children_backup (self, backup):
        for c, child in enumerate(self.children):
            child.boundary = backup[c]

    # Add a new frame in the display with the current segments of this room and its children
    # Also an 'extra' argument may be passed with extra segments to be represented
    def update_display (self, extra : list = []):
        if not display_solving_process:
            return
        # Find the root room
        root = self.get_root_room()
        # Get all children rooms recursively
        rooms = root.get_rooms_recuersive()
        # Display all current rooms together
        elements_to_display = [ *rooms, *extra ]
        add_frame(elements_to_display)

# Exception for when user input is wrong
class InputError(Exception):
    pass

# Solve room distributions
# The display flag may be passed in order to generate a dynamic graph to display the solving process
def solve(room : Room, display : bool = False):
    global display_solving_process
    display_solving_process = display
    room.set_children_boundaries(recursive=True)

# -----------------------------------------------------

# Auxiliar functions

# Given a list of segments (path) return the sum of their lengths
def get_path_length (path : List['Segment']) -> number:
    return sum([ segment.length for segment in path ])