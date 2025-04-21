from typing import List, Tuple, Dict, Union, Optional

from auxiliar import *
from scheme_display import add_frame, display_options
from vectorial_base import *

import random
from math import sqrt, inf, tan, atan, degrees, radians

# Set the seed and print it
seed = None

#seed = 679238 # Escalera separada de la pared y dejando discarded grid
#seed = 898828 # Queda free space una vez resuleto todo, sin fallos
#seed = 562106 # Queda free space una vez resuleto todo, sin fallos
seed = 322825

if not seed:
    seed = round(random.random() * 999999)
print(f'Seed {seed}')
random.seed(seed)

# A door is a segment in a boundary
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
        # Set the length of the required space to each side of the door
        length : Optional[number] = None,
        # Set the direction the door is open thorugh
        direction : Optional[Vector] = None,
        # Set the pivot point for the door to open
        pivot : Optional[Point] = None,
        # Set if the door is rigid
        # i.e. its point may not change as a result of the solving process
        rigid : bool = False,
        # Set if the doors direction points outside the room instead of inside
        reverse : bool = False,
        # Set the parent room
        room : Optional['Room'] = None
    ):  
        # These values are usually None at this point
        # They are usually set further from the door room 'door_args' value
        self._width = width
        self._margin = margin
        self._length = length
        self._point = point
        self._margined_width = None
        self._segment = None
        self._margined_segment = None
        self._direction = direction
        self._pivot = pivot
        self.rigid = rigid
        if self.rigid and not self.point:
            # DANI: Esto es solo de momento, pero podemos hacer que el rigid tenga efecto una vez se asigna el punto
            raise InputError('A point must be defined if the rigid flag is passed')
        self.reverse = reverse
        # The room this door belongs to
        self.room = room

    def __repr__ (self):
        point = str(self.point) if self.point else 'No point'
        width = str(self.width) if self.width else 'No width'
        margin = str(self.margin) if self.margin else 'No margin'
        return f'<Door {point} {width} ({margin})>'

    # Get the width
    def get_width (self) -> number:
        # If we have a stored value already then return it
        if self._width != None:
            return self._width
        # Otherwise we must get it from the parent room args
        # If there is no parent room then we have nothing to do
        if not self.room:
            return None
        self._width = self.room.door_args['width']
        return self._width

    # Set the width (regular setter)
    def set_width (self, new_width : number):
        self._width = new_width
    
    # The door width
    width = property(get_width, set_width, None, "The door width")

    # Get the margin
    def get_margin (self) -> number:
        # If we have a stored value already then return it
        if self._margin != None:
            return self._margin
        # Otherwise we must get it from the parent room args
        # If there is no parent room then we have nothing to do
        if not self.room:
            return None
        self._margin = self.room.door_args['margin']
        return self._margin

    # Set the margin (regular setter)
    def set_margin (self, new_margin : number):
        self._margin = new_margin
    
    # The door margin
    margin = property(get_margin, set_margin, None, "The door margin")

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

    # Get the length
    def get_length (self) -> number:
        # If we have a stored value already then return it
        if self._length != None:
            return self._length
        # Otherwise we must get it from the parent room args
        # If there is no parent room then we have nothing to do
        if not self.room:
            return None
        self._length = self.room.door_args['length']
        return self._length

    # Set the length (regular setter)
    def set_length (self, new_length : number):
        self._length = new_length
    
    # The door length
    length = property(get_length, set_length, None, "The door required space length on each side")

    # Get the door point
    def get_point (self) -> Optional[Point]:
        return self._point

    # If the door point is set then reset its segment and margined segment
    def set_point (self, point : Optional[Point]):
        if self.rigid:
            raise RuntimeError(f'Can not change rigid door point in door {self}')
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
        # If the door has no point then complain
        if not self.point:
            raise RuntimeError('Trying to get door margined segment when no point is defined')
        self._margined_segment = self.generate_segment(self.margined_width)
        return self._margined_segment

    # The door margined segment
    margined_segment = property(get_margined_segment, None, None, "The door margined segment")

    # Given a segment width, generate a new
    # The new segment will be centered in the door point
    # The new segment will be overlaped with the boundary segment where the door point is
    def generate_segment (self, width : number) -> Segment:
        # If width is 0 then the segment can not exist
        if width == 0:
            raise ValueError('Cannot generate a segment for a door of width 0')
        # If the door point is not assigned then we can not generate the segment
        if not self.point:
            return None
        # If we can not retrieve the boundary then we can not generate the segment
        room_boundary = self.get_room_boundary()
        if not room_boundary:
            return None
        # Otheriwse, generate the margined segment
        boundary_segment = next(( segment for segment in room_boundary.segments if self.point in segment ), None)
        if not boundary_segment:
            raise ValueError(f'The door point {self.point} is not over its room boundary ({self.room.name})')
        direction = boundary_segment.direction
        half_width = width / 2
        a = self.point - direction * half_width
        b = self.point + direction * half_width
        if a not in boundary_segment or b not in boundary_segment:
            raise ValueError(f'The door segment ({Segment(a,b)}) does not fit in its room boundary ({self.room.name})')
        return Segment(a,b)

    # Make a backup of the current door
    def make_backup (self) -> dict:
        return {
            'point': self._point,
            'segment': self._segment,
            'margined_segment': self._margined_segment,
            'direction': self._direction,
            'pivot': self._pivot
        }

    # Restore a backup
    def restore_backup (self, backup : dict):
        self._point = backup['point']
        self._segment = backup['segment']
        self._margined_segment = backup['margined_segment']
        self._direction = backup['direction']
        self._pivot = backup['pivot']

    # Generate a rect containing the minimum required space for this door
    # If inside is true (default) then the rect in the inside side of the room is returned
    # Otherwise the outside side rect is returned
    def get_required_space (self, inside : bool = True) -> Optional[Rect]:
        margined_segment = self.margined_segment
        if not margined_segment:
            return None
        if margined_segment not in self.room.boundary:
            return None
        inside_direction = self.get_inside_direction()
        direction = inside_direction if inside else -inside_direction
        second_segment = margined_segment.translate(direction * self.length)
        return Rect.from_segments([margined_segment, second_segment])

    # Get the door inside direction
    # i.e. the direction towards the door boundary inside side
    def get_inside_direction (self) -> Optional[Vector]:
        if not self.segment:
            return None
        room_boundary = self.get_room_boundary()
        if not room_boundary:
            return None
        inside_direction = room_boundary.get_border_inside(self.segment)
        return inside_direction

    # Get the door direction
    # i.e. the direction the door is open thorugh
    # By default the direction points to the door boundary inside side
    def get_direction (self) -> Optional[Vector]:
        # Return internal value if it exists
        if self._direction:
            return self._direction
        # Otherwise, find the direction
        inside_direction = self.get_inside_direction()
        self._direction = -inside_direction if self.reverse else inside_direction
        return self._direction
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
        if not self.segment or not self.room:
            return None
        room_boundary = self.get_room_boundary()
        if not room_boundary:
            return None
        door_points = self.segment.points
        boundary_segment = next( segment for segment in room_boundary.segments if self.point in segment )
        boundary_segment_points = boundary_segment.points
        outside_corners = [ corner for corner in room_boundary.corners if corner in boundary_segment_points and not corner.inside ]
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

    # Get the boundary where the door is meant to be
    def get_room_boundary (self) -> Optional[Boundary]:
        room = self.room
        if not room:
            return None
        boundary = room.boundary
        if not boundary:
            return None
        return boundary

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
        # If not available segments are passed then we use the room boundary segments after substracting other doors
        if available_segments == None:
            # The door must have a boundary
            room_boundary = self.get_room_boundary()
            if not room_boundary:
                raise ValueError('Cannot find a suitable region for a room without boundary')
            # Get segments in its boundary which are wide enought for the door
            candidate_segments = [ segment for segment in room_boundary.segments if segment.length >= minimum_segment_length ]
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
            if lower(segment.length, minimum_segment_length):
                continue
            # Cut the margins of all available segments for the reminaing segments to be available to store the door point (center)
            elif equal(segment.length, minimum_segment_length):
                suitable_point = segment.get_middle_point()
                suitable_points.append(suitable_point)
            else:
                suitable_segment = segment.get_margined_segment(minimum_segment_length / 2)
                suitable_segments.append(suitable_segment)
        return suitable_segments, suitable_points

    # Find suitable unique points and extreme points from suitable segment regions
    def find_suitable_points (self, available_segments : Optional[List[Segment]] = None) -> List['Point']:
        # Find suitable regions
        suitable_segments, suitable_points = self.find_suitable_regions(available_segments)
        if len(suitable_segments) == 0 and len(suitable_points) == 0:
            return []
        # Reduce to suitable points
        suitable_points += sum([ list(segment.points) for segment in suitable_segments ], [])
        return suitable_points

    # Check if a point in a boundary is siutable for placing this door
    # If no point is provided then self point is used by default
    # If no boundary is provided then self room boundary is used by default
    # i.e. check if the point is in the door room boundary and there is space enough around for its margined width
    def is_point_suitable (self, point : Optional[Point] = None, boundary : Optional[Boundary] = None) -> bool:
        # If input point is missing set self point
        if point == None:
            # If there is no self point either then we have nothing to do
            if self.point == None:
                raise ValueError('No point was provided and the door has no point already')
            point = self.point
        # If width is 0 then the segment can not exist
        if self.width == 0:
            raise ValueError('Cannot generate a segment for a door of width 0')
        # If boundary is not assigned we use the room boundary
        if not boundary:
            boundary = self.get_room_boundary()
            # If there is not room boundary either then we have nothing to do
            if not boundary:
                raise ValueError('No boundary was provided and there is no room boundary to place the door')
        # Otheriwse, generate the margined segment
        boundary_segment = next(( segment for segment in boundary.segments if self.point in segment ), None)
        if not boundary_segment:
            return False
        direction = boundary_segment.direction
        half_width = self.width / 2
        a = self.point - direction * half_width
        b = self.point + direction * half_width
        if a not in boundary_segment or b not in boundary_segment:
            return False
        return True

    # Relocate self door in a suitable region in contact with the parent corridor
    # Find suitable regions in the room boundaries which contact the corridor
    # The boundary to place the door is the room boundary by default, but a custom boundary may be passed
    def relocate (self, boundary : Optional[Boundary] = None) -> bool:
        # If the door has not room or width then it makes not sense to relocate the door
        if not self.room or not self.width:
            return False
        # If the door is rigid then it may not be relocated
        if self.rigid:
            return False
        parent_room = self.room.parent
        # If there is no parent (i.e. room is the root) then the parent corridor is not a restriction
        parent_corridor = parent_room.corridor_grid
        if not parent_corridor:
            return False
        parent_corridor_boundaries = parent_corridor.boundaries if parent_room else []
        parent_corridor_segments = sum([ boundary.segments for boundary in parent_corridor_boundaries ],[])
        # Get the available segments to place the door
        # Regions where the exterior boundary of the door room and any parent corridor converges
        available_boundary = boundary if boundary else self.room.boundary
        available_segments = available_boundary.get_segments_overlap_segments(parent_corridor_segments)
        # If there are not available segments at this point then we can not relocate the door
        if len(available_segments) == 0:
            return False
        # Find suitable regions to place the door according to its size and margins among the available segments
        suitable_segments, suitable_points = self.find_suitable_regions(available_segments)
        if len(suitable_segments) == 0 and len(suitable_points) == 0:
            return False
        # Reduce to suitable points
        suitable_points += sum([ list(segment.points) for segment in suitable_segments ], [])
        # Now priorize those suitable regions which are a better placement for the door
        # The ideal place is which makes the open door to stay next to a wall
        def sort_by_distance_to_wall (point : Point) -> number:
            boundary_segment = available_boundary.get_border_element(point)
            corners = [ available_boundary.get_corner(p) for p in boundary_segment.points ]
            outside_corners = [ corner for corner in corners if not corner.inside ]
            if len(outside_corners) == 0:
                return inf
            corner_distances = [ point.get_distance_to(corner) for corner in outside_corners ]
            return min(corner_distances)
        suitable_points.sort(key=sort_by_distance_to_wall)
        most_suitable_point = suitable_points[0]
        # If the current point is the most suitable point already then exit here
        if self.point == most_suitable_point:
            return False
        # Set the most suitable point as the current door point
        # Note that setting the point already resets the segment, the pivot, etc.
        self.point = most_suitable_point
        return True

    # Check if the door is placed in a suitable point already and, if not, try to relocate it
    # Return Ture both if it was fine already or it could be relocated
    def check_and_relocate (self, boundary : Optional[Boundary] = None) -> bool:
        if self.is_point_suitable(boundary=boundary):
            return True
        return self.relocate(boundary=boundary)
        
# A room is a smart boundary that may contain other boundaries with conservative areas and size restrictions
class Room:
    def __init__ (self,
        # A start boundary may be passed. If no boundary is passed it is assigned automatically
        boundary : Optional[Boundary] = None,
        # Set the minimum and maximum area which must be accomplished at the end of the solving process
        # Failure to set the room within this area range will result in process kill
        # Area may be a number (absolute value) or a string (percent) with the format 'XX%'
        # If no area is passed then the original boundary area will be used
        min_area : Optional[Union[number, str]] = None,
        max_area : Optional[Union[number, str]] = None,
        # Minimum size in both x and y dimensions
        min_size : Optional[number] = None,
        # Set if the room boundary may me modified
        # It is true by default when an input boundary is passed
        # It is false by default otherwise
        rigid : Optional[bool] = None,
        # Set the maximum number of corners in the final room boundary
        # Note that this parameter may be passed only when the boundary is not passed
        # Root only
        max_corners : Optional[number] = None,
        # Corridor size (the minimum size by default)
        corridor_size : Optional[number] = None,
        # The doors to enter the room and start corridors
        doors : Optional[List[Door]] = None,
        # The default inputs for all doors whose inputs are not specified
        # These options are inherited by children rooms whose options are not specified
        door_args : Optional[Dict] = None,
        # Set the room height
        # This height is the used to raise the 3D scheme
        # However, it may also affect the 2D distribution, since stairs rely on this height to calculate their lengths
        height : Optional[number] = None,
        # The name and color parameters are only representation parameters and they have no effect in the logic
        name : str = 'Unnamed',
        segments_color : str = 'black',
        fill_color : str = 'white',
        # Set other rooms inside this room
        children : List['Room'] = [],
        # Set the parent room from the inputs
        # Note that this can be done only when the parent has been previously initialized
        parent : Optional['Room'] = None,
        # Set the parent building (if any)
        # The parent building may be used to guess some default values
        # However the parent building is not mandatory even in the root room
        parent_building : Optional['Room'] = None,
    ):
        # WARNING: Note that the declaration of many variables in the init must be in a strict order
        # An internal boundary and save the input boundary
        self._boundary = None
        self.input_boundary = boundary
        # Set if the boundary is rigid
        if rigid == None:
            self.rigid = bool(boundary)
        else:
            if rigid == True and not boundary:
                raise InputError('Cannot set a rigid room and not provide a boundary')
            self.rigid = rigid
        # The forced grid is the space which must be occupied by a room whose boundary is not fixed (child adaptable)
        # Note that this value is totally redundant for rooms with fixed boundaries
        # Note also that the input boundary becomes the forced grid by default
        # Thus it makes total sense passing a boundary and flagging rigid as false
        self.forced_grid = boundary.grid if boundary else Grid()
        self._grid = Grid()
        self._free_grid = None
        self._available_grid = None
        self._rigid_grid = None
        # Set representation parameters
        self.name = name
        self.segments_color = segments_color
        self.fill_color = fill_color
        # Save the height
        self._height = height
        # Set up the hierarchy of rooms
        # Parent is never assigned from the instance itself, but it is assigned by the parent
        self._parent = None
        self.parent = parent
        # Save the parent building
        self.parent_building = parent_building
        # Set an atrribute to sabe the corridor grid when it is set
        self._corridor_grid = Grid()
        # Set the boundary
        # If the boundary has been forced then update the display with the initial segments
        if boundary:
            # We must skip first display since some internal variables are not yet declared and it may fail
            self.set_boundary(boundary, skip_update_display=True)
        # Set if this room boundary is to be stablished by its children
        # This is set True when the root room has not a prestablished boundary
        # This variable may be true only in the root room (if this is not the root then it will become False further)
        # Note that all the logic changes drastically whether the root room has or not a predefined boundary
        self._child_adaptable_boundary = boundary == None and min_area == None and max_area == None
        # Set the children rooms
        self._children = None
        self.children = children
        # Check input minimum and maximum area formats to make sense
        for value in [ min_area, max_area ]:
            # A numeric value is accepted
            # Also a None is accepted, even if there is no boundary
            # Note that this may be the root and have a child adaptable boundary
            if is_number(value) or value == None:
                continue
            # String values are allowed as long as they have a percent format (i.e. 'XX%')
            elif type(value) == str:
                if not value[-1] == '%':
                    raise InputError(f'Area range ({value}) has a non-supported format in room {name}')
            else:
                raise InputError(f'Area range ({value}) has a non-supported format in room {name}')
        # Set the internal values for minimum and maximum areas
        self.input_min_area = min_area
        self.input_max_area = max_area
        # Set the internal value for the expected final area
        self._target_area = None
        # Set the maximum number of corners
        self.max_corners = max_corners
        if max_corners != None:
            if boundary:
                raise InputError('Maximum number of corners is to be passed only when there is no predefined boundary')
            if max_corners < 4:
                raise InputError('It is not supported having less than 4 corners in the final boundary')
        # Set size limits
        self._min_size = min_size
        # Set the internal max size value
        self._max_size = None
        # Parent free limit is set by the parent while setting the child boundary
        # Parent free limit is the maximum min size of all parent children but this child
        self.parent_free_limit = None
        # Save input door args
        self._door_args = door_args
        # Set the doors
        self._doors = None
        self.input_doors = doors
        # If there is no input doors then set a single defualt door, by default
        # Note that to set no doors the input must be '[]' instead of 'None'
        self.doors = doors if doors != None else [ Door() ]
        # Set the room corridor size
        self._corridor_size = corridor_size
        # Set an empty grid as the discarded space
        self._discarded_grid = None
        # Children handling:
        if self.children:
            # Now set some parameters in children
            for child in self.children:
                child._post_init(self)
            # If a boundary was passed then check children to fit on it
            if self.area:
                # Check areas of all children rooms to do not sum up more than the parent area
                # In addition check if any of the children room has boundary and, if so, check the boundary is inside the parent
                children_minimum_area = sum([ child.min_area for child in self.children ])
                if greater(children_minimum_area, self.area):
                    raise InputError('Children together require more area than the parent has')

    # Set a function which sets other initial values once the parent has been set
    # This function is called by the parent once it has been initiated
    def _post_init (self, parent):
        # Set the child adaptable boundary as false, since this is not the root room
        self._child_adaptable_boundary = False
        # In case the parent has a boundary, check this room may fit on it
        if parent.boundary:
            # Check this room is inside the parent boundary, if they have a predefined boundary
            if self.boundary and self.boundary not in parent.grid:
                raise InputError('The child room "' + self.name + '" is out of the parent boundary')
            # Check if this rooms minimum size fits in the parent boundary
            if not self.boundary and not parent.does_room_fit(self, force=True):
                raise InputError('The child room "' + self.name + '" minimum size does not fit in the parent boundary')
        # If this room (not the root) has a max_corners input parameter
        if self.max_corners:
            raise InputError('Parameter max_corners is supported only in the root room')

    def __str__(self):
        return '<Room "' + str(self.name) + '">'

    def __repr__(self):
        return '<Room "' + str(self.name) + '">'

    # Get the parent room (normal getter)
    def get_parent (self) -> 'Room':
        return self._parent
    # Set the parent room
    def set_parent (self, parent : 'Room'):
        self._parent = parent
        # If there is no parent then there is nothing else to do here
        if parent == None:
            return
        # Update parent children to maintain coherence
        parent._children.append(self)
        # Reset the corresponding parent grids
        if self.grid:
            parent.reset_free_grid()
            parent.reset_available_grid()
            if self.rigid:
                parent.reset_rigid_grid()
                parent.reset_discarded_grid()

    # The parent room
    parent = property(get_parent, set_parent, None, "The parent room")

    # Get the children rooms (normal getter)
    def get_children (self) -> List['Room']:
        return self._children

    # Set the children rooms
    # Update hierarchy
    def set_children (self, children : List['Room']):
        # Make sure there are no duplicates
        unique_children = unique(children)
        if len(children) != len(unique_children):
            raise RuntimeError(f'Duplicated child in room {self.name}: ' + ', '.join([ child.name for child in children]))
        # Set every child parent as self to maintain coherence
        for child in children:
            child._parent = self
        self._children = children
        # If new children are added then some of them may have overlapping boundaries
        # e.g. when the stairs rooms are added
        # For this reason we must reset the free grid and some other inner grids
        self.reset_free_grid()
        self.reset_available_grid()
        self.reset_rigid_grid()
        self.reset_discarded_grid()

    # The children rooms
    children = property(get_children, set_children, None, "The children rooms")

    # Get the brother rooms
    def get_brother_rooms (self) -> List['Room']:
        if not self.parent:
            return []
        return [ room for room in self.parent.children if room != self ]

    # Get the boundary
    def get_boundary (self) -> Boundary:
        # Return the internal stored value if already exists
        if self._boundary:
            return self._boundary
        # Otherwise, get the grid boundary
        # WARNING: Note that we are accessing the public grid value, not the private
        # WARNING: The grid does never read the public boundary value because we would fall in a loop
        # The grid leads over the boundary in adaptable boundary scenarios
        if self.grid:
            boundaries = self.grid.boundaries
            if len(boundaries) > 1:
                raise ValueError('A single room can not have more than 1 boundary')
            self._boundary = boundaries[0]
            return self._boundary
        return None

    # Set the boundary
    # Reset the own grid, free grid and parent free grid
    def set_boundary (self, boundary : Boundary, skip_update_display : bool = False):
        self.set_grid(boundary.grid, skip_update_display)

    # The room boundary
    boundary = property(get_boundary, set_boundary, None, "The room boundary")

    # Get the grid
    # Just return the internal boundary grid value
    def get_grid (self) -> Optional[Grid]:
        # Note that we explicitly evaluate the grid not to be None
        # This is because the free grid may be an empty grid (e.g. child adaptable boundary when there are not children yet)
        # An empty grid evaluates as false, but it has a meaning and thus it must be returned
        if self._grid != None:
            return self._grid
        if self._boundary:
            return self._boundary.grid
        if self._child_adaptable_boundary:
            self._grid = self.get_provisional_grid()
            return self._grid
        return None

    # Set the grid
    # Reset the boundary
    # DANI: Esto, aunque debería funcionar, nunca se ha provado y no se usa actualmente
    def set_grid (self, grid : Grid, skip_update_display : bool = False):
        # Show areas for a better debug
        previous_area = self._grid.area if self._grid else 0
        new_area = grid.area if grid else 0
        area_difference = round_to_hundredths(new_area - previous_area)
        # Set the new grid
        self._grid = grid
        # Reset the boundary
        self._boundary = None
        # Reset the corresponding inner grids
        self.reset_free_grid()
        self.reset_available_grid()
        self.reset_discarded_grid()
        # Update the parent
        parent_room = self.parent
        if parent_room:
            # Reset the corresponding inner grids
            parent_room.reset_free_grid()
            # In case the parent boundary was provisional
            if parent_room._child_adaptable_boundary:
                parent_room.reset_grid()
            # If this is a rigid room then we must check if we have produced discarded regions in the parent
            if self.rigid:
                # Update parent rigid and discarded grids
                parent_room.reset_available_grid()
                parent_room.reset_rigid_grid()
                parent_room.reset_discarded_grid()
        # Update current display if the flag to skip the update is not passed
        if not skip_update_display:
            self.update_display(title=f'Modified grid in room {self.name} (area difference {area_difference})')

    # The room boundary
    grid = property(get_grid, set_grid, None, "The room grid")

    # Reset room grid
    # This must be done, for instance, when the room is child adaptable and some child grid is modified
    def reset_grid (self):
        self._boundary = None
        self._grid = None

    # Get the minimum area
    def get_min_area (self) -> number:
        # Return the value as it is, in case it is numeric
        if is_number(self.input_min_area):
            return self.input_min_area
        # If area is a string then it means it is a percent and this it must be parsed
        if type(self.input_min_area) == str:
            if not self.parent:
                raise RuntimeError(f'Can not solve a percent area since parent room is not defined in room {self.name}')
            portion = float(self.input_min_area[0:-1]) / 100
            return portion * self.parent.area
        # If not then calculate the minimum area
        # If no input minimum area was passed then a boundary must have been passed
        if self.input_min_area == None:
            if self.area and self.rigid:
                return self.area
            raise InputError(f'Minimum area can not be guessed in room {self.name}. Please specify an area range or provide a rigid boundary')
        raise ValueError(f'Not processable minimum area {self.input_min_area} in room {self.name}')

    # The minimum area to be reached at the end of the solving process (read only)
    min_area = property(get_min_area, None, None, "Minimum area to be reached after the solving process")

    # Get the maximum area
    def get_max_area (self) -> number:
        # Return the value as it is, in case it is numeric
        if is_number(self.input_max_area):
            return self.input_max_area
        # If area is a string then it means it is a percent and this it must be parsed
        if type(self.input_max_area) == str:
            if not self.parent:
                raise RuntimeError(f'Can not solve a percent area since parent room is not defined in room {self.name}')
            portion = float(self.input_max_area[0:-1]) / 100
            return portion * self.parent.area
        # If not then calculate the maximum area
        # If no input maximum area was passed then a boundary must have been passed
        if self.input_max_area == None:
            if self.area and self.rigid:
                return self.area
            raise InputError(f'Maximum area can not be guessed in room {self.name}. Please specify an area range or provide a rigid boundary')
        raise ValueError(f'Not processable maximum area {self.input_max_area} in room {self.name}')

    # The maximum area to be reached at the end of the solving process (read only)
    max_area = property(get_max_area, None, None, "Maximum area to be reached after the solving process")

    # Get the target area
    def get_target_area (self) -> number:
        # Return the internal value, if any
        if self._target_area != None:
            return self._target_area
        # If the room is rigid then the target area is the already occupied area
        if self.rigid:
            return self.area
        # We must calculate target areas
        if not self.parent:
            raise RuntimeError('Can not calculate target area since parent room is not defined in room '  + self.name)
        self.parent.set_children_target_areas()
        # Before returning the target area make a fast checking
        # If target area does not cover the minimum size then it makes no sense
        if self.min_size and self.target_area and self.target_area < self.min_size**2:
            print(f'Target area {self.target_area} is not reachable with a minimum size of {self.min_size}')
            raise InputError('Room area is not sufficient for the minimum size in room ' + self.name)
        return self._target_area

    # Room target area (read only)
    target_area = property(get_target_area, None, None, "Target area to be reached after the solving process")

    # Find out the donable area
    # i.e. the area over the minimum area
    def get_donable_area (self) -> number:
        if self.area < self.min_area:
            return 0
        return self.area - self.min_area
    donable_area = property(get_donable_area, None, None, "Donable area over the minimum area")

    # Find out the donable area proportion
    # i.e. the area over the minimum area in proportion to the area range between the minimum and the maximum areas
    def get_donable_area_proportion (self) -> number:
        area_range = self.max_area - self.min_area
        return self.donable_area / area_range
    donable_area_proportion = property(get_donable_area_proportion, None, None, "Donable area over the minimum area in proportion to area range between the minimum and the maximum areas")

    # Set a function to set all children target areas at once
    # Note that all children areas have to be assigned together
    def set_children_target_areas (self, verbose : bool = False):
         # We assign target areas for non rigid rooms only
        non_rigid_children = [ child for child in self.children if not child.rigid ]
        if verbose: print('Setting children target areas of ' + ', '.join([ child.name for child in non_rigid_children ]))
        # First of all check all children area ranges to make sense
        # Note that this is calculated here and not in the init for a reason:
        # One area could be a percent and the other a number, for instance
        for child in non_rigid_children:
            if lower(child.max_area, child.min_area):
                raise InputError('Maximum area must be higher than minimum area in room ' + child.name)
        # If the parent (self) is child adaptable then children target areas are set randomly
        if self._child_adaptable_boundary:
            for child in non_rigid_children:
                child._target_area = random.uniform(child.min_area, child.max_area)
                if verbose: print(f'  {child.name} target area: {child._target_area} (random)')
            return
        # Get the available area to fit the non-rigid children
        available_area = self.available_area
        if verbose: print(f'  Available area: {available_area}')
        # Get all minimum areas and calculate the overall minimum area needed
        min_areas = [ child.min_area for child in non_rigid_children ]
        overall_min_area = sum(min_areas)
        # Check we have enough area to allocate all children
        if greater(overall_min_area, available_area):
            raise RuntimeError(f'Parent {self.name} has not enough available area {available_area} to allocate all its children minimum areas {overall_min_area}')
        # Get all maximum areas and calculate the overall maximum area which may be occupied
        max_areas = [ child.max_area for child in non_rigid_children ]
        overall_max_area = sum(max_areas)
        # If the parent is able to allocate all children in their maximum areas then we set all target areas as the maximum
        if lower(overall_max_area, available_area) or equal(overall_max_area, available_area):
            for child in non_rigid_children:
                child._target_area = child.max_area
                if verbose: print(f'  {child.name} target area: {child._target_area} (maximum)')
            return
        # If the parent is not able to allocate all children in their maximum areas then target areas must be calculated
        overall_range = overall_max_area - overall_min_area
        overall_percent = (available_area - overall_min_area) / overall_range
        for child in non_rigid_children:
            child_range = child.max_area - child.min_area
            child._target_area = child.min_area + child_range * overall_percent
            if verbose: print(f'  {child.name} target area: {child._target_area} (proportion)')

    # Area inside the room boundary (read only)
    def get_area(self) -> number:
        boundary = self.boundary
        if not boundary:
            return 0
        return boundary.area
    area = property(get_area, None, None, "Area inside the room boundary")

    # Find different configurations of inner grids
    # This is the region not occupied by inner elements such as:
    # All children rooms, rigid children rooms only, the corridor or discarded regions
    def find_inner_grid (self, remove_non_rigid_children : bool = False) -> Grid:
        # Return none if there is not boundary yet
        if not self.boundary:
            return None
        # Start including all the room regions
        inner_grid = self.grid
        # Then remove non free/available regions
        # Remove all children regions if we are strict
        # Remove only rigid children regions otherwise
        for child in self.children:
            # Make sure the child has a grid
            if not child.grid:
                continue
            if remove_non_rigid_children or child.rigid:
                inner_grid -= child.grid
        # Remove the corridor grid
        if self.corridor_grid:
            inner_grid -= self.corridor_grid
        # Remove the discarded grid
        if self.discarded_grid:
            inner_grid -= self.discarded_grid
        return inner_grid

    # Get the free space inside the room boundary as a rectangles grid
    # i.e. space not filled by children rooms or the corridor
    def get_free_grid (self) -> Optional[Grid]:
        # If free grid is previously calculated then return it
        # Note that we explicitly evaluate the free grid not to be None
        # This is because the free grid may be also an empty grid (i.e. there is no free space)
        # An empty grid evaluates as false, but it has a meaning and thus free grid must not be calculated again
        if self._free_grid != None:
            return self._free_grid
        # Calculate the free grid otheriwse
        # Then save it and return it
        self._free_grid = self.find_inner_grid(remove_non_rigid_children=True)
        return self._free_grid

    # Set the free grid as None for it to be recalculated again next time it is needed
    def reset_free_grid (self):
        self._free_grid = None

    # Free space grid (read only)
    free_grid = property(get_free_grid, None, reset_free_grid, "The room free space grid")

    # Free space area (read only)
    def get_free_area (self) -> number:
        return self.free_grid.area
    free_area = property(get_free_area, None, None, "Free space area (read only)")

    # Get the available grid
    def get_available_grid (self) -> Grid:
        # If available grid is previously calculated then return it
        if self._available_grid != None:
            return self._available_grid
        # Calculate the available grid otheriwse
        # Then save it and return it
        self._available_grid = self.find_inner_grid()
        return self._available_grid

    # Set the available grid as None for it to be recalculated again next time it is needed
    def reset_available_grid (self):
        # If the grid is alredy None then there is nothing to do here
        if self._available_grid == None:
            return
        # Set the grid as none so the next time it is getted it will be recalculated
        self._available_grid = None
        # Reset children target areas as well since they will need to be recalculated
        for child in self.children:
            child._target_area = None

    # Available grid is the grid not occupied by rigid children, the corridor grid or the discarded grid
    available_grid = property(get_available_grid, None, reset_available_grid, "The room available space grid, excluding discarded regions")

    # Available space area (read only)
    def get_available_area (self) -> number:
        return self.available_grid.area
    available_area = property(get_available_area, None, None, "Available space area (read only)")

    # Get the corridor size
    def get_corridor_size (self) -> number:
        # Return the internal value, if any
        if self._corridor_size != None:
            return self._corridor_size
        # Otherwise, find valid corridor size
        # Try to inherit the parent corridor size
        if self.parent:
            self._corridor_size = self.parent.corridor_size
            return self._corridor_size
        # If we are the root but there is a parent building then inherit its value
        if self.parent_building:
            self._corridor_size = self.parent_building.room_args.get('corridor_size', None)
            if self._corridor_size != None:
                return self._corridor_size
        # Otherwise, rely in the minimum size
        self._corridor_size = self.min_size
        return self._corridor_size

    # Set the corridor size (regular setter)
    def set_corridor_size (self, new_corridor_size : number):
        self._corridor_size = new_corridor_size

    # The corridor size
    corridor_size = property(get_corridor_size, set_corridor_size, None, "The corridor size")

    # Get the corridor grid
    def get_corridor_grid (self) -> Optional[Grid]:
        return self._corridor_grid
    def set_corridor_grid (self, new_corridor_grid : Optional[Grid], skip_update_display : bool = False):
        # Save the new corridor grid
        self._corridor_grid = new_corridor_grid
        # Reset the corresponding inner grids
        self.reset_free_grid()
        self.reset_available_grid()
        # Update rigid and discarded grids
        self.reset_rigid_grid()
        self.reset_discarded_grid()
        # If the room is child adaptable then reset its grid
        # Corridor may be outside its boundaries
        if self._child_adaptable_boundary:
            self.reset_grid()
        # Update current display if the flag to skip the update is not passed
        if not skip_update_display:
            self.update_display(title=f'Modified corridor grid in room {self.name}')

    # Free space grid (read only)
    corridor_grid = property(get_corridor_grid, set_corridor_grid, None, "The room corridor space grid")

    # Get the rigid grid
    def get_rigid_grid (self) -> Grid:
        # Return the internal value, if any
        if self._rigid_grid != None:
            return self._rigid_grid
        # Calculate the value otherwise
        self._rigid_grid = Grid()
        # Add rigid children grids
        for child in self.children:
            if child.rigid and child.grid:
                self._rigid_grid += child.grid
        # Add the corridor grid
        if self.corridor_grid:
            self._rigid_grid += self.corridor_grid
        return self._rigid_grid
    # Reset the rigid grid
    def reset_rigid_grid (self):
        self._rigid_grid = None
    # The rigid grid is the one made of rigid room grids and the corridor grid
    rigid_grid = property(get_rigid_grid, None, reset_rigid_grid, "The rigid space grid")

    # Simple getter
    def get_discarded_grid (self) -> Grid:
        # Return the internal value, if any
        if self._discarded_grid != None:
            return self._discarded_grid
        # Otherwise we must calculate it
        self._discarded_grid = Grid()
        # Find regions between the whole grid and the rigid grid respecting the minimum size
        non_rigid_grid = self.grid - self.rigid_grid
        safe_grid = non_rigid_grid.keep_minimum(self.min_size)
        # The grid which is removed from the truncated grid by not respecting the minimum must be flaged as discarded
        lost_grid = non_rigid_grid - safe_grid
        if lost_grid:
            self._discarded_grid += lost_grid
        # Reset the available grid
        self.reset_available_grid()
        # Return the new calculated discarded grid
        return self._discarded_grid
    # Reset the discarded grid
    def reset_discarded_grid (self):
        self._discarded_grid = None
    # Set a grid for discarded space
    # This is space which is not suitable to be claimed and thus it is not set free to avoid the solver to try it
    # This space may be generated when setting the corridor and it may be impossible to recover
    discarded_grid = property(get_discarded_grid, None, reset_discarded_grid, "The discarded regions which are not possible to be claimed")

    # Get the minimum size (read only)
    def get_min_size (self) -> number:
        # If we have a stored value already then return it
        if self._min_size:
            return self._min_size
        # Get the parent minimum size
        if self.parent:
            self._min_size = self.parent.min_size
            return self._min_size
        # If we are the root but we are part of a building then get the minimum size from the room args
        if self.parent_building:
            self._min_size = self.parent_building.room_args.get('min_size', None)
            return self._min_size
        raise
        raise InputError(f'Cannot find minimum size for room {self.name}. Please, set the minimum size in the root room or the buildings room args')
    min_size = property(get_min_size, None, None, "The room minimum space")

    # Get the maximum size (read only)
    def get_max_size (self) -> number:
        # If we have a stored value already then return it
        if self._max_size:
            return self._max_size
        # Calculate the maximum size: the maximum possible size in case the boundary was the thinest rectangle
        # Then save it and return it
        if not self.min_size:
            return inf
        self._max_size = self.target_area / self.min_size
        return self._max_size
    max_size = property(get_max_size, None, None, "The room maximum space")

    # Get the door arguments
    def get_door_args (self) -> dict:
        # If we have a stored value already then return it
        if self._door_args:
            return self._door_args
        # If we are not the root, inherit parent default door arguments
        if self.parent:
            self._door_args = self.parent.door_args
            return self._door_args
        # If we are the root but there is a parent building then inherit its value
        if self.parent_building:
            self._door_args = self.parent_building.room_args['door_args']
            return self._door_args
        # DANI: Este trozo de código está repetido y por defecto es el otro el que correr
        # DANI: Habría que replantearselo
        # If we are the root and there is no parent building, we have to guess the most suitable door arguments
        # Make the margined width of all doors equal to the minimum room minimum size
        # Make the width of all doors the 80% of the margined width
        margined_width = self.corridor_size
        # If we have not a margined width at this point me must complain about the inputs
        if not margined_width:
            raise InputError('Cannot guess the door width. Please set a minimum size or provide explicit "door_args" in the root room or in the parent building "room_args"')
        width = margined_width * 0.8
        margin = margined_width * 0.1
        length = margined_width
        self._door_args = {
            'width': width,
            'margin': margin,
            'length': length
        }
        return self._door_args

    # Set the door arguments (regular setter)
    def set_door_args (self, new_door_args : dict):
        self._door_args = new_door_args

    # Arguments to set doors by default in this room
    door_args = property(get_door_args, set_door_args, None, "Arguments to set doors by default in this room")

    # Get the doors (regular getter)
    def get_doors (self) -> List[Door]:
        return self._doors

    # Set the doors
    # As soon as doors are set, set self room as the parent room of each door
    def set_doors (self, new_doors : List[Door]):
        self._doors = new_doors
        for door in self._doors:
            door.room = self
            # If any door is rigid then we can already assign a corridor region to this door
            if door.rigid:
                self.parent.corridor_grid += Grid([door.get_required_space(inside=False)])

    # The room doors
    doors = property(get_doors, set_doors, None, "The room doors")

    # Get the height
    def get_height (self) -> number:
        # If there is a stored value already then return it
        if self._height != None:
            return self._height
        # Otherwise, the parent height as the current room height
        # Note that we are reading the public height, not the private
        # Thus it is a recursive setter until we find a parent with a height value
        if self.parent:
            return self.height
        # If we are the root and there is a parent building we can inherit the height from it
        if self.parent_building:
            return self.parent_building.room_args['height']
        # If there is no parent (i.e. this room is the root) there is noheight to guess, so we complain
        raise InputError('Cannot guess ' + self.name + ' height. Please set a height at least in the root room/building')
        
    # Set the height (regular setter)
    def set_height (self, new_height : number):
        self._height = new_height

    # The room height
    height = property(get_height, set_height, None, "The room height")

    # Set children grids according to the room configuration
    # This function triggers the logic to solve room distributions
    # If the recursive flag is passed then solve each child's children grids and so on recursively
    # All children rooms must have their boundary fully set before solving the next generation of children
    def solve_children (self, recursive : bool = False, verbose : bool = False) -> bool:
        rooms = self.children
        # If there are not children then we have nothing to do here
        if len(rooms) == 0:
            return True
        # Before allocating place for any children room, set corridor spaces in front of parent doors
        # This way this space is never claimed since it has to be occupied by the corridor anyway
        # Note that the root room may have non defined doors
        already_set_doors = [ door for door in self.doors if door.point ]
        for door in already_set_doors:
            self.corridor_grid += Grid([door.get_required_space(inside=True)])
        # Sort children rooms by minimum size, with the biggest sizes first
        if self._child_adaptable_boundary:
            random.shuffle(rooms)
        else:
            def sort_by_size (room):
                return room.min_size
            rooms.sort( key=sort_by_size, reverse=True )
        # Set up each room by giving them a position and correct size to match the forced area
        for room in rooms:
            # Configure the child room to respect the parent free min size limit according to its brothers
            other_rooms = [ self ] + [ other for other in rooms if other != room ]
            parent_free_limit = max([ other.min_size for other in other_rooms ])
            room.parent_free_limit = parent_free_limit
            # If the room already has boundary at this point then ignore it
            if room.boundary: continue
            # If the children has no boundary it must be built
            # Use the proper function depending on if this is a child adaptable boundary parent
            child_grid_setter = self.set_child_grid_adaptable if self._child_adaptable_boundary else self.set_child_grid
            if not child_grid_setter(room):
                if verbose: print(f'Child {room.name} has failed to be set')
                return False

        # Reshaped children to reduce unnecessary corners as well
        self.reduce_children_corners()

        # Now that all children bondaries are set we must set the corridor

        # Set the corridor
        if len(rooms) > 0:
            if not self.set_corridor():
                return False

        # At this point the boundary is no longer adaptable to child boundaries, in case it was
        # This is because then the reducing corneres process requires real free space to work
        self._child_adaptable_boundary = False

        # If there is a limit of corners in the room (i.e. this is the root room) then reshape children now
        if self._child_adaptable_boundary and self.max_corners:
            self.reduce_corners()
        # Reshape children to reduce unnecessary corners as well
        self.reduce_children_corners()

        # Relocate the doors to the most suitable placement now that boundaries will change no more
        any_relocated_door = False
        for child in self.children:
            for door in child.doors:
                if door.rigid:
                    continue
                if door.relocate():
                    any_relocated_door = True

        # Show the relocated doors
        if any_relocated_door:
            self.update_display(title='Relocated doors')

        # Solve children recurisvely if the recursive flag was passed
        if recursive:
            for room in rooms:
                if not room.solve_children(recursive=True):
                    return False

        # We successfully solved it
        return True

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
    # Search in free space by default and all space if the argument 'forced' is passed
    def get_room_fitting_rects (self, room : 'Room', forced : bool = False) -> List[Rect]:
        size = room.min_size
        # Get all fitting rects
        grid = self.available_grid if forced else self.free_grid
        fitting_rects = list(grid.generate_fitting_rects(size, size))
        return fitting_rects

    # Set up a child room boundary with no parent boundary
    def set_child_grid_adaptable (self, room : 'Room') -> bool:
        # There are no boundary restrictions
        # All children will be set in one single step, as the maximum initial grids
        # To do so, we must find a suitable corner and space (rect) for the boundary to be set
        # This space will be as long as the double of the room maximum size, for it to be able to expand without problems
        huge_size = room.max_size * 2
        # First of all, find all children with already set boundaries/grid
        # They are the only reference to place the new room
        brother_boundaries = [ child.boundary for child in self.children if child.boundary ]
        # If this is the first room then there is no reference at all, so we set the first corner as the (0,0) point
        # Then the "available space" will be a fake rect which is huge enough for this room to set its whole boundary
        if len(brother_boundaries) == 0:
            space = Rect(0, 0, huge_size, huge_size)
            # This is the "real" corner, which include its segments
            initial_point = Point(0,0)
            corner = space.get_corner(initial_point)
            room.grid = room.build_maximum_initial_grid(corner, space)
            return True
        # If there are other childs with boundaries already then we must find the best location for the new room
        # It has to be next to its brother rooms, in a point which minimizes the corridor length further
        # If there is one brother only then it makes not sense to find a corridor
        if len(brother_boundaries) == 1:
            space = Rect(-huge_size, 0, 0, huge_size)
            reference_point = Point(0,0)
            corner = space.get_corner(reference_point)
            room.grid = room.build_maximum_initial_grid(corner, space)
            return True
        # Otherwise, calculate the "pre-corridor" and find the closest point in the current parent "pre-exterior-polygon"
        # First, get the current corridor points
        corridor_backbone = self.find_corridor_backbone()
        corridor_segments = corridor_backbone['segments'] if corridor_backbone else None
        corridor_nodes = corridor_backbone['nodes'] if corridor_backbone else None
        corridor_points = list(set(sum([ segment.points for segment in corridor_segments ], ()))) + corridor_nodes
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
        # If we find a point in the polygon which is already in the corridor then we are done
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
            room.grid = room.build_maximum_initial_grid(corner, space)
            return True
        # In case the point is in a segment we have only one direction to expand and we can choose the other
        polygon_segment = exterior_polygon.get_border_element(minimum_distance_point)
        direction_a = -exterior_polygon.get_border_inside(polygon_segment)
        direction_b = direction_a.rotate(90)
        space_segment_a = Segment(minimum_distance_point + direction_a * huge_size, minimum_distance_point)
        space_segment_b = Segment(minimum_distance_point, minimum_distance_point + direction_b * huge_size)
        corner = Corner(minimum_distance_point.x, minimum_distance_point.y, space_segment_a, space_segment_b)
        space = Rect.from_corner(corner)
        room.grid = room.build_maximum_initial_grid(corner, space)
        return True

    # Set up a child room boundary between parent and brother bundaries
    def set_child_grid (self, room : 'Room', forced : bool = False, verbose : bool = False) -> bool:
        if verbose: print(f'Setting {self.name} child room {room.name} boundary')
        # Find a suitable maximum free rectangle to deploy a starting base grid
        # The minimum base grid is a square with both sides as long as the room minimum size
        # If there are no suitable rects it means this child fits nowhere in the free space
        # This may happen with the last childs, when previous childs have take almost all space
        # In this case we take as availbale space all the room space and then we invade overlapped children
        if forced:
            suitable_rects = self.get_room_fitting_rects(room, forced=True)
            if len(suitable_rects) == 0:
                # DANI: Esto no debería pasar nunca. Debería preveerse de antes
                raise RuntimeError(f'The room {room.name} fits nowhere')
        else:
            suitable_rects = self.get_room_fitting_rects(room)
            if len(suitable_rects) == 0:
                if verbose: print(f'  Room {room.name} fits nowhere in the free space')
                return self.set_child_grid(room, forced=True)
        # Shuffle the suitable rects
        random.shuffle(suitable_rects)
        # Sort the suitable rects by minimum size
        def sort_by_size(rect):
            return min(rect.get_size())
        sorted_suitable_rects = sorted( suitable_rects, key=sort_by_size )
        # Make a backup in case we have to force since other rooms will be modified
        if forced:
            backup = self.make_grid_backup(children_only=True)
        # Try to set up the new room in all possible sites until we find one
        # Each 'site' means each corner in each suitable rectangle
        # Check each site to allow other rooms to grow
        # Pick only corners wich are already in the free grid (no matter if interior or exterior)
        sites = [ (corner, rect) for rect in sorted_suitable_rects for corner in rect.get_corners() ]
        # Save previous grid to not repeat already checked (and failed) grids
        previous_initial_grid = None
        # Iterate the different corners
        for corner, rect in sites:
            # In case we force the base grid we must check which children were invaded
            # In addition, the base grid must be as small as possible
            initial_grid = None
            if forced:
                initial_grid = room.build_minimum_initial_grid(corner)
            # Otherwise we set freely the maximum possible grid
            else:
                initial_grid = room.build_maximum_initial_grid(corner, rect)
            # There must be always an initial grid at this point
            if not initial_grid:
                raise RuntimeError('Something went wrong with initial grid')
            # A few tests to avoid inconvinient but not fatal scenarios
            # They are not mandatory and thus they will be skipped if this is forced
            if not forced:
                # If we can, should check parent free grid is not getting more splitted than it is
                current_splits = len(self.free_grid.boundaries)
                new_hypothetical_free_grid = self.free_grid - initial_grid
                new_splits = len(new_hypothetical_free_grid.boundaries)
                if new_splits > current_splits:
                    if verbose: print(f'  Free space is getting more splitted than it was')
                    continue
                # If we can, should check parent doors are not gettin occupied
                if not self.check_doors_side(initial_grid, inside=True):
                    if verbose: print(f'  Not respecting some {self.name} door space')
                    continue
            # At this point, check if the grid is equal to a previous tried (and failed) grid
            # This may happend with the 4 corners of the same rect
            if previous_initial_grid == initial_grid:
                if verbose: print(f'  Grid already tried and failed')
                continue
            # Save the current grid to skip it in further iterations in case it fails
            previous_initial_grid = initial_grid
            # Now actually set the inital grid
            if not room.expand_grid(initial_grid, check_parent_free_grid=False):
                if verbose: print(f'  Failed to set {room.name} initial grid')
                continue
            # Proceed with the expansion of this child room until it reaches its forced area
            if not room.fit_to_required_area():
                if verbose: print(f'  Failed to fit {room.name} required area')
                # If the expansion failed then clean the grid (i.e. the initial grid)
                room.grid = None
                if forced:
                    self.restore_grid_backup(backup)
                continue
            return True
        # If at this point we still have no boundary it means there is no available place to set the perimeter
        # If this was not forced then retry forcing
        if not forced:
            return self.set_child_grid(room, forced=True)
        return False

    # Build a provisional grid from the children room grids
    # This function is meant to be used only when the parent (self) room has an adaptable boundary
    def get_provisional_grid (self) -> Grid:
        # Get children room grids
        children_grids = [ child.grid for child in self.children if child.grid ]
        # Include also the forced grid, if any
        # Note that this is the canonical way to add 'free space' in a child adaptable boundary room
        # You must never get the free grid from here, since the free grid is calculated from the grid and this function calculates the grid
        if self.forced_grid:
           children_grids.append(self.forced_grid)
        # Include also the corridor grid, if any
        if self.corridor_grid:
           children_grids.append(self.corridor_grid)
        # If there are no grids at this point then we are done
        if len(children_grids) == 0:
            return Grid()
        # Merge all grids
        provisional_grid = merge_grids(children_grids)
        return provisional_grid

    # Build a provisional exterior grid
    # This function is meant to be used only when the parent room has not boundary
    # Create a fake square which contains all the children and then add a margin as long as the distance argument
    # This is used by children to expand outwards
    def generate_exterior_free_grid (self, margin_size : number) -> 'Grid':
        exterior_polygon = self.boundary.exterior_polygon
        exterior_box = exterior_polygon.get_box()
        expanded_box = exterior_box.expand_margins(margin_size)
        exterior_free_boundary = Boundary(Polygon.from_rect(expanded_box), [exterior_polygon])
        return exterior_free_boundary.grid

    # Set the initial room grid as the maximum possible rectangle
    # It is useful to set a whole room at the begining, when there is plenty of free space
    def build_maximum_initial_grid (self, corner : Point, space : Rect) -> Grid:
        x_space, y_space = space.get_size()
        # If the room area is greater than the space then return the whole space as a permeter
        if space.area <= self.target_area:
            # If both space sizes are bigger than the maximum size we can not return the whole space
            # The size in one of both sides must be limited (the biggest size is filled)
            if x_space > self.max_size and y_space > self.max_size:
                if x_space >= y_space:
                    maximum_rect = Rect.from_corner(corner, x_space, self.max_size)
                else:
                    maximum_rect = Rect.from_corner(corner, self.max_size, y_space)
                return Grid([maximum_rect])
            return Grid([space])
        # Else, fit the room in the space
        # Try to create the most "squared" possible rectangle
        # Calculate how long would be the side of a perfect square
        square_side_length = sqrt(self.target_area)
        # Set how long will be the short (restricted) side of the rectangle
        # Then calculate the other side length
        # The limit may come from the square side limit, the space limit or the own room limit
        minimum_space = min(x_space, y_space)
        first_side_length = min(square_side_length, minimum_space, self.max_size)
        # In case the first side length is between the limit space and the margined limit space we must fit it to the limit
        margined_minimum_space = minimum_space - self.parent_free_limit
        if minimum_space > first_side_length and first_side_length > margined_minimum_space:
            first_side_length = margined_minimum_space
        # Once we have calculated the first length we can calculate the second one
        # The second length will be over the maximum space
        maximum_space = max(x_space, y_space)
        maximum_side_length = self.target_area / first_side_length
        second_side_length = min(maximum_space, maximum_side_length)
        # In case the second side length is between the limit space and the margined limit space we must fit it to the limit as well
        margined_maximum_space = maximum_space - self.parent_free_limit
        if maximum_space > second_side_length and second_side_length > margined_maximum_space:
            second_side_length = margined_maximum_space
        # Create the new rect fitting the biggest size in the biggest space and the opposite
        # For each size of the new rect use the calculated size only if it is not longer than the space
        if x_space >= y_space:
            new_x_size = second_side_length
            new_y_size = first_side_length
        else:
            new_x_size = first_side_length
            new_y_size = second_side_length
        maximum_rect = Rect.from_corner(corner, new_x_size, new_y_size)
        return Grid([maximum_rect])

    # Set the initial room grid as the minimum possible rectangle
    # This is used when the initial grid must be forced over other children rooms
    def build_minimum_initial_grid (self, corner : Point) -> Grid:
        minimum_rect = Rect.from_corner(corner, self.min_size, self.min_size)
        return Grid([minimum_rect])

    # Using both parent and children boundaries, calculate which segments make the shortest path to connect parent doors and all children rooms
    # Missing children doors are also set during this process and they are placed according to make the corridor as short as possible
    # Note that this function may be used several times to find a "pre-corridor" before definitely setting it
    # NEVER FORGET: The algorithm to solve the corridor finds all possible combinations of nodes to know the final result is the shortest possible path
    # NEVER FORGET: Once I tried to add as many nodes as 'candidate doors' and as a result the algorithm was taking hours to solve the problem
    def find_corridor_backbone (self, verbose : bool = False) -> Optional[dict]:
        if verbose: print(f'Setting {self.name} corridor')

        # Set the corridor nodes
        # This is a preprocessing step which is required to find the shortest corridor and other similar calculations

        # Get the exterior polygon
        exterior_polygon = self.boundary.exterior_polygon
        
        # Use both self boundary and all children room boundary segments
        available_segments = []
        # Set the doors which must be reached by the corridor, both from self and children
        # Doors may have a point (i.e. they are already set) or not
        # Doors which have a point must be reached by the corridor
        # Doors which do not have a point must be set after the corridor and then the corridor may be expanded to reach them
        # Save already stablished door points in a set for nodes to further know if they are doors or not
        already_set_door_points = set()
        # Save rooms which have door nodes since there is no need to include them in the required rooms list
        # Door nodes are mandatory so these rooms will be included anyway
        already_doored_rooms = set()
        # Split self and children boundaries at the doors
        for room in [ self, *self.children ]:
            if not room.boundary:
                continue
            # Get points from already set doors
            current_door_points = [ door.point for door in room.doors if door.point ]
            # Split each of the room segments at the door points
            for room_segment in room.boundary.segments:
                available_segments += room_segment.split_at_points(current_door_points)
            # Update the set with all already set door points
            already_set_door_points.update(current_door_points)
            # Add the current room to the set of already doored rooms
            if len(current_door_points) > 0:
                already_doored_rooms.add(room)
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
        # Set the rooms which are rigid for the path solving
        # i.e. if a path is surrounded by rigid rooms then it is discarded
        # Note that if parent is not child adaptable then its boundaries are also rigid
        path_rigid_rooms = [ room for room in self.children if room.rigid ]
        if not self._child_adaptable_boundary:
            path_rigid_rooms.append(self)
        # Find which rooms are in contact to each node and if nodes are in the exterior boundary
        # Find also which nodes are doors
        for node_point, node_data in nodes.items():
            # First find the node rooms
            rooms = set()
            for child in self.children:
                if not child.boundary:
                    continue
                if node_point in child.boundary.exterior_polygon:
                    rooms.add(child)
            if node_point in exterior_polygon:
                rooms.add(self)
            node_data['rooms'] = list(rooms)
            # Now find out if it is a door
            is_door = node_point in already_set_door_points
            node_data['is_door'] = is_door
            # Note that a door node should always have 2 rooms
            # The exception is an scenario where the corridor is set while there is parent free space available yet
            # Note that if a door node is surrounded by 2 rigid rooms then it will be not reachable by the corridor
            if is_door and len(rooms) == 2 and all([ room in path_rigid_rooms for room in rooms ]):
                raise ValueError(f'Door node {node_point} is not reachable by the corridor. Is it between 2 rigid rooms?')
        # Now find "non-redundant" nodes and the "paths" between them
        # Redundant nodes are those whose contact rooms are already included in all connected nodes
        # Knwoing this is useful when we are expanding our corridor since a reundant node will never solve the puzzle
        # When we expand thorugh redundant nodes we can claim several nodes until we find a non-redundant node
        # Note that nodes with 2 connected segments will always be redundant while others will be non-redundant
        # IMPORTANT: Doors are also non-redundant nodes
        for node_data in nodes.values():
            node_data['is_redundant'] = len(node_data['connected_segments']) == 2 and not node_data['is_door']
        # Save non-redundant nodes apart
        non_redundant_nodes = { k:v for k, v in nodes.items() if v['is_redundant'] == False }
        # If there is no redundant nodes then it makes not sense trying to set a corridor at all
        # DANI: No se si esto es posible, pero just in case
        if len(non_redundant_nodes) == 0:
            return None
        # If there is only one non redundant node then it makes not sense trying to set a corridor. See figure 03
        # However this node is to be set as a door for implicated rooms
        # Fix this very specific situtation by making the only door to be the door of all rooms
        if len(non_redundant_nodes) == 1:
            rooms = [ self, *self.children ]
            doors = sum([ room.doors for room in rooms ], [])
            doors = [ door for door in doors if door.point ]
            if len(doors) != 1:
                raise RuntimeError('Vete tu a saber. Si no se parece a la figura 3 hay que replantearse esto')
            only_door = doors[0]
            for room in rooms:
                if len(room.doors) == 0:
                    continue
                if len(room.doors) > 1:
                    raise RuntimeError('No podemos dar soporte a este escenario. Igual hay que restringirlo todo a una puerta por habitación')
                room.doors[0] = only_door
            # DANI: Esto de aquí abajo en lugar de asignar la misma puerta hace que la puerta ya existente tenga el mismo punto
            # only_node = list(non_redundant_nodes.keys())[0]
            # for room in self.children:
            #     main_door = room.doors[0]
            #     if not main_door.point:
            #         main_door.point = only_node
            return None
        # Then a path between non-redundant nodes is a list of segments, which are connected by redundant nodes
        for node_point, node_data in non_redundant_nodes.items():
            # Find all current node paths to other non-redundant nodes
            paths = []
            # Save the non-redundant node each path is leading to
            path_nodes = []
            for starting_segment in node_data['connected_segments']:
                # Get the path rooms
                # A path must always have 2 and only 2 rooms in a scenario where children have fully consumed parent area
                # However, if the corridor is set while there is still free space it may happen that a node has only 1 room
                path_rooms = [ room for room in node_data['rooms'] if starting_segment in room.boundary.exterior_polygon ]
                # Check if both rooms from this path are rigid rooms
                # In that case we discard the path rigth now since we can not build a corridor here
                if len(path_rooms) == 2 and all([ room in path_rigid_rooms for room in path_rooms ]):
                    continue
                last_segment = starting_segment
                last_point = next(point for point in last_segment.points if point != node_point)
                last_node = nodes[last_point]
                path = [ last_segment ]
                # Set if the room is isolated
                is_isolated = False
                # Keep moving forward until we find a non-redundant node
                while last_node['is_redundant']:
                    # If the node is redundant then there will be always only 2 connected segments and one of them will be the last
                    last_segment = next( segment for segment in last_node['connected_segments'] if segment != last_segment )
                    last_point = next(point for point in last_segment.points if point != last_point)
                    last_node = nodes[last_point]
                    path.append(last_segment)
                    # If we closed a path not finding any new redundant point then this may be an isolated room with a set door
                    # If the room is isolated then it means it has no contact with any brother room neither the parent
                    # The already set door may be in contact with free space or a preliminar corridor grid
                    # In any case there is no need to expand the corridor to reach it and thus we can exclude all this path from the corridor
                    if last_point == node_point:
                        is_isolated = True
                        break
                        #add_frame(path)
                        #raise ValueError('This should never happen. Do we have a room separated from parent and brothers with an already set door?')
                # If the last path belongs to an isolated region then do not include it in the corridor available paths
                if is_isolated:
                    continue
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
        required_children_rooms = [ child for child in self.children if child.boundary and len(child.doors) > 0 and child not in already_doored_rooms ]
        required_rooms = ([ self ] if len(self.doors) > 0 and self not in already_doored_rooms else []) + required_children_rooms
        # Set the variables to store the current corridor
        current_corridor = []
        current_corridor_nodes = []
        current_corridor_length = None
        # At this point we must consider both the already stablished corridor grid and the remaining free grid as corridor contributors
        # However we want to keep this grids independent and not just make the whole free grid the corridor grid
        # Free grid is still space to be claimed and it may be required further while the corridor grid is not to be claimed
        functional_corridor_grid = self.corridor_grid + self.free_grid
        # Get the nodes and segments which are in the free space
        # Note that if free space is splitted then each region is treated independently
        # WARNING: Treating all free regions as one is a problem since it may result in a splitted corridor
        # If a free space region is reached during the solving then include its nodes and segments in the corridor
        # They will be removed further during the corridor expansion, but they must be considered as corridors for colliding rooms to be included
        corridor_free_regions = []
        free_region_nodes = {}
        if functional_corridor_grid:
            # Define each region of the free space separately
            for i, boundary in enumerate(functional_corridor_grid.boundaries):
                # Now find which of the splitted segments are in the free region boundary
                # Note that all the previous segments should be covered by the splitted segments, but we need these segments splitted
                free_corridor_segments = []
                for segment in splitted_segments:
                    for free_space_segment in boundary.segments:
                        if segment in free_space_segment:
                            free_corridor_segments.append(segment)
                # Get those non-redundant nodes which are in the free space
                free_corridor_nodes = []
                for point, node in non_redundant_nodes.items():
                    for free_space_segment in boundary.segments:
                        if point in free_space_segment:
                            # If we already associate the point to one free space segment then go to the next
                            free_corridor_nodes.append(point)
                            break
                # Set a region object with the already set segments and nodes
                corridor_free_region = {
                    'corridor_segments': free_corridor_segments,
                    'corridor_nodes': free_corridor_nodes
                }
                corridor_free_regions.append(corridor_free_region)
                # Now associate this free region data to each of its nodes so it is easier to find later
                # Also create a fake 'node room' to be associated to this free regions and make this has a required room
                # This way we make sure the solver reaches all free regions
                # Note that we can assign a string as the node room since this value is just used as a key
                free_region_room_hash = 'free_region_room_hash_' + str(i)
                required_rooms.append(free_region_room_hash)
                for point in free_corridor_nodes:
                    # Set current free region to this node point
                    free_region_nodes[point] = corridor_free_region
                    # Add the free region room hash to the node data rooms list
                    nodes[point]['rooms'].append(free_region_room_hash)
            # Set the first of the free regions as part of the current corridor already
            sample_free_region = next(corridor_free_region for corridor_free_region in free_region_nodes.values())
            # Update the current corridor values
            # WARNING: The current_corridor_length is not set since it must be only updated when the corridor is complete
            current_corridor = sample_free_region['corridor_segments']
            current_corridor_nodes = sample_free_region['corridor_nodes']

        # Set a function to check if the corridor is finished, given a list of rooms and nodes
        def is_corridor_finished (corridor_rooms : List['Room'], corridor_nodes :List[Point]) -> bool:
            # Check all required rooms are in the corridor
            contains_all_rooms = all(room in corridor_rooms for room in required_rooms)
            if not contains_all_rooms:
                return False
            # Check all required doors are in the corridor
            contains_all_doors = all(door_point in corridor_nodes for door_point in already_set_door_points)
            if not contains_all_doors:
                return False
            return True

        # We check just in case we already have covered all nodes with the current free space corridor
        # In this case there is no need to build a corridor at all, we can stop here
        current_rooms = set(sum([ nodes[point]['rooms'] for point in current_corridor_nodes ],[]))
        if is_corridor_finished(current_rooms, current_corridor_nodes):
            if verbose: print('  There is no need to build a corridor')
            # Set the free grid as the corridor grid
            # DANI: Esto tal vez se puede replantear. No es necesario y perdemos la free grid (espacio claimable)
            self.corridor_grid = functional_corridor_grid
            return None

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
                # Get the remaining available paths/nodes after substracting the current next path
                following_available_paths = available_paths[0:i] + available_paths[i+1:]
                following_available_path_nodes = available_path_nodes[0:i] + available_path_nodes[i+1:]
                # The following node will be the other next path's node
                following_node = available_path_nodes[i]
                # Get the new following path after adding the last path while getting the next node
                following_path = current_path + next_path
                following_path_nodes = current_path_nodes + [ following_node ]
                # Get the corresponding node data
                following_node_data = nodes[following_node]
                # Get the following node paths which are not already included in the current path and its nodes
                following_node_paths = [ *following_node_data['paths'] ]
                following_node_path_nodes = [ *following_node_data['path_nodes'] ]
                # In case we find a free region node we immediately add all its segments and nodes to the corridor
                corridor_free_region = free_region_nodes.get(following_node, None)
                if corridor_free_region:
                    following_path += corridor_free_region['corridor_segments']
                    following_path_nodes += corridor_free_region['corridor_nodes']
                    for node in corridor_free_region['corridor_nodes']:
                        node_data = nodes[node]
                        following_node_paths += node_data['paths']
                        following_node_path_nodes += node_data['path_nodes']
                # Add the following node paths/nodes to the remaning available paths/nodes
                # Then we get the available paths/nodes for the following path
                #print(len(following_node_paths))
                for path, node in zip(following_node_paths, following_node_path_nodes):
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
                # Get the following path covered rooms
                following_rooms = current_rooms.union(set(following_node_data['rooms']))
                # DANI: Usa esto para ver los pasos intermedios
                #elements_to_display = [ segment.get_colored_segment('red') for segment in following_path ]
                #self.update_display(extra=elements_to_display, title='Corridor solver step')
                # If following path includes all rooms then it is a candidate to be the corridor
                if is_corridor_finished(following_rooms, following_path_nodes):
                    # Check if this path is shorter than the current corridor
                    # The shorter path will remain as the current corridor
                    # Also the current corridor length may be none if this is the first attempt
                    following_path_length = get_path_length(following_path)
                    if current_corridor_length == None or following_path_length < current_corridor_length:
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
        # Check if we already have any corridor
        # If not, try to find a starting point (e.g. an already set door)
        if len(current_corridor_nodes) == 0:
            if len(already_set_door_points) > 0:
                # Get a sample set door in case we have doors and append it to the list of nodes
                sample_set_door_point = next(iter(already_set_door_points))
                current_corridor_nodes.append(sample_set_door_point)
        # In case we already have a node to start, we can solve the rest of the corridor from it
        if len(current_corridor_nodes) > 0:
            start_path = current_corridor # It may contain segments already, from the free space
            start_path_points = current_corridor_nodes
            start_nodes = [ nodes[point] for point in start_path_points ]
            start_non_redundant_nodes = [ node for node in start_nodes if not node['is_redundant'] ] # Maybe this is redundant? (ironically)
            start_rooms = set(sum([ node['rooms'] for node in start_non_redundant_nodes ],[]))
            start_available_paths = sum([ node['paths'] for node in start_non_redundant_nodes ],[])
            start_available_path_nodes = sum([ node['path_nodes'] for node in start_non_redundant_nodes ],[])
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
                node_point for node_point, node_data in non_redundant_nodes.items() if room_with_less_nodes in node_data['rooms']
            ]
            # It may happen that a node alone is enough to cover all rooms (and there are not set doors yet)
            # We must check it at this point, or it will add a random segment (path) which may be not necessary and take much space
            # If following path includes all rooms then it is a candidate to be the corridor
            # Note that here we do not check doors. This is because if we are here then it means there are not doors
            for node in nodes_from_room_with_less_nodes:
                contains_all_rooms = all(room in nodes[node]['rooms'] for room in required_rooms)
                if contains_all_rooms:
                    current_corridor = []
                    current_corridor_nodes = [node]
                    break
            # If we failed to find a single node corridor then try so set the corridor from every node
            if len(current_corridor_nodes) == 0:
                for node in nodes_from_room_with_less_nodes:
                    start_point = node
                    start_node = nodes[start_point]
                    start_rooms = set(start_node['rooms'])
                    start_path = current_corridor # It may contain segments already, from the free space
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

        # At this point we have the backbone of the corridor
        # Check there is something
        if len(current_corridor) == 0 and len(current_corridor_nodes) == 0:
            raise ValueError('Empty corridor')

        # Check the corridor is a unified path, and not more than one (i.e. all corridor segments are connected)
        # If the corridor is not connected, include the minimum amount of additional segments to connect them all
        # This may happen when the parent (self) free space is splitted by children rooms

        # Instead of checking all segments to be connected (which would take longer) check all nodes to be connected
        # To do so, check the connected segment of each node to be in the current corridor
        # If the connected segment is in the corridor then consider its corresponding node to be connected
        # connected_nodes = [ corridor_nodes[0] ]
        # already_checked_nodes = []
        # is_corridor_splitted = False
        # while len(connected_nodes) < len(corridor_nodes):
        #     # Get a new node which has not be checked yet
        #     sample_node = next((node for node in connected_nodes if node not in already_checked_nodes), None)
        #     # If there is no more nodes to check while there are still nodes to be connected then it means the corridor is broken
        #     if not sample_node:
        #         are_nodes_connected = True
        #         break
        #     # Find all connected nodes according to if their corresponding connected segment is in the current corridor
        #     sample_node_data = nodes[sample_node]
        #     for connected_segment, connected_node in zip(sample_node_data['connected_segments'], sample_node_data['path_nodes']):
        #         if connected_segment in corridor_segments:
        #             connected_nodes.append(connected_node)
        #     already_checked_nodes.append(sample_node)

        # # In case we have separated paths we must find the shortest path(s) to connect them
        # # Note that extra nodes may be not required, but only segments. See figure 5
        # if is_corridor_splitted:
        #     pass

        # Check the current corridor contains al nodes at this point
        current_rooms = set(sum([ nodes[point]['rooms'] for point in current_corridor_nodes ],[]))
        if not is_corridor_finished(current_rooms, current_corridor_nodes):
            raise ValueError('Failed to set the corridor')

        # Display the current corridor
        elements_to_display = [ segment.get_colored_segment('red') for segment in current_corridor ]
        self.update_display(extra=elements_to_display, title='Display the corridor backbone')

        # Return all the results inside a dict
        # It will include not only the backbone but also other intenral varibales which my be used further
        corridor_backbone = {
            'segments': current_corridor,
            'nodes' : current_corridor_nodes,
            # Save expensive internal vlaues to be forwarded to the set_corridor function
            'internal': (required_children_rooms, corridor_free_regions)
        }
        return corridor_backbone

    # Find the corridor backbone and then claim the area around the selected corridor segments to build the actual corridor
    # Then set the doors to connect the corridor with every child
    # Sometime additional area is claimed to grant door placement
    def set_corridor (self, verbose : bool = False) -> bool:
        
        # Find the corridor backbone
        corridor_backbone = self.find_corridor_backbone(verbose=verbose)

        # If the returned corridor is none then it menas there is no need to stablish a corridor
        if corridor_backbone == None:
            return True

        # Set other variables to be used
        # LORE: These variables are defined in the 'find_corridor_backbone' function
        # LORE: This function and the set_corridor function were one single function before

        # Get the exterior polygon
        exterior_polygon = self.boundary.exterior_polygon
        # Get all grids which may contribute to the corridor
        functional_corridor_grid = self.corridor_grid + self.free_grid
        # Get the actual corridor backbone
        current_corridor = corridor_backbone['segments']
        current_corridor_nodes = corridor_backbone['nodes']
        # Inherit some internal variables from the corridor backbone function
        required_children_rooms, corridor_free_regions = corridor_backbone['internal']

        # Build the corridor by claiming area around the corridor path
        # To do so we must set a boundary
        # With the current implementation there should never be cyclic (closed) corridors
        # For this reason the boundary should only have an exterior polygon
        # However if this happens in the future there should be no problem since boundaries have interior polygons as well

        # Set the corridor size (width)
        corridor_size = self.corridor_size

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

        # In case we have rigid rooms, set the rigid grid which must never be modified
        # Include these overlaps as out regions
        rigid_grid = None
        rigid_rooms = [ room for room in self.children if room.rigid ]
        if len(rigid_rooms) > 0:
            # Merge all rigid room grids
            rigid_grid = rigid_rooms[0].grid
            for room in rigid_rooms[1:]:
                rigid_grid += room.grid

        # Set a function to generate the corridor boundary
        # This process is wrapped in a function because we may have to change the corridor and redo the boundary further
        # e.g. a door can not be relocated in the boundary so it must be relocated now and the corridor will change
        def make_corridor_grid () -> Grid:
            corridor = Path(current_corridor)
            # Now we must substract segments in the free space (fake corridors)
            # Always leave a minimal corridor fragment from the free space for the grid to be continuous
            if functional_corridor_grid:
                # Iterate over the different free grids and remove its segments from the corridor
                for corridor_free_region in corridor_free_regions:
                    for segment in corridor_free_region['corridor_segments']:
                        corridor -= segment
                # Now that we have the corridor alone we can search for corridor fragments to connect
                for corridor_free_region in corridor_free_regions:
                    for segment in corridor_free_region['corridor_segments']:
                        # Skip segments not connected to the remaining corridor
                        if not corridor.is_connected_to_segment(segment):
                            continue
                        # Set a fragment to connect the corridor
                        connection_node = corridor.get_segment_connected_node(segment)
                        connection_segment = segment.get_side_segment(connection_node, corridor_size)
                        corridor += connection_segment
                        break

            # It may happen that there are no corridor segments at this point when all rooms are connected by a point
            if len(corridor) == 0:
                # Check that there is just one node, as expected
                if len(current_corridor_nodes) != 1:
                    raise SystemExit('There is something very wrong with the corridor backbone')
                corridor_grid = Grid([ generate_point_rect(current_corridor_nodes[0], corridor_size) ])
            else:
                # Display the current corridor
                elements_to_display = [ segment.get_colored_segment('red') for segment in corridor.segments ]
                self.update_display(extra=elements_to_display, title='Display the corridor to build the boundary')
                # Generate a boundary around the current corridor path
                corridor_grid = corridor.get_margined_grid(corridor_size)
            # Display the very first corridor boundary
            elements_to_display = corridor_grid.get_perimeter_segments('blue')
            self.update_display(extra=elements_to_display, title='Displaying first corridor boundary')

            # Set the regions of the corridor which must be removed
            # e.g. regions out of the parent boundary, in case it is not adaptable
            # e.g. regions over rigid rooms which must not be modified
            excluding_regions = Grid()

            # Find the regions of the corridor which are out of the parent (self)
            current_grid = self.grid if self.grid else exterior_polygon.grid
            out_regions = corridor_grid - current_grid

            # In case we have out regions...
            if out_regions and not self._child_adaptable_boundary:
                excluding_regions += out_regions

            # Check also if the corridor overlaps with rigid rooms (rooms which must not be modified)
            rigid_regions = rigid_grid.get_overlap_grid(corridor_grid) if rigid_grid else None

            # In case we have rigid regions we must exclude them
            if rigid_regions:
                excluding_regions += rigid_regions

            # If there are not excluding regions (not the usual case) then we are done
            if not excluding_regions:
                self.corridor_grid += corridor_grid
                return
            
            # Now we must substract excluding regions from the current corridor
            corridor_grid -= excluding_regions

            # Check the corridor has not been fully consumed
            if not corridor_grid:
                raise ValueError('Corridor was all in excluding regions')

            # Get the corridor boundary segments
            corridor_boundaries = corridor_grid.boundaries
            corridor_boundary_segments = sum([ boundary.segments for boundary in corridor_boundaries ], [])

            # Display the corridor boundaries after excluding regions removal
            elements_to_display = [ segment.get_colored_segment('blue') for segment in corridor_boundary_segments ]
            self.update_display(extra=elements_to_display, title='Displaying corridor boundaries after removing the excluding regions')

            # And now we must expand the corridor regions where we substracted the excluding regions
            # Otherwise the corridor would have regions which do not respect the minimum size
            for corridor_boundary in corridor_grid.boundaries:
                # IMPORTANT: Use the overlap with the path segments instead of the excluding region boundary segments
                # IMPORTANT: Otherwise, we may expand the corridor over unnecessary space. See figure 7
                # The region to be expanded is deducted from the segments in the path which overlap the already truncated corridor boundary
                excluded_reference_segments = corridor_boundary.get_segments_overlap_segments(current_corridor)
                # Get the corridor exterior polygon
                corridor_polygon = corridor_boundary.exterior_polygon
                # Once we have these segments we must "project" a corridor from them
                # This is like creating a corridor along the exterior polygon, which is fully inside of the polygon
                def all_inside (segment : Segment, direction : Vector) -> number:
                    # For the dead ends
                    # Note that for dead ends direction will always be equal to segment.direction, and not -segment.direction
                    if direction == segment.direction:
                        return 0
                    # For the inside
                    if direction == corridor_polygon.get_border_inside(segment):
                        return corridor_size
                    # For the outside
                    return 0
                # Display the segments used as reference for the excluded region
                elements_to_display = [ segment.get_colored_segment('green') for segment in excluded_reference_segments ]
                self.update_display(extra = elements_to_display, title = 'DEBUG: Excluded reference segments')
                # Generate the extension boundary from the excluded reference segments
                extension_corridor = Path(excluded_reference_segments)
                extension_grid = extension_corridor.get_margined_grid(all_inside)
                # Display the segments used as reference for the excluded region
                elements_to_display = extension_grid.get_perimeter_segments('purple')
                self.update_display(extra = elements_to_display, title = 'DEBUG: Extension segments')
                # elements_to_display = [ segment.get_colored_segment('purple') for segment in extension_boundaries[0].segments ]
                # self.update_display(extra=elements_to_display, title='Debug 2')
                # Now add the extended grid to the corridor grid
                # Note that both grids will always overlap
                corridor_grid += extension_grid

            # Display the corridor boundaries
            corridor_boundaries = corridor_grid.boundaries
            corridor_boundary_segments = sum([ boundary.segments for boundary in corridor_boundaries ], [])
            elements_to_display = [ segment.get_colored_segment('blue') for segment in corridor_boundary_segments ]
            self.update_display(extra=elements_to_display, title='Displaying corridor boundaries after expanding to compensate the removal of excluding regions')

            # Now add the free grid to the corridor grid
            # DANI: Esto tal vez se puede replantear. No es necesario y perdemos la free grid (espacio claimable)
            if self.free_grid:
                corridor_grid += self.free_grid
                # Display the corridor boundaries
                corridor_boundaries = corridor_grid.boundaries
                corridor_boundary_segments = sum([ boundary.segments for boundary in corridor_boundaries ], [])
                elements_to_display = [ segment.get_colored_segment('blue') for segment in corridor_boundary_segments ]
                self.update_display(extra=elements_to_display, title='Displaying corridor boundaries after adding back the parent free space')

            # Check we have only one boundary at this point
            if len(corridor_grid.boundaries) > 1:
                raise ValueError('We have a splitted corridor')

            # Remove regions from the corridor grid which are not respecting the minimum size
            # This may happen in some scenarios but there should be no problem in removing them. See figure 8
            corridor_grid = corridor_grid.keep_minimum(self.corridor_size)

            # Check the corridor has not been fully consumed
            if not corridor_grid:
                raise ValueError('Corridor was not respecting minimum size')

            # Now set the corridor grid officially
            self.corridor_grid += corridor_grid

            # Display the corridor boundaries
            corridor_boundaries = self.corridor_grid.boundaries
            corridor_boundary_segments = sum([ boundary.segments for boundary in corridor_boundaries ], [])
            elements_to_display = [ segment.get_colored_segment('blue') for segment in corridor_boundary_segments ]
            self.update_display(extra=elements_to_display, title='Displaying corridor boundaries after removing regions not respecting the corridor minimum size')

        # Run the function above to generate the corridor
        make_corridor_grid()

        # At this point there must be a corridor grid
        if not self.corridor_grid:
            raise SystemExit('Failed to set a corridor grid')

        # Set a function to substract the corridor region from the rest of child rooms
        # Note that, at this point, other rooms do not expand to compensate the lost area yet
        for child in self.children:
            if not child.truncate(self.corridor_grid, force=True, skip_update_display=True):
                raise ValueError(f'The space required by the corridor cannot be claimed from {child.name}')

        # Show the relocated doors
        self.update_display(title='Displaying early corridor')

        # Discarded spots may be generated while placing the early corridor
        # We must try to save those spaces which may be reclaimed by other rooms
        # It is not always possible to retrieve these regions however
        if self.discarded_grid:
            self.reassign_discarded_regions()

        # Get all children doors, both the already stablished and the not stablished ones
        # Add the room door also if this room is the root
        # Otherwise the door room must never be moved since it is already placed according to the parent corridor
        # The solving must be top-down, and it may be a node in the recently solved corridor so it is better not to move it
        door_rooms = self.children if self.parent else [ self ] + self.children
        doors = sum([ room.doors for room in door_rooms ], [])
        
        # Now expand the corridor as much as we need to make space for doors between the corridor and each connected room
        while True:
            # Track if there has been a corridor expansion
            # If so, we must restart the process and check again each door
            expanded = False
            # Locate doors in the new corridor boundary
            # Note that an already placed door may need to be relocated since its room may have been truncated by the corridor
            # In case a door has no space to be located we expand the corridor grid a bit in order to make space
            for door in doors:
                # Get the polygon where the door must be placed
                door_polygon = door.room.boundary.exterior_polygon
                # Set a couple of useful values
                corridor_polygon = self.corridor_grid.boundaries[0].exterior_polygon
                # Set a function to expand the corridor in reference to a given segment
                # The expansion segment is a segment in contact with the corridor
                def expand_corridor (expansion_segment : Segment) -> bool:
                    #print('EXPANDING CORRIDOR FOR ' + door.room.name + ' -> ' + str(expansion_segment))
                    nonlocal door_polygon
                    # Get a perpendicular segment as long as the corridor size and pointing inside the corridor (outside the door room)
                    # Do not rely on the other corner segment. In a diagonal supported scenario it may not be perpendicular
                    corridor_reference_segment = None
                    for segment in corridor_polygon.segments:
                        overlap = segment.get_overlap_segment(expansion_segment)
                        if overlap:
                            corridor_reference_segment = segment
                            break
                    if not corridor_reference_segment:
                        elements_to_display = [ s.get_colored_segment('blue') for s in corridor_polygon.segments ]
                        elements_to_display.append(expansion_segment.get_colored_segment('red'))
                        add_frame(elements_to_display, 'DEBUG')
                        raise ValueError('Cannot find the reference segment to expand corridor')
                    corridor_inside_direction = corridor_polygon.get_border_inside(corridor_reference_segment)
                    reference_point = expansion_segment.points[0] # Just any point in the segment
                    perpendicular_point = reference_point + corridor_inside_direction * corridor_size
                    perpendicular_segment = Segment(reference_point, perpendicular_point)
                    # Create the new Rect and then the new grid from it
                    # Note that this is not fully diagonal supported friendly beacuse the two segments may not create a classical 'Rect'
                    # However it should be easy to port once diagonals are supported
                    expansion_rect = Rect.from_segments([expansion_segment, perpendicular_segment])
                    expansion_grid = Grid([expansion_rect])
                    # Check the expansion grid to be inside the parent boundary
                    if not self._child_adaptable_boundary and expansion_grid not in self.grid:
                        return False
                    # If the expansion grid is colliding with the rigid grid we must stop
                    if rigid_grid and rigid_grid.get_overlap_grid(expansion_grid):
                        return False
                    # Now merge the expansion grid with the current corridor grid
                    self.corridor_grid += expansion_grid
                    # Substract the new corridor region from the children rooms
                    for child in self.children:
                        if not child.truncate(expansion_grid, force=True, skip_update_display=True):
                            raise ValueError(f'The space required by the corridor cannot be claimed from {child.name}')
                    # If there was a door in the expansion segment then we must relocate it
                    # DANI: Si no recolocas las puertas afectadas aquí se hace después igualmente
                    # DANI: Pero si te falla en algún punto antes de que se recoloque verás una puerta fuera de sitio
                    # DANI: Esto puede hacer parecer que hay un problema, así que mejor recolocar la puerta aquí
                    affected_doors = [ d for d in doors if d.point in expansion_rect ]
                    for affected_door in affected_doors:
                        affected_door.check_and_relocate()
                    # Update the corridor boundary segments and corners now that the corridor has been modified
                    corridor_boundary_segments = sum([ boundary.segments for boundary in self.corridor_grid.boundaries ], [])
                    # Also update the door polygon now that the door room may have been truncated
                    door_polygon = door.room.boundary.exterior_polygon
                    # Now that we expanded the grid there should be suitable space for the door to be located
                    available_segments = door_polygon.get_segments_overlap_segments(corridor_boundary_segments)
                    suitable_points = door.find_suitable_points(available_segments)
                    if len(suitable_points) == 0:
                        self.update_display(title='Debug')
                        raise RuntimeError(f'No suitable door space for room {door.room.name} after corridor expansion')
                    # Place the door in a suitable position if not already
                    # It should be already in a rigid door scenario
                    if not door.point or not door.is_point_suitable():
                        door.point = random.choice(suitable_points)
                    return True
                # If the door was placed already check it is still in a right place
                if door.point:
                    # If the door is already in both the corridor and its room boundaries then we do not need to relocate
                    door_in_room = door.margined_segment in door_polygon
                    door_in_corridor = any(door.margined_segment in boundary for boundary in self.corridor_grid.boundaries)
                    if door_in_room and door_in_corridor:
                        continue
                    # In case the door is rigid we must check it is over the door room and then expand the corridor directly
                    if door.rigid:
                        # If door is not over the room then there is something very wrong
                        if not door_in_room:
                            raise ValueError(f'Door from room {door.room.name} is rigid and not over the room')
                        # Expand corridor to reach the door
                        if expand_corridor(door.margined_segment):
                            # If we expanded successfully then skip checking other available segments and proceed
                            expanded = True
                            break
                # If the door has not been located yet or it needs to be relocated then we find the available space
                corridor_boundary_segments = sum([ boundary.segments for boundary in self.corridor_grid.boundaries ], [])
                available_segments = door_polygon.get_segments_overlap_segments(corridor_boundary_segments)
                # If the door may be relocated and there are avaibale segments then try to relocate the door
                if len(available_segments) > 0:
                    # Then find suitable space among the available space
                    suitable_points = door.find_suitable_points(available_segments)
                    if len(suitable_points) > 0:
                        door.point = random.choice(suitable_points)
                        continue
                # Before trying to expand the corridor, try to expand the room if it is not rigid
                # Corridor space is never returned while room space may be compensated further
                if not door.room.rigid:
                    room_expanded = False
                    # If there is at least one segment overlap between the corridor and the door room then try to expand these
                    if len(available_segments) > 0:
                        # Set a function to expand the room in case it is possible and return a boolean
                        # This function is called just once in the code, but we do it this way to get the boolean
                        def expand_room (available_segment : Segment) -> bool:
                            # Get the whole corridor segment where the available segment belongs
                            corridor_reference_segment = None
                            for segment in corridor_polygon.segments:
                                overlap = segment.get_overlap_segment(available_segment)
                                if overlap:
                                    corridor_reference_segment = segment
                                    break
                            if not corridor_reference_segment:
                                raise ValueError('Cannot find the reference segment to expand room')
                            if corridor_reference_segment.length < door.margined_width:
                                return False
                            # Now find in which direction we must expand this segment to do not get outside the corridor
                            # There should be always a point in common with the reference segment
                            common_point = next((point for point in available_segment.points if point in corridor_reference_segment.points), None)
                            if common_point == None:
                                raise ValueError(f'There is no common point between {door.room.name} expansion segment and the corridor reference segment')
                            other_point = next(point for point in available_segment.points if point != common_point)
                            # Now find the segment over the corridor which must be covered by the room expansion to fit the door
                            expansion_direction = (common_point + other_point).normalized()
                            expansion_length = door.margined_width - available_segment.length
                            expansion_segment = Segment(other_point, other_point + expansion_direction * expansion_length)
                            # Now find the whole region to be expanded 
                            corridor_inside_direction = corridor_polygon.get_border_inside(corridor_reference_segment)
                            perpendicular_point = other_point - corridor_inside_direction * door.room.min_size
                            perpendicular_segment = Segment(other_point, perpendicular_point)
                            expansion_rect = Rect.from_segments([expansion_segment, perpendicular_segment])
                            expansion_grid = Grid([expansion_rect])
                            # If the expansion grid is colliding with the rigid grid we must stop
                            if rigid_grid and rigid_grid.get_overlap_grid(expansion_grid):
                                return False
                            # If we can not expand the grid at this point then we must stop
                            # It may happen, for instance, when we try to expand over a door which can not be relocated anywhere else
                            if not door.room.expand_grid(expansion_grid, compensate_invaded=False, check_parent_free_grid=False):
                                return False
                            door.point = common_point + expansion_direction * (door.margined_width / 2)
                            return True
                        # ---- End of expand_room function
                        # Iterate over available segments trying to expand the room
                        for available_segment in available_segments:
                            if expand_room(available_segment):
                                room_expanded = True
                                break
                    # If there is no segment overlap between corridor and door room then we must rely in a point overlap
                    # This point will be a corner for both the corridor and the door room and it must always be there
                    else:
                        # If there is no segment overlap between corridor and door room then we must rely in a point overlap
                        # This point will be a corner for both the corridor and the door room and it must always be there
                        corridor_outside_corners = sum([ boundary.outside_corners for boundary in self.corridor_grid.boundaries ], [])
                        overlap_point = next((corner for corner in corridor_outside_corners if corner in door_polygon.corners), None)
                        if not overlap_point:
                            extra = [ s.get_colored_segment('blue') for s in corridor_boundary_segments ] + [ s.get_colored_segment('green') for s in door_polygon.segments ]
                            self.update_display(title='Debug', extra=extra)
                            raise RuntimeError(f'There is no overlap point between the corridor and the {door.room.name}. This should never happen.')
                        # Now find the two segments in the corridor boundary which include the overlap corner
                        corridor_reference_segments = [ segment for segment in corridor_polygon.segments if overlap_point in segment ]
                        # Set the expansion size as the whole door margined width
                        expansion_size = door.margined_width
                        for corridor_reference_segment in corridor_reference_segments:
                            # Get the other point (i.e. the point in the reference segment which is not overlapping with the door room)
                            other_point = next(point for point in corridor_reference_segment.points if point != overlap_point)
                            # Set the segment which covers the region that the corridor must cover
                            # Note that the expansion segment includes also the whole original corridor segment
                            # Most of the new space may overlap with already existing space
                            # However, in a diagonal supported scenario it may be critical to fill required space
                            expansion_direction = (overlap_point + other_point).normalized()
                            expansion_point = overlap_point + expansion_direction * expansion_size
                            expansion_segment = Segment(overlap_point, expansion_point)
                            # Set a function to expand the room in case it is possible and return a boolean
                            # This function is called just once in the code, but we do it this way to get the boolean
                            def expand_room (expansion_segment : Segment) -> bool:
                                # Now find the whole region to be expanded 
                                corridor_inside_direction = corridor_polygon.get_border_inside(corridor_reference_segment)
                                perpendicular_point = expansion_point - corridor_inside_direction * door.room.min_size
                                perpendicular_segment = Segment(expansion_point, perpendicular_point)
                                expansion_rect = Rect.from_segments([expansion_segment, perpendicular_segment])
                                expansion_grid = Grid([expansion_rect])
                                # If the expansion grid is colliding with the rigid grid we must stop
                                if rigid_grid and rigid_grid.get_overlap_grid(expansion_grid):
                                    return False
                                # If we can not expand the grid at this point then we must stop
                                # It may happen, for instance, when we try to expand over a door which can not be relocated anywhere else
                                if not door.room.expand_grid(expansion_grid, check_parent_free_grid=False, compensate_invaded=False):
                                    return False
                                door.point = overlap_point + expansion_direction * (door.margined_width / 2)
                                return True
                            # Now we have all we need to set the expansion using the previous function
                            if expand_room(expansion_segment):
                                # If we expanded successfully then skip checking other reference segments and proceed
                                room_expanded = True
                                break
                    # If there was a room expansion then we have to restart the process and check every door again
                    # If there was no expansion then we let the logic try to expand the corridor instead of the room
                    if room_expanded:
                        continue
                # If there is not suitable space then we have to expand a bit the corridor
                # At this point there should always be an overlap between the corridor and the door room
                # It may be just a point, in a few ocassions, but normally it will be a segment
                # In case we have segments in common, it normally will be one
                # Theorically, we should be always able to expand this segment to cover the minimum length in the door room
                # However, it may be not possible to the corridor to expand in some direction because of a conflict with a strict boundary
                # For this reason, we must try to expand the corridor in all possible directions before we surrender
                # Theorically* there should always work in at least one direction
                # * The only exception is a bottle neck made by a hand-set strict room
                # In case we cannot get the suitable space by expanding in any direction it is a fatal scenario and we must stop here
                if len(available_segments) > 0:
                    for available_segment in available_segments:
                        # Set the segment which covers the region that the corridor must cover
                        # Note that the expansion segment includes also the whole original corridor segment
                        # Most of the new space may overlap with already existing space
                        # However, in a diagonal supported scenario it may be critical to fill required space
                        expansion_segment = None
                        # Set the expansion size
                        # The expansion must cover the missing segment length to be wide enought to cover the door margined width
                        expansion_size = door.margined_width - available_segment.length
                        # We must find the direction to expand the corridor
                        # The available segment will always have one point which is an outside corner in the corridor boundary
                        # The expansion direction will be towards this point (from the other one)
                        corridor_outside_corners = sum([ boundary.outside_corners for boundary in self.corridor_grid.boundaries ], [])
                        towards_point = next((point for point in available_segment.points if point in corridor_outside_corners), None)
                        # This should never happen, theorically
                        if not towards_point:
                            raise ValueError('There is not a point which is in the corridor segment and another which is not, as expected')
                        # Get the other point
                        other_point = next(point for point in available_segment.points if point != towards_point)
                        # In case both points are in outside corners (it may happen if the corridor has been truncated) we must reconsider
                        if other_point in corridor_outside_corners:
                            # This is an all-proof solution
                            # We find the whole expansible segment in the door room where our available segment intersects
                            expansible_segment = None
                            # In case the door belongs to the parent (self) room...
                            if door.room == self:
                                # If the parent is child adaptable boundary just throw the segment we need and check it does not cut self further
                                if door.room._child_adaptable_boundary:
                                    # Create a fake grid all around self grid with a margin equal to the expansion size
                                    surrounding_rect = door_polygon.get_box(margin=expansion_size)
                                    surrounding_grid = Grid([surrounding_rect])
                                    expansible_grid = surrounding_grid - self.grid
                                    expansible_candidate_segments = expansible_grid.get_line_overlap_segments(available_segment.line)
                                    expansible_segment = next(( segment for segment in expansible_candidate_segments if towards_point in segment or other_point in segment), None)
                                # Set as expansible segment the whole segment in the polygon
                                else:
                                    expansible_segment = next(segment for segment in door_polygon if available_segment in segment)
                            # If the door room is a child then the whole space inside of its polygon is suitable for expansion
                            else:
                                expansible_candidate_segments = door_polygon.grid.get_line_overlap_segments(available_segment.line)
                                expansible_segment = next(( segment for segment in expansible_candidate_segments if towards_point in segment or other_point in segment), None)
                            # Check the expansible segment is suitable
                            if not expansible_segment:
                                raise ValueError('No expansible segment (this should never happen)')
                            if expansible_segment.length < door.margined_width:
                                raise ValueError('The expansible size is not enougth to allocate the door. Is room minimum size lower than the door margined size?')
                            # Now allocate a segment to be expanded with the required size in the expansible segment
                            # Easy: try to expand in one of the points. If it is enougth then we are done, if not try to expand by the other
                            first_point = expansible_segment.points[0]
                            second_point = expansible_segment.points[1]
                            def by_distance(point):
                                return first_point.get_distance_to(point)
                            sorted_available_points = sorted(available_segment.points, key=by_distance)
                            first_closer_point = sorted_available_points[0]
                            second_closer_point = sorted_available_points[1]
                            first_distance = first_point.get_distance_to(first_closer_point)
                            # If the first point is far enought to cover the expansion required then we expan at this site
                            if first_distance > expansion_size:
                                expansion_direction = (first_closer_point + first_point).normalized()
                                expansion_point = first_closer_point + expansion_direction * expansion_size
                                expansion_segment = Segment(expansion_point, second_closer_point)
                            # Otherwise, we claim the whole first segment and we expand the second as much as needed
                            else:
                                remaining_expansion_size = expansion_size - first_distance
                                expansion_direction = (second_closer_point + second_point).normalized()
                                expansion_point = second_closer_point + expansion_direction * remaining_expansion_size
                                expansion_segment = Segment(first_point, expansion_point)
                        # If we have the canonical scenario we simply push the outside corner
                        else:
                            expansion_direction = (other_point + towards_point).normalized()
                            expansion_point = towards_point + expansion_direction * expansion_size
                            expansion_segment = Segment(other_point, expansion_point)
                        # Now we have all we need to set the expansion using the previous function
                        if expand_corridor(expansion_segment):
                            # If we expanded successfully then skip checking other available segments and proceed
                            expanded = True
                            break
                    # If we expanded the corridor then stop here
                    # Note that the expansion may require to relocate again some doors which were already located
                    # We must start again to check all doors
                    if expanded:
                        # Show the expanded corridor
                        self.update_display(title='Displaying expanded corridor')
                        break
                    else:
                        raise RuntimeError(f'Failed to expand corridor for reaching room {door.room.name}')
                # If there is no segment overlap between corridor and door room then we must rely in a point overlap
                # This point will be a corner for both the corridor and the door room and it must always be there
                corridor_outside_corners = sum([ boundary.outside_corners for boundary in self.corridor_grid.boundaries ], [])
                overlap_point = next((corner for corner in corridor_outside_corners if corner in door_polygon.corners), None)
                if not overlap_point:
                    raise RuntimeError(f'There is no overlap point between the corridor and the {door.room.name} door. This should never happen.')
                # Now find the two segments in the corridor boundary which include the overlap corner
                reference_segments = [ segment for segment in corridor_boundary_segments if overlap_point in segment ]
                # Set the expansion size as the whole door margined width
                expansion_size = door.margined_width
                for reference_segment in reference_segments:
                    # Get the other point (i.e. the point in the reference segment which is not overlapping with the door room)
                    other_point = next(point for point in reference_segment.points if point != overlap_point)
                    # Set the segment which covers the region that the corridor must cover
                    # Note that the expansion segment includes also the whole original corridor segment
                    # Most of the new space may overlap with already existing space
                    # However, in a diagonal supported scenario it may be critical to fill required space
                    expansion_direction = (other_point + overlap_point).normalized()
                    expansion_point = overlap_point + expansion_direction * expansion_size
                    expansion_segment = Segment(other_point, expansion_point)
                    # Now we have all we need to set the expansion using the previous function
                    if expand_corridor(expansion_segment):
                        # If we expanded successfully then skip checking other reference segments and proceed
                        expanded = True
                        break
            # If there was a corridor expansion then we have to restart the process and check every door again
            # If there was no expansion then it means all doors al well placed and we can exit the process
            if not expanded:
                break
        # Show the relocated doors
        self.update_display(title='Displaying relocated doors')

        # Discarded spots may be generated while truncating rooms to place the corridor
        # We must try to save those spaces which may be reclaimed by other rooms
        # It is not always possible to retrieve these regions however
        if self.discarded_grid:
            self.reassign_discarded_regions()

        # Check there is enought space to fit children after corridor area substraction
        children_available_area = self.area - self.corridor_grid.area - self.discarded_grid.area
        children_required_area = sum([ child.min_area for child in self.children ])
        if children_available_area < children_required_area:
            if verbose: print(f'There is not enought space {children_available_area} to fit all children {children_required_area} after corridor area substraction')
            return False

        # Now, relocate and reshape children rooms
        for child in self.children:
            current_behaviour = 'conformist' if equal(self.free_area, 0) else 'greedy'
            if not child.fit_to_required_area(behaviour=current_behaviour, verbose=True):
                if verbose: print(f'{child.name} failed to fit to required area after corridor area truncation')
                return False

        # Show redistribution after reshaping child rooms
        self.update_display(title='Displaying redistributed rooms')
        return True

    # Reduce the number of corners in this room exterior polygon by reshaping self and children boundaries
    # This function is meant to run in the root room only
    # Return True if the reduction suceed or False if it failed
    def reduce_corners (self, verbose : bool = False) -> bool:
        if verbose: print(f'Reducing {self.name} corners')
        if self.rigid:
            raise RuntimeError('Trying to reduce corners in a rigid room')
        # Get the exterior polygon
        exterior_polygon = self.boundary.exterior_polygon
        while len(exterior_polygon.corners) > self.max_corners:

            # Get self zigzags
            zigzags = get_polygon_zigzags(exterior_polygon)
            
            # Check zigzag segments to be modified are not in contact with rigid rooms which must not be modified
            forbidden_segments = list(set(sum([ child.boundary.exterior_polygon.segments for child in self.children if child.rigid ], [])))
            # Add corridor contact segments to the list of forbidden segments
            for boundary in self.corridor_grid.boundaries:
                forbidden_segments += boundary.segments
            def has_rigid_conflict (zigzag : dict) -> bool:
                modified_segments = [ zigzag['inside_segment'], zigzag['middle_segment'], zigzag['outside_segment'] ]
                for modified_segment in modified_segments:
                    for forbidden_segment in forbidden_segments:
                        overlap = modified_segment.get_overlap_segment(forbidden_segment)
                        if overlap:
                            return True
                return False

            # Filter out those zigzags with conflicts
            zigzags = [ zigzag for zigzag in zigzags if not has_rigid_conflict(zigzag) ]

            # Try to remove corners in the most suitable zigzag
            # If it fails, try with the next one
            # WARNING:
            # Functions which may fail are smart enought to recover boundary backups in case of failure
            # However we must backup self and children boundaries, since we may succeed in several steps and then fail later
            succeed = False
            for zigzag in zigzags:
                backup = self.make_grid_backup()
                # Push the inside segment and pull the outside segment
                # The push and pull lengths must be calculated to make both segments match while the parent area is kept
                area = zigzag['middle_segment'].length * zigzag['outside_segment'].length
                new_segment_length = zigzag['inside_segment'].length + zigzag['outside_segment'].length
                inside_push_length = area / new_segment_length
                outside_pull_length = zigzag['middle_segment'].length - inside_push_length
                # Push before pull, so we have enought area to recover after the pull
                if not self.push_boundary_segment(zigzag['inside_segment'], inside_push_length):
                    if verbose: print(f'  Failed to reduce corner: Failed to push {zigzag["inside_segment"]}')
                    # There is no need to recover the backup at this point
                    continue
                # Now pull the boundary
                # Note tha this will truncate any children in the pulled area automatically
                # In this case we do not check the parent minimum size to be respected
                # Note that the pulled segment will always belong to a non-rigid brother
                # Thus the free space created will be claimed right away
                if not self.pull_boundary_segment(zigzag['outside_segment'], outside_pull_length, force_child_truncation=True):
                    if verbose: print(f'Failed to reduce corner: Failed to pull {zigzag["outside_segment"]}')
                    # Revert the previous push
                    self.restore_grid_backup(backup)
                    continue
                # Relocate children to fit in the new boundary
                truncated_children = [ child for child in self.children if not child.is_fit_to_required_area()  ]
                child_conflict = False
                for child in truncated_children:
                    if not child.fit_to_required_area():
                        child_conflict = True
                        break
                # If there was a failure during children relocation then restore the boundary backups and proceed to the next zigzag
                if child_conflict:
                    if verbose: print(f'Failed to reduce corner: Something went wrong while refitting {child.name}')
                    self.restore_grid_backup(backup)
                    continue
                # If everything was fine then stop here
                # Only 1 zigzag may be solved at once
                succeed = True
                break

            # If we succeeded to remove one of the corners then we continue
            if succeed:
                # Recalculate the exterior polygon
                exterior_polygon = self.boundary.exterior_polygon
                continue

            # If we failed to remove the corner with all the available zigzags then we must surrender at this point
            current_corners = len(exterior_polygon.corners)
            if verbose: print(f'Failed to reduce corners to {self.max_corners}. Current number: {current_corners}')
            return False
        return True

    # Reassign discarded regions to colliding rooms
    # Ignore rooms area budget
    def reassign_discarded_regions (self):
        # Set if there was at leas 1 reasignment
        any_reasigned = False
        # Keep trying to find new regions to reasssign until there is no more discarded grid
        # This loop will also break if we failed to reassign any region
        while self.discarded_grid:
            # Set if there was any reassignment
            reasigned = False
            # Iterate over each regions in the discarded grid
            for discarded_boundary in self.discarded_grid.boundaries:
                if self._find_and_reassign_discarded_boundary_region(discarded_boundary):
                    reasigned = True
                    # Show the progress
                    # self.update_display(title='Displaying discarded grid reasignation step')
            # If there was a reasignment then set the global as true
            if reasigned:
                any_reasigned = True
            # If there was no reasignment then stop here
            else:
                break
        return any_reasigned

    # Given a boundary which is meant to belong to the discarded grid,
    # Find one region which may be claimed by the corridor and claim it
    # DANI: Esta función NO está pensada para diagonales
    # LORE: This function had more functionalities when the discarded grid was supervised but now it is self calculated
    def _find_and_reassign_discarded_boundary_region (self, discarded_boundary : Boundary) -> bool:
        # Get the discarded exterior polygon
        discarded_polygon = discarded_boundary.exterior_polygon
        # Get all segments between corridor and discarded grid
        corridor_segments = []
        for corridor_boundary in self.corridor_grid.boundaries:
            corridor_segments += corridor_boundary.exterior_polygon.get_polygon_overlap_segments(discarded_polygon)
        # If any segment is wide enough to be claimed then go ahead
        corridor_wide_segment = next((segment for segment in corridor_segments if equal_or_greater(segment.length, self.corridor_size)), None)
        if corridor_wide_segment:
            # Set the grid to be claimed by the corridor and claim it
            corridor_region = discarded_boundary.get_border_projection(corridor_wide_segment)
            self.corridor_grid += corridor_region
            return True
        # Find corners along corridor segments
        corridor_corner = None
        for current, others in afterwise(corridor_segments):
            for other in others:
                # Find out if these segments connect
                connecting_point = current.get_connecting_point(other)
                if not connecting_point:
                    continue
                # If so, make sure the corner is well oriented
                corner = discarded_polygon.get_corner(connecting_point)
                if not corner.inside:
                    corridor_corner = corner
                    break
            # If we already found a corner then stop here
            if corridor_corner:
                break
        # If we find a corner then claim it
        if corridor_corner:
            # Set the grid to be claimed by the corridor
            # Note that this is not as simple as creating a rect from these 2 segments
            # We could ignore and interior polygon of the boundary and create a free spot which does not respect the minimum
            corridor_corner_segments = corridor_corner.segments
            corridor_region = discarded_boundary.get_border_projection(corridor_corner_segments[0])
            corridor_region += discarded_boundary.get_border_projection(corridor_corner_segments[1])
            # Claim it to the corridor
            self.corridor_grid += corridor_region
            return True
        # If we failed to find any of the previous features then we return false
        return False

    # Reduce the number of corners in children exterior polygons by reshaping children boundaries
    # This can be done only in a very specific situation:
    #   There must be a zigzag whose segments are fully covered by only one room on each side
    #   These rooms must be not rigid and not the parent room
    # Note that this function works in an environment which is not much flexible, so its changes will be limited
    # It will not return true or false, since this function is expected to fail at some point before reaching the 4 corners
    def reduce_children_corners (self, verbose : bool = False):
        # Iterate over children rooms
        for child in self.children:
            if verbose: print(f'Reducing corners in {child.name}')
            # Skip rigid children
            if child.rigid:
                continue
            exterior_polygon = child.boundary.exterior_polygon
            # DANI: Podría ser > child.max_corners en lugar de > 4, pero eso implicaría que el 4 fuese el por defecto
            # DANI: O sino implicaría tener que especificar que quieres 4 corners en todos los children
            while len(exterior_polygon.corners) > 4:
                # Get self zigzags
                zigzags = get_polygon_zigzags(exterior_polygon)
                # Check zigzag segments to be suitable
                # i.e. all segments fully overlap with the same brother non rigid room
                free_frontiers, brother_frontiers, parent_frontiers = child.get_frontiers()
                suitable_frontiers = free_frontiers + brother_frontiers
                def is_suitable (zigzag : dict) -> bool:
                    # Check contact with parent or rigid rooms
                    modified_segments = [ zigzag['inside_segment'], zigzag['middle_segment'], zigzag['outside_segment'] ]
                    # Check the three segments are among free/conflict frontiers and their frontier rooms match
                    common_brother_room = None
                    for modified_segment in modified_segments:
                        equivalent_frontier = next(( frontier for frontier in suitable_frontiers if frontier == modified_segment ), None)
                        if not equivalent_frontier:
                            return False
                        if len(equivalent_frontier.rooms) > 1:
                            return False
                        frontier_room = equivalent_frontier.rooms[0]
                        if common_brother_room:
                            if common_brother_room != frontier_room:
                                return False
                        else:
                            common_brother_room = frontier_room
                    # Add the common room to the dit content
                    zigzag['room'] = common_brother_room
                    return True

                # Filter out those zigzags with conflicts
                zigzags = [ zigzag for zigzag in zigzags if is_suitable(zigzag) ]

                # Try to remove corners in the most suitable zigzag
                # If it fails, try with the next one
                # WARNING:
                # Functions which may fail are smart enought to recover boundary backups in case of failure
                # However we must backup the self boundary, since we may succeed in several steps and then fail later
                succeed = False
                for zigzag in zigzags:
                    if verbose: print(f'Reducing corner with middle segment {zigzag["middle_segment"]}')
                    backup = self.make_grid_backup()
                    # Get the zigzag segments
                    inside_segment = zigzag['inside_segment']
                    middle_segment = zigzag['middle_segment']
                    outside_segment = zigzag['outside_segment']
                    # Push the inside segment and pull the outside segment
                    # The push and pull lengths must be calculated to make both segments match while the parent area is kept
                    area = middle_segment.length * outside_segment.length
                    new_segment_length = inside_segment.length + outside_segment.length
                    inside_push_length = area / new_segment_length
                    outside_pull_length = middle_segment.length - inside_push_length
                    # Pull the outside segment
                    # Pull before push, so we have enought area to recover after the push
                    if not child.pull_boundary_segment(outside_segment, outside_pull_length, check_parent_free_grid=False, verbose=verbose):
                        if verbose: print(f'  Failed to pull corner outside segment {outside_segment}')
                        # Before we give up we try to pull a more conservative distance
                        # It may happen that the resulting grid is not respecting the minimum size
                        # This will never happen if we pull the segment to align it to its closer paralel neighbour segment
                        other_outside_corner_point = next( point for point in outside_segment.points if point != zigzag['outside_corner'] )
                        other_outside_corner = exterior_polygon.get_corner(other_outside_corner_point)
                        # If the other corner is an inside corner then we can not do the trick
                        if other_outside_corner.inside:
                            # There is no need to recover the backup at this point
                            continue
                        other_outside_segment = next( segment for segment in other_outside_corner.segments if segment != outside_segment )
                        # The safe pull distance cannot be larger that the default distance if this was the problem
                        # In this case, it means the previous pull failed for other reason and not the minimum size
                        safe_pull_length = other_outside_segment.length
                        if safe_pull_length >= outside_pull_length:
                            # There is no need to recover the backup at this point
                            continue
                        # Now if we reduce the pull length we must also reduce the push length
                        reduction = safe_pull_length / outside_pull_length
                        inside_push_length = inside_push_length * reduction
                        # Now we are ready to try to pull again
                        # If it fails again then the problem was not the minimum size
                        if verbose: print('     Retrying safe pull')
                        if not child.pull_boundary_segment(outside_segment, safe_pull_length):
                            if verbose: print('     Safe pull failed as well')
                            # There is no need to recover the backup at this point
                            continue
                    # Now push the inside segment
                    if not child.push_boundary_segment(inside_segment, inside_push_length, easy=True, check_parent_free_grid=False, verbose=verbose):
                        if verbose: print(f'  Failed to push corner inside segment {inside_segment}')
                        # Revert the previous push
                        self.restore_grid_backup(backup)
                        continue
                    # Relocate children to fit in the new boundary
                    truncated_children = [ child for child in self.children if not child.is_fit_to_required_area()  ]
                    child_conflict = False
                    for child in truncated_children:
                        # Save a backup of the current child in case we have to recover its boundary later
                        child_boundary_backup = child.boundary
                        if not child.fit_to_required_area():
                            if verbose: print(f'Something went wrong while refitting {child.name}')
                            child_conflict = True
                            break
                        # Now add the child boundary to the backup
                        # Note that this is not done before since in case of failure the current children is backuped already
                        backup[child] = child_boundary_backup
                    # If there was a failure during children relocation then restore the boundary backups and proceed to the next zigzag
                    if child_conflict:
                        self.restore_grid_backup(backup)
                        continue
                    # If everything was fine then stop here
                    # Only 1 zigzag may be solved at once
                    succeed = True
                    break

                # If we succeeded to remove one of the corners then we continue
                if succeed:
                    # Recalculate the exterior polygon
                    exterior_polygon = child.boundary.exterior_polygon
                    continue

                # If there are not more suitable zigzags or the current ones failed to collapse then we stop here
                break

    # Get the minimum of all children minimum sizes
    # Get to the the root and the check all children min sizes recuersively in order to get the minimum
    # This function is meant to be used only once by the root, so its value is not stored
    # WARNING: 0 values are removed
    def get_children_min_min_size (self) -> number:
        all_min_sizes = [ child.min_size for child in self.children if child.min_size > 0 ]
        if len(all_min_sizes) == 0:
            return 0
        return min(all_min_sizes)

    # Get the minimum size of this room and its children recurisvely
    # Return the lowest minimum size
    def get_min_size_recursive (self) -> number:
        all_rooms = self.get_rooms_recuersive()
        all_min_sizes = [ room.min_size for room in all_rooms if room.min_size > 0 ]
        if len(all_min_sizes) == 0:
            return 0
        return min(all_min_sizes)

    # Calculate how much area we need to expand or contract
    # The behaviour argument sets if the area is conformist or exigent:
    # - exigent: aim for the exact target area
    # - conformist: required area is 0 as long as the current area is between the maximum and minimum areas
    # - greedy: become conformist only if the current area is above the target area
    # - humble: become conformist only if the current area is below the target area
    def get_required_area (self, behaviour : str, verbose : bool = False) -> number:
        if verbose:
            print(f'Setting required area for {self.name}')
            print(f'  Area range: {self.max_area} - {self.min_area}')
            print(f'  Target area: {self.target_area}')
            print(f'  Current area: {self.area}')
        # If the behaviour is exigent then aim for the exact target area
        # This is the usual behaviour when the corridor is not yet set
        if behaviour == 'exigent':
            return resolute(self.target_area - self.area)
        # If the behaviour is not exigent then target area is not a specific value but a range of values
        # This is the usual behaviour when the corridor is already set
        if behaviour == 'conformist':
            # If we are over the maximum we must reduce the area so we return a negative area
            if self.area > self.max_area:
                return resolute(self.max_area - self.area)
            # If area is below the minimum we must expand its area
            if self.area < self.min_area:
                return resolute(self.min_area - self.area)
            # If the current area is exactly within the target area range then we are done
            return 0
        # If it is greedy then request more area to reach the target area only
        if behaviour == 'greedy':
            # If we are over the maximum we must reduce the area so we return a negative area
            if self.area > self.max_area:
                return resolute(self.max_area - self.area)
            # If area is below the minimum we must expand its area
            if self.area < self.target_area:
                return resolute(self.target_area - self.area)
            # If the current area is exactly within the target area range then we are done
            return 0
        # If it is humble then request more area to reach the target area only
        if behaviour == 'humble':
            # If we are over the maximum we must reduce the area so we return a negative area
            if self.area > self.target_area:
                return resolute(self.target_area - self.area)
            # If area is below the minimum we must expand its area
            if self.area < self.min_area:
                return resolute(self.min_area - self.area)
            # If the current area is exactly within the target area range then we are done
            return 0
        raise ValueError(f'Not supported behaviour: {behaviour}')        

    # Check if this room is already fit to its required area
    def is_fit_to_required_area (self, behaviour : str) -> bool:
        required_area = self.get_required_area(behaviour=behaviour)
        # print(f'Fitting {self.name} to target area: {self.min_area} - {self.max_area}')
        # print(f'  Current area: {self.area} -> Required area: {required_area}')
        # print(f'  Is fit? {abs(required_area) < minimum_resolution * self.max_area}')
        return abs(required_area) < minimum_resolution * self.max_area

    # Expand or contract this room until it reaches the forced area
    # In case it is not able to fit at some point recover the original situation
    # Restricted segments are segments which must remain as are
    def fit_to_required_area (self, restricted_segments : list = [], behaviour : str = 'exigent', verbose : bool = False) -> bool:
        if verbose: print(f'Fitting {self.name}')
        # A room with no initial boundary/grid is not to be fitted
        if not self.grid:
            raise RuntimeError(f'Trying to fit {self.name} but it has no initial grid/boundary')
        if self.rigid:
            raise RuntimeError(f'Trying to fit {self.name} but it is rigid')
        # Calculate how much area we need to expand
        required_area = self.get_required_area(behaviour=behaviour, verbose=verbose)
        # If the area is already satisfied then stop here
        if self.is_fit_to_required_area(behaviour=behaviour):
            return True
        # If more area is required and this is an exigent call then make sure there is free space available to claim
        if required_area > 0 and behaviour == 'exigent' and not self.parent.free_grid:
            raise RuntimeError(f'Trying to fit {self.name} with {required_area} required area (exigent behaviour) but there is no free space available')
        # Make a boundary backup of this room and all its brothers
        backup = self.parent.make_grid_backup(children_only=True)
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
        while not self.is_fit_to_required_area(behaviour=behaviour):
            if verbose: print(f' Required area: {required_area}')
            # If the required are is positive it means we must expand our room
            if required_area > 0:
                # Set a function to supervise if each expansion step is succesful or not
                def expand_step () -> bool:
                    # Get the most suitable frontier to expand and try to expand it
                    # If the expansions fails, try with the next one
                    for frontier, loan_permission in self.get_best_frontiers(restricted_segments, verbose=verbose):
                        if self.expand_frontier(frontier, required_area, loan_permission, behaviour=behaviour, verbose=verbose):
                            return True
                    return False
                # Expand
                if expand_step():
                    # Refresh how much area we need to expand and go for the next step
                    required_area = self.get_required_area(behaviour=behaviour, verbose=verbose)
                    continue
            # If the required are is negative it means we need to contract our room
            if required_area < 0:
                # Set a function to supervise if each expansion step is succesful or not
                def contract_step () -> bool:
                    # Get the most suitable frontier to contract and try to contract it
                    # If the contraction fails, try with the next one
                    # DANI: De momento uso la misma lógica que la de la expansión porque no va mal
                    # DANI: i.e. evitar contraer fronteras del padre y priorizar fronteras libres es bueno
                    for frontier, loan_permission in self.get_best_frontiers(restricted_segments, contraction=True, verbose=verbose):
                        if self.contract_frontier(frontier, -required_area, verbose=verbose):
                            return True
                    return False
                # Contract
                if contract_step():
                    # Refresh how much area we need to expand and go for the next step
                    required_area = self.get_required_area(behaviour=behaviour, verbose=verbose)
                    continue
            # In case the fitting failed restore the backup and exit
            if verbose: print(f' Failed to fit {self.name}')
            self.parent.restore_grid_backup(backup)
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
        verbose : bool = False,
        ) -> Generator[ Tuple[ Segment, bool ], None, None ]:
        if verbose: print(f'Finding best frontiers for {self.name}')

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
            return frontier.length >= self.min_size or frontier.a in inside_corners or frontier.b in inside_corners

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
                    next_point = next(p for p in connected_frontier.points if p != next_point)
                # If the combined frontier reaches the minimum length add it to the list to be returned
                if combined_frontier.length >= self.min_size:
                    # Get unique implicated rooms
                    implicated_rooms = unique(implicated_rooms)
                    # Create a fit combined frontier which takes the minimum possible part from the last added frontier
                    # This way the expansion takes as much possible from the main frontier room
                    extended_direction = (point + next_point).normalized()
                    oposite_point = next(p for p in frontier.points if p != point)
                    fit_combined_frontier = Segment(oposite_point, oposite_point + extended_direction * self.min_size)
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
        # If there is no free space then the room with more area to give away in proportion
        # Then use this value to set the "score" of each brother room frontiers and sort them
        def sort_by_shortest_path (frontiers_group : list, verbose : bool = False) -> list:
            if verbose: print('Sorting best frontiers')
            # Find out if there is free space
            is_free_space_available = bool(self.parent.free_grid)
            if verbose: print(f'  There is free space: {is_free_space_available}')
            # If there is no free space then we must find the most suitable room to give away some area
            donnor_room = None
            if not is_free_space_available:
                print('WARNING: We are running a part of the code which is implemented but not tested deeply')
                # Find the brother room with the highest donable area proportion
                all_brothers = self.get_brother_rooms()
                donable_area_proportions = [ room.donable_area_proportion for room in all_brothers ]
                max_donable_area_proportion = max(donable_area_proportions)
                donnor_room_index = donable_area_proportions.index(max_donable_area_proportion)
                donnor_room = all_brothers[donnor_room_index]
                if verbose: print(f'  Donnor room: {donnor_room.name}')
            # Find the immediate colliding rooms
            colliding_rooms = unique([ frontier.rooms[0] for frontier in frontiers_group ])
            # WARNING: We must shuffle the colliding rooms at this point
            # Otherwise frontiers from the same room are always returned first
            # This has been observed experimentally
            random.shuffle(colliding_rooms)
            meaningful_frontiers = []
            for colliding_room in colliding_rooms:
                # If this is the donnor room already then set the scores here already
                if colliding_room == donnor_room:
                    colliding_room_frontiers = [ frontier for frontier in frontiers_group if frontier.rooms[0] == colliding_room ]
                    for frontier in colliding_room_frontiers:
                        frontier.score = 0
                        meaningful_frontiers.append(frontier)
                    continue
                # Otherwise run the recursive logic
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
                        # If there is free space then we finish when we find free space
                        if is_free_space_available:
                            if len(free_frontiers) > 0:
                                searching_free = False
                                break
                        # If there is no free space then we finish when we find the target room
                        else:
                            # Check if any of the brother frontiers is the donnor room
                            brother_rooms = set(sum([ frontier.rooms for frontier in brother_frontiers ], []))
                            if donnor_room in brother_rooms:
                                searching_free = False
                                break
                        current_frontiers += brother_frontiers
                    if not searching_free:
                        if verbose: print(f'  Room { colliding_room.name} -> Score {counter}')
                        break
                    previous_rooms = previous_rooms + current_rooms
                    current_rooms = unique([ frontier.rooms[0] for frontier in current_frontiers if frontier.rooms[0] not in previous_rooms ])
                    if len(current_rooms) == 0:
                        # Froniers which do not lead to free space are meaningless
                        # Trying to expand through them is a waste of time so they will be excluded
                        counter = None
                        if verbose: print(f'  Room { colliding_room.name} -> Dead end')
                        break
                    counter += 1
                # Set the scores for all frontiers
                colliding_room_frontiers = [ frontier for frontier in frontiers_group if frontier.rooms[0] == colliding_room ]
                if counter:
                    for frontier in colliding_room_frontiers:
                        frontier.score = counter
                        meaningful_frontiers.append(frontier)
            # Now sort frontiers using previous scores
            def by_score (frontier):
                return frontier.score
            return sorted(meaningful_frontiers, key=by_score)

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
            brother_frontiers = sort_by_shortest_path(brother_frontiers, verbose=verbose)
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
        if verbose: print('WARNING: There are no more frontiers available')

    # Push a segment in the boundary
    # Check everything is fine after the push and, if so, return True
    # In case there is any problem the push is not done and this function returns False
    # If the 'is_loaned' argument is passed then this room area is fitted after brother rooms truncation, before compensation
    # The behaviour argument is forwarded to the 'fit_to_required_area' function
    # The easy argument is forwarded to the 'truncate' function
    def push_boundary_segment (self,
        segment : Segment,
        push_length : number,
        is_loaned : bool = False,
        behaviour : str = 'exigent',
        easy : bool = False,
        check_parent_free_grid : bool = True,
        verbose : bool = False,
    ) -> bool:
        if verbose: print(f'Pushing boundary segment {segment}')
        # Get the push direction
        direction = -self.boundary.exterior_polygon.get_border_inside(segment)
        # If the push length at this point is 0 or close to it then we can not push
        if push_length < minimum_resolution:
            if verbose: print('  Failed to push: Push length does not reach the minimum resolution')
            return False
        # Get the parent room
        parent_room = self.parent
        # Check if the segment to be pushed is in the parent boundary
        # If it is, then we can not push it
        # You must push the parent boundary first
        if parent_room and next(parent_room.boundary.exterior_polygon.get_segment_overlap_segments(segment), None):
            if verbose: print('  Failed to push: Push segment is in the parent boundary. Push the parent first.')
            return False
        # Create the new rect with the definitive length
        new_point = segment.a + direction.normalized() * push_length
        new_side = Segment(segment.a, new_point)
        new_rect = Rect.from_segments([segment, new_side])
        # Make a grid out of the new rect and expand the grid
        new_grid = Grid([new_rect])
        return self.expand_grid(new_grid, is_loaned=is_loaned, behaviour=behaviour, easy=easy, check_parent_free_grid=check_parent_free_grid)
    
    # Pull a segment in the boundary
    # Check everything is fine after the pull and, if so, return True
    # In case there is any problem the pull is not done and this function returns False
    def pull_boundary_segment (self,
        segment : Segment,
        pull_length : number,
        force_child_truncation : bool = False,
        check_parent_free_grid : bool = True,
        verbose : bool = False,
    ) -> bool:
        if verbose: print(f'Pulling boundary segment {segment}')
        # Get the pull direction
        direction = self.boundary.exterior_polygon.get_border_inside(segment)
        # If the pull length at this point is 0 or close to it then we can not pull
        if pull_length < minimum_resolution:
            if verbose: print('  Failed to pull: Pull length does not reach the minimum resolution')
            return False
        # Create the new rect with the definitive length
        new_point = segment.a + direction.normalized() * pull_length
        new_side = Segment(segment.a, new_point)
        new_rect = Rect.from_segments([segment, new_side])
        # Make a grid out of the new rect
        removed_region = Grid([new_rect])
        # Truncate self grid
        return self.truncate(removed_region, force=force_child_truncation, check_parent_free_grid=check_parent_free_grid, verbose=verbose)

    # Try to expand a specific room frontier
    # Note that the frontier must contain the room it belongs to
    # Set the required (maximum) area it can expand
    # Return True if the expansion was succesful or False if there was no expansion
    def expand_frontier (self,
        frontier : Segment,
        required_area : number,
        allowed_loan_push : bool = False,
        behaviour : str = 'exigent',
        verbose : bool = False,
    ) -> bool:
        if verbose: print(f'Expanding frontier {frontier}')
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
            # Get the next room grid
            if next_room == parent_room:
                # If the room is the parent and it is child adaptable we must generate a new free grid around the parent (exterior grid)
                if parent_room._child_adaptable_boundary:
                    next_grid = parent_room.generate_exterior_free_grid(margin_size=self.max_size*2)
                else:
                    next_grid = next_room.free_grid + next_room.discarded_grid
            else:
                next_grid = next_room.grid
            # Add the current grid to the previous accumulated grid, if any
            if not grid:
                grid = next_grid
                continue
            grid = grid.get_merge_grid(next_grid, check_overlaps=False)
        # Find the maximum rectangles which are in contact with our frontier
        max_rects = grid.max_rects
        contact_max_rects = [ max_rect[0] for max_rect in max_rects if frontier in max_rect[0] ]
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
            print(f'Problem when expanding {self.name} frontier at {frontier}')
            raise RuntimeError(f'Empty grid: {grid}')

        # One and only one of the rows/columns will always include the segment
        space = next((rect for rect in rects if frontier in rect), None)
        if space == None:
            # If this happens it may mean there is a problem with the grid
            add_frame(rects +  [ frontier.get_colored_segment('red') ], 'Debug')
            room_names = ', '.join([ room.name for room in rooms ])
            raise RuntimeError('Frontier ' + str(frontier) + ' has no space in rooms ' + room_names)
        space_contact = next(segment for segment in space.segments if frontier in segment)
        space_forward_limit = space.get_size()[forward]

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
        def push_segment (pushed_segment : Segment, protocol : int = 1, verbose : bool = False) -> bool:
            if verbose:
                if protocol == 1:
                    print(f'Greedy push of {pushed_segment}')
                elif protocol == 2:
                    print(f'Moderate push of {pushed_segment}')
                elif protocol == 3:
                    print(f'Loaned push of {pushed_segment}')
            # Set the starting push length according to the required area and the segment length
            push_length = required_area / pushed_segment.length
            # If the push length at this point is 0 or close to it then we can not push
            # Try to reduce the segment
            if push_length < minimum_resolution:
                if verbose: print(f'WARNING: The push length is too small for segment {pushed_segment}')
                new_point_a = pushed_segment.a
                new_point_b = pushed_segment.a + pushed_segment.direction * self.min_size
                pushed_segment = Segment(new_point_a, new_point_b)
                if verbose: print('WARNING: segment has been reduced to {pushed_segment}')
                push_length = required_area / self.min_size
            # In case this is an insider segment which is not wide enought to be pushed alone,
            # Find out how much we can push this segment
            # i.e. find the connected frontier/s and get the maximum length of these segments
            corner_push_limit = None
            if lower(pushed_segment.length, self.min_size):
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
                # If the push length is "equal" to the maximum forward limit then we set the push length as the limit
                # This may seem redundant, but in some ocasions it solves resolution problems
                # The push length may be not identical to the maximum forward limit, but very similar
                elif equal(push_length, maximum_forward_limit):
                    push_length = maximum_forward_limit
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
            # Now that we have the definitive push length, we actually push the segment
            is_loaned = protocol == 3
            if not self.push_boundary_segment(pushed_segment, push_length, behaviour=behaviour,
                is_loaned=is_loaned, verbose=verbose, check_parent_free_grid=False):
                # If the push failed with a greedy push then try a conservative push
                if protocol == 1:
                    return push_segment(pushed_segment, 2, verbose=verbose)
                # Otherwise we surrender to push this segment
                return False
            # If the push succeeded
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
            return push_segment(frontier, push_protocol, verbose=verbose)

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

        return push_segment(reduced_frontier, push_protocol, verbose=verbose)

    # Try to contract a specific room frontier
    # Note that the frontier must contain the room it belongs to
    # Set the required (maximum) area it can contract
    # Return True if the contraction was succesful or False if there was no contraction
    def contract_frontier (self, frontier : Segment, required_area : number, verbose : bool = False) -> bool:
        if verbose: print(f'Contracting frontier {frontier}')
        # Get the parent room
        parent_room = self.parent
        # Get the exterior polygon of the room boundary
        exterior_polygon = self.boundary.exterior_polygon
        # Get the inside corners
        inside_corners = [ corner for corner in exterior_polygon.corners if corner.inside == True ]
        # Find the maximum rectangles which are in contact with our frontier
        grid = self.grid
        max_rects = grid.max_rects
        contact_max_rects = [ max_rect[0] for max_rect in max_rects if frontier in max_rect[0] ]
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
        margin_limit = self.min_size
        margined_space_forward_limit = space_forward_limit - margin_limit

        # Get the forward expansion limit according to maximum rectangles
        maximum_forward_limit = max([ rect.get_size()[forward] for rect in contact_max_rects ])
        margined_maximum_forward_limit = maximum_forward_limit - margin_limit

        # Now create a function to make a rectangle by pulling a segment (which may change)
        # The forward length of this rectangle will depend on the available area and limits
        # This segment may not be the original frontier, but a variation of it
        # The pull may happen using 2 different protocols, which are tried in the following order
        # - Protocol 1 (greedy): It tries to pull as much as it can, respecting the required area
        # - Protocol 2 (moderate): It tries to pull the minimum possible according to maximum rects
        def pull_segment (pulled_segment : Segment, protocol : int = 1, already_tried_length : Optional[number] = None) -> bool:
            pull_length = required_area / pulled_segment.length
            # First, try to pull using the maximum available space
            # This is risky since it may split the invaded room in two parts or make regions which do not respect the minimum size
            if protocol == 1:
                # If the length exceeds the maximum limit then stay in the maximum limit
                if pull_length > maximum_forward_limit:
                    pull_length = maximum_forward_limit
                # If the length is between the maximum limit and the margined maximum limit then stay at the margin
                elif pull_length > margined_maximum_forward_limit:
                    pull_length = margined_maximum_forward_limit
            # If the greedy try fails use the moderate try
            # Pull only using the closest row/column rectangle
            # This expansion is smaller but safe
            elif protocol == 2:
                # If at this point the length exceeds the space limit then stay in the space limit
                if pull_length > space_forward_limit:
                    pull_length = space_forward_limit
                # If the length is between the space limit and the margined space limit then stay at the margin
                elif pull_length > margined_space_forward_limit:
                    pull_length = margined_space_forward_limit
                # In case the length is bigger than the margined maximum we stay at the margin
                if pull_length > margined_maximum_forward_limit:
                    pull_length = margined_maximum_forward_limit
            # If the pull length was already tried, and thus failed, then stop here
            if pull_length == already_tried_length:
                return False
            # If the push length at this point is 0 or close to it then we can not push
            # WARNING: This usually happens because of forward limits, not because the area was not big enought
            # For this reason, trying to reduce the segment and push again will have no effect almost always
            if pull_length < minimum_resolution:
                #print('WARNING: The pull length is too small: ' + str(pull_length))
                return False
            # We must substract the new rect from this room and check everything is fine after
            if not self.pull_boundary_segment(pulled_segment, pull_length, check_parent_free_grid=False):
                # If the pull failed with the greedy protocol then try it again with the moderate protocol
                # Pass the already tried length so the pull is not repeated
                if protocol != 2:
                    return pull_segment(pulled_segment, 2, already_tried_length=pull_length)
                # Otherwise we surrender here
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
            if new_point.get_distance_to(other_point) < self.min_size:
                print('WARNING: The reduced frontier is not wide enought')
                return False
            reduced_frontier = Segment(new_point, other_point)

        return pull_segment(reduced_frontier)

    # Check all doors are respected given a new (truncated) grid
    def check_doors_required_spaces (self, truncated_grid : Grid) -> bool:
        # Get the doors which have been placed already
        stablished_doors = [ door for door in self.doors if door.point ]
        # Use this to visually check the contacts
        # truncated_bounday_segments = sum([ boundary.segments for boundary in truncated_grid.boundaries ], [])
        # truncated_bounday_segments = [ segment.get_colored_segment('blue') for segment in truncated_bounday_segments ]
        # door_segments = []
        # for door in stablished_doors:
        #     inside_segments = [ segment.get_colored_segment('green') for segment in door.get_required_space(inside=True).segments ]
        #     outside_segments = [ segment.get_colored_segment('red') for segment in door.get_required_space(inside=False).segments ]
        #     door_segments += inside_segments + outside_segments
        # elements_to_display = truncated_bounday_segments + door_segments
        # self.update_display(extra=elements_to_display, title='Doors checking')
        # For each stablished door,
        for door in stablished_doors:
            # Check the inside space to be fully inside
            inside_rect = door.get_required_space(inside=True)
            # This may happen when the door is not even in the boundary
            if not inside_rect:
                return False
            inside_space = Grid([inside_rect])
            uncovered_inside_space = inside_space.get_substract_grid(truncated_grid)
            if uncovered_inside_space:
                return False
            # Check the outside space to be fully outside
            outside_rect = door.get_required_space(inside=False)
             # This may happen when the door is not even in the boundary
            if not outside_rect:
                return False
            outside_space = Grid([outside_rect])
            covered_outside_space = outside_space.get_overlap_grid(truncated_grid)
            if covered_outside_space:
                return False
        return True

    # Check all doors are respected in the inside or outside region given a new (invader) grid
    def check_doors_side (self, invaded_grid : Grid, inside : bool) -> bool:
        # Get the doors which have been placed already
        stablished_doors = [ door for door in self.doors if door.point ]
        # If there are not doors yet then there is no problem at all
        if len(stablished_doors) == 0:
            return True
        # For each stablished door,
        for door in stablished_doors:
            # Check the outside space to be fully outside
            rect = door.get_required_space(inside=inside)
            # This may happen when the door is not even in the boundary
            if not rect:
                return False
            space = Grid([rect])
            covered_space = space.get_overlap_grid(invaded_grid)
            if covered_space:
                return False
        return True

    # Relocate doors in a new (truncated) grid
    # Return True if the relocation succed or False if it failed
    def relocate_doors (self, truncated_boundary : Boundary) -> bool:
        # Make a backup in case any door fails to relocate
        door_backups = { door: door.make_backup() for door in self.doors }
        # Get the truncated grid polygon
        truncated_polygon = truncated_boundary.exterior_polygon
        # Relocate each door
        succeed = True
        for door in self.doors:
            if not door.relocate(truncated_polygon):
                succeed = False
                break
        # Restore the backups in case something went wrong
        if not succeed:
            for door, backup in door_backups.items():
                door.restore_backup(backup)
        return succeed

    # Remove part of the room space and check everything is fine after
    # Use the force argument to skip all checkings and simply truncate the boundary
    # Use the easy argument to simply fail when minimum size is not respected in the truncated grid
    # Otherwise it tries to compenaste the grid before surrending, but this is not desired sometimes (e.g. when reducing corners)
    # Use the skip_update_display to avoid this truncate to generate a display frame
    # This is useful when truncating several rooms at the same time, so we do not have to recalculate the parent free grid every time
    def truncate (self,
        region : Grid,
        force : bool = False,
        easy : bool = False,
        check_parent_free_grid : bool = False,
        skip_update_display : bool = False,
        verbose : bool = False
    ) -> bool:
        if verbose: print(f'Truncating {self.name}')
        grid = self.grid
        # In case the room has not grid there is no problem at all in the truncation
        if not grid:
            return True
        # Save the removed region in a new variable wich may be further mutated
        removed_region = region
        # Calculate the boundary of this room after substracting the invaded region
        truncated_grid = grid - removed_region
        # If the grid has not been truncated then we have nothing to check
        if truncated_grid == grid:
            if verbose: print(' No truncation, no problem')
            return True
        # If the grid has been fully consumed then go back
        if not truncated_grid:
            if verbose: print(' Failed to truncate: The room would be fully consumed')
            return False
        # In case the truncate was forced remove all regions in the truncated grid which do not respect the minimum size
        if force:
            # Find regions in the truncated grid respecting the minimum size
            safe_grid = truncated_grid.keep_minimum(self.min_size)
            # The grid which is removed from the truncated grid by not respecting the minimum must be flaged as discarded
            lost_grid = truncated_grid - safe_grid
            # LORE: At this point the lost region was manually added to the discarded grid
            removed_region += lost_grid
            # Now update the truncated grid
            truncated_grid = safe_grid
        # If the grid has been fully consumed then give up
        if not truncated_grid:
            raise RuntimeError(f'We fully removed a whole room ({self.name}) by force-truncating')
        truncated_boundaries = truncated_grid.find_boundaries()
        # In case the room has been splitted in 2 parts as a result of the invasion we go back
        if not force and len(truncated_boundaries) > 1:
            if verbose: print(' Failed to truncate: The room would be splitted')
            return False
        truncated_boundary = truncated_boundaries[0]
        # Check doors to be respected
        if not force and not self.check_doors_required_spaces(truncated_grid):
            if verbose: print(' Door conflict when truncating: Doors will be relocated')
            if not self.relocate_doors(truncated_boundary):
                if verbose: print(' Failed to truncate: The room would not respect door spaces')
                return False
        # DANI: Es posible que se coma toda la habitación??
        if not force and len(truncated_boundaries) == 0:
            if verbose: print(' Failed to truncate: The room has been fully consumed')
            return False
        # Start saving a backup of self and children boundaries
        backup = {}
        # Now we must substract the removed region from every child with an already existing boundary
        # Note that it is not usual truncating a room with already set children, but we try to support this
        for child in self.children:
            # If th child is not yet set then skip it
            if not child.grid:
                continue
            # Save a backup of the current child in case we have to recover its boundary later
            child_grid_backup = child.grid
            # Truncate the child
            if not child.truncate(removed_region, force=force, check_parent_free_grid=check_parent_free_grid):
                # If the truncate process failed then restore backups and return True
                self.restore_grid_backup(backup)
                return False
            # Now add the child grid to the backup
            # Note that this is not done before since in case of failure the current children is backuped already
            backup[child] = child_grid_backup
        # Save a backup of the current grid
        backup[self] = self.grid
        # Set the truncated grid as the current grid
        self.set_grid(truncated_grid, skip_update_display=skip_update_display)
        # Check the truncated grid to still respecting the minimum size
        if not self.grid.check_minimum(self.min_size):
            # Surrender when the minimum check fails if the easy argument was passed
            if easy:
                if verbose: print(' Minimum size failed and grid compensation will not be tried -> Restoring backup')
                self.restore_grid_backup(backup)
                return False
            if verbose: print(' Minimum size conflict when truncating: Grid will be compensated')
            # If not then try to expand the truncated room accordingly to the truncated region to make it respect the minimum
            compensation_grid = self.grid.get_compensation_grid(region, self.min_size)
            if not self.expand_grid(compensation_grid, check_parent_free_grid=check_parent_free_grid):
                if verbose: print(' Failed to truncate: The room would not respect minimum size and it can not be compensated -> Restoring backup')
                self.restore_grid_backup(backup)
                return False
        # Check the new parent free grid to be respecting the minimum size in case it is requested
        # It may happen that we pull slightly a frontier which was in contact to a parent/brother frontier and create a small space
        # Note that parent grid is not checked in situtations where the truncated region is to be claimed by other room rigth away
        # DANI: Esto tiene un problema y es resolver un pull cuando ya está toda el area consumida -> se hace eterno
        # DANI: No se me ocurre solución sencilla así que de momento lo quito y ya
        if check_parent_free_grid:
            if not self.parent.free_grid.check_minimum(self.parent_free_limit):
                if verbose: print(f' Minimum size ({self.parent_free_limit}) is not respected in parent free space -> Restoring backup')
                # add_frame(self.parent.free_grid, title='Debug -> ' + self.parent.name)
                self.restore_grid_backup(backup)
                return False
        # If we made it this far then the truncation succeeded
        if verbose: print(f'Truncating {self.name}, succeeded')
        return True

    # Fuse a grid to self room grid
    # Check the new grid can be joined to the current grid as a single room
    # Run additional checkings to make sure other spaces are respected
    # The 'check_parent_free_grid' makes sure the new parent grid respects the minimum size
    # The 'compensate_invaded' flag runs a fitting for all brother rooms after truncation
    # The 'is_loaned' flag runs a fitting for self room after brother rooms truncation 
    # The 'conformist' flag is forwarded to the fitting functions
    # The 'easy' flag is forwarded to the truncate function
    def expand_grid (self,
        expansion_grid : Grid,
        check_parent_free_grid : bool = True,
        check_overlaps : bool = True,
        compensate_invaded : bool = True,
        is_loaned : bool = False,
        behaviour : str = 'exigent',
        easy : bool = False,
        verbose : bool = False
    ) -> bool:
        if verbose: print(f'Expanding room {self.name} at {expansion_grid}')
        new_grid = None
        # If the room has not a grid already then set the expansion grid as the current grid
        if self.grid == None:
            new_grid = expansion_grid
        # Otherwise create a new grid by mergin the current and the expansion grid
        else:
            new_grid = self.grid.get_merge_grid(expansion_grid, check_overlaps=check_overlaps)
        # In case there is a minimum size restriction check that the colliding segment is as long as required
        if not new_grid.check_minimum(self.min_size):
            if verbose: print(f' Grid expansion failed: Minimum size of the expanded grid would be not respected')
            return False
        # Make a backup of the current grid
        backup = { self: self.grid }
        # Set the new grid as the current grid
        self.grid = new_grid
        # In case this room has a parent room a few extra step must be done
        if self.parent:
            # Check the minimum size in the parent free grid to be respected
            if check_parent_free_grid and not self.parent.free_grid.check_minimum(self.parent_free_limit):
                if verbose: print('  Grid expansion failed: Parent minimum size would be not respected -> Restoring backup')
                # add_frame(self.parent.free_grid, 'Debug')
                self.restore_grid_backup(backup)
                return False
            # Check if we are exceeding parent boundaries
            parent_outer_grid = expansion_grid - self.parent.grid
            if parent_outer_grid:
                # If parent boundaries are flexible then we must expand the parent grid as well
                if self.parent._child_adaptable_boundary:
                    # If parent grid can not be expanded then go back
                    if not self.parent.expand_grid(expansion_grid):
                        if verbose: print(f' Grid expansion failed: Parent grid failed to expand -> Restoring backup')
                        self.restore_grid_backup(backup)
                        return False
                # If parent boundaries are rigid then we can not expand
                else:
                    if verbose: print(f' Grid expansion failed: Grid was expanded outer the parent grid -> Restoring backup')
                    self.restore_grid_backup(backup)
                    return False
            # Check if we are overlapping other spaces
            if check_overlaps:
                # Check if we are expanding over the corridor and, if so, abort the expansion
                if self.parent.corridor_grid and self.parent.corridor_grid.get_overlap_grid(expansion_grid):
                    if verbose: print(f' Grid expansion failed: Grid was expanded over the corridor -> Restoring backup')
                    self.restore_grid_backup(backup)
                    return False
                # If we have a adaptable parent then check we did not expand over any parent door outside space
                if self.parent._child_adaptable_boundary and not self.parent.check_doors_side(new_grid, inside=False):
                    # DANI: Tal vez se podría provar a recolocar la door, pero no quiero perder más tiempo con lo child adaptable
                    if verbose: print('  Grid expansion failed: Parent doors would be not respected')
                    return False
                # Now truncate all brothers
                # Note that doors which are not respected after the expansion are relocated here if possible
                brother_rooms = self.get_brother_rooms()
                for brother_room in brother_rooms:
                    if verbose: print(f'  Truncating {brother_room.name}')
                    # Make a backup grid of the current brother
                    # We do not save it yet to the bacup object since the grid would be backed up automatically if failed to truncate
                    brother_grid_backup = brother_room.grid
                    # Try to truncate the brother room
                    if not brother_room.truncate(expansion_grid, easy=easy):
                        if verbose: print(f' Grid expansion failed: Grid was expanded over a brother room which failed to truncate -> Restoring backup')
                        self.restore_grid_backup(backup)
                        return False
                    # Now save the previous backup
                    backup[brother_room] = brother_grid_backup
                # If this is a loaned push then we must retrieve the area debt at this point
                # Otherwise, if we invade other rooms it may happen that there is not free space enought for them to expand after
                if is_loaned and self.get_required_area(behaviour=behaviour) < 0:
                    # Set the new segments to respect during this fitting
                    # Note that if we pull these new segments then we woul be the same as before
                    new_segments = []
                    for boundary in expansion_grid.boundaries:
                        new_segments += boundary.segments
                    # Now run the fitting
                    if not self.fit_to_required_area(restricted_segments=new_segments, behaviour=behaviour):
                        if verbose: print(f' Grid expansion failed: Failed to fit self room after loaned expansion -> Restoring backup')
                        self.restore_grid_backup(backup)
                        return False
                # Now if we must compensate then fit truncated brothers
                if compensate_invaded:
                    for brother_room in brother_rooms:
                        # If the child is not yet set then skip it
                        if not brother_room.grid or brother_room.rigid:
                            continue
                        if verbose: print(f'  Compensating {brother_room.name} area')
                        # Try to fit the brother room
                        if not brother_room.fit_to_required_area(behaviour=behaviour):
                            if verbose: print(f' Grid expansion failed: Grid was expanded over a brother room which failed to comepnsate area after truncation -> Restoring backup')
                            self.restore_grid_backup(backup)
                            return False
        # At this point the expand succeed
        if verbose: print(f' Grid expansion succeed')
        return True

    # Get all overlapped segments between the current room and other
    def get_frontiers_with_room (self, other : 'Room') -> list:
        self_segments = self.boundary.segments
        other_segments = other.boundary.segments
        overlap_segments = []
        for self_segment in self_segments:
            for other_segment in other_segments:
                overlap_segment = self_segment.get_overlap_segment(other_segment)
                if overlap_segment:
                    # Save the room this segment belongs to inside the segment object
                    overlap_segment.rooms = [other]
                    overlap_segments.append(overlap_segment)
        return overlap_segments

    # Get all overlapped segments between the current room and a boundary
    def get_frontiers_with_boundary (self, boundary : 'Boundary') -> list:
        self_segments = self.boundary.exterior_polygon.segments
        other_segments = boundary.exterior_polygon.segments
        overlap_segments = []
        for self_segment in self_segments:
            for other_segment in other_segments:
                overlap_segment = self_segment.get_overlap_segment(other_segment)
                if overlap_segment:
                    overlap_segments.append(overlap_segment)
        return overlap_segments

    # Get all self frontiers segments
    # They come inside a tuple separated by free frontiers, conflict frontiers and forbidden frontiers respectively
    # i.e. free parent space frontiers, brother room frontiers and parent (if parent boundary) or rigid frontiers respectively
    def get_frontiers (self) -> tuple:
        # Now set which frontiers are free to expansion, which ones are in conflict and which ones are forbbiden
        forbidden_frontiers = []
        conflict_frontiers = []
        free_frontiers = []
        # Get the parent room
        parent_room = self.parent
        parent_frontiers = self.get_frontiers_with_room(parent_room)
        # In case we have a parent which is child adaptable, frontiers with its boundary are allowed to expansion
        # However, if there are free spaces inside the boundary then avoid expanding the perimeter until free spaces are covered
        if parent_room._child_adaptable_boundary and not parent_room.free_grid:
            free_frontiers += parent_frontiers
        # In case we have a parent which is not child adaptable, frontiers with its boundary are totally forbidden to expansion
        else:
            forbidden_frontiers += parent_frontiers
        # Other rooms inside the same parent may be displaced if there is no free space available
        brother_rooms = [ room for room in parent_room.children if room is not self and room.boundary ]
        for room in brother_rooms:
            # This function assign the frontier.rooms already
            brother_frontiers = self.get_frontiers_with_room(room)
            if not brother_frontiers:
                continue
            # If the room is flagged as 'rigid' then we add frontiers to the forbidden list
            if room.rigid:
                forbidden_frontiers += brother_frontiers
            else:
                conflict_frontiers += brother_frontiers
        # Find out which frontiers are connected to the parent corridor
        # These frontiers are forbidden
        if parent_room.corridor_grid:
            corridor_frontiers = []
            for boundary in parent_room.corridor_grid.boundaries:
                corridor_frontiers += self.get_frontiers_with_boundary(boundary)
            forbidden_frontiers += corridor_frontiers
        # The prefered limits to expand are those connected to free space inside the current parent
        # First of all get all already found frontier segments
        already_covered_frontiers = forbidden_frontiers + conflict_frontiers + free_frontiers
        # Now susbstract all previous frontiers from the rest of the external polygon to find free frontiers
        remaining_frontiers = []
        for segment in self.boundary.exterior_polygon.segments:
            remaining_frontiers += segment.substract_segments(already_covered_frontiers)
        free_frontiers += remaining_frontiers
        # Set the room all free segment belong to as None
        for frontier in free_frontiers:
            frontier.rooms = [parent_room]
        # Use this to display final frontiers
        # free_lines = [segment.get_colored_segment('green') for segment in free_frontiers]
        # conflict_lines = [segment.get_colored_segment('yellow') for segment in conflict_frontiers]
        # forbidden_lines = [segment.get_colored_segment('red') for segment in forbidden_frontiers]
        # self.update_display(title='Frontiers display', extra = free_lines + conflict_lines + forbidden_lines)
        # Return the calssified frontiers
        return free_frontiers, conflict_frontiers, forbidden_frontiers

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

    # Make backup of self and children room grids
    # A backup is a dict where keys are rooms and values are grids
    def make_grid_backup (self, children_only : bool = False) -> dict:
        backup = {}
        if not children_only:
            backup[self] = self.grid
        for child in self.children:
            backup[child] = child.grid
        return backup

    # Restore backup of room grids
    # A backup is a dict where keys are rooms and values are grids
    def restore_grid_backup (self, backup : dict):
        for room, grid_backup in backup.items():
            room.set_grid(grid_backup, skip_update_display=True)
        self.update_display(title='Restored grid backup')

    # Add a new frame in the display with the current segments of this room and its children
    # Also an 'extra' argument may be passed with extra segments to be represented
    def update_display (self, extra : list = [], title : Optional[str] = None):
        global display_options
        if not display_options['enabled']:
            return
        if display_options['frames_limit']:
            display_options['frames_limit'] -= 1
            if display_options['frames_limit'] == 0:
                raise SystemExit('Reached frames limit')
        # Find the root room
        root = self.get_root_room()
        # Get all children rooms recursively
        # This includes self room
        rooms = root.get_rooms_recuersive()
        # Display all current rooms together
        elements_to_display = [ *rooms, *extra ]
        add_frame(elements_to_display, title)

    # Solve room distributions
    # The display flag may be passed in order to generate a dynamic graph to display the solving process
    def solve (self) -> bool:
        return self.solve_children(recursive=True)

# The element which connects diferent floors of a building
class Stairs:
    def __init__ (self,
        # Set vairables to force the stairs location
        polygon : Optional[Polygon] = None,
        lower_door : Optional[Door] = None,
        upper_door : Optional[Door] = None,
        # Set variables to randomly generate stairs
        # Set the slope of the stairs
        # Default: 45°
        # WARNING: Note that length and slope are dependent and thus you can not pass both arguments
        slope : Optional[number] = None,
        # Set the stairs length
        # Deafult: length enough to match the slope according to the building height
        # WARNING: Note that length and slope are dependent and thus you can not pass both arguments
        length : Optional[number] = None,
        # Set the stairs width
        # Default: corridor size
        width : Optional[number] = None,
        # Note that height is not an argument, since we depend on the building height
        # In case stairs have to be generated randomly there are 4 possible configurations
        # 0 - Vertical: Just a hole for vertical stairs or elevators (slope has no effect and length is equal to width)
        # 1 - One line: Regular straight stairs. The room is expected to have a rectangular polygon
        # 2 - Two lines: Stairs with a corner at some point. The room is expected to have a 'L' polygon
        # 3 - Three lines: Stairs with two corners. The room is expected to have a squared polygon
        # Default: 1
        configuration : Optional[int] = None,
        # Set a flag to define if corners are flat or 'staired'
        # DANI: aquí haría falta una imágen: flat son las de la escalera de casa y staired las de casa de Juan Luís y Mariajosé
        # Note that this value only makes sense for stairs with more than 1 line
        # Note that this argument is valid for both forced and random stairs
        # (When forced) Note that when corners are flat they don't count for the stairs length and thus the slope required is higher
        # (When random) Note that when corners are flat they don't count for the stairs length and thus the space required is bigger
        flat_corners : bool = True,
        # Set the parent building
        parent_building : Optional['Building'] = None,
        # Set the parent building floors connected by the stairs
        # Note that each stairs connect only one floor to another
        # Stairs connecting multiple floors are actually several starts stacked one over the other
        lower_floor_number : Optional[int] = None,
        upper_floor_number : Optional[int] = None,
        # Set if stairs are to be stacked one upon the other
        # This means that the upper room of a lower floor stairs will be identical to the lower room of an upper floor stairs with the same configuration, when possible
        # Note that rooms will have the same polygon and will totally overlap but they will be not the same room since they will have different doors
        # This is meant so save space and it is a very common feature in realistic buildings
        stacking_stairs : bool = True,
        
    ):
        # Save the initiation arguments
        self.polygon = polygon
        self._lower_door = lower_door
        self._upper_door = upper_door
        self._slope = slope
        self._length = length
        self._width = width
        self.configuration = configuration
        self.flat_corners = flat_corners
        self.parent_building = parent_building
        self.lower_floor_number = lower_floor_number
        self.upper_floor_number = upper_floor_number
        self.stacking_stairs = stacking_stairs
        # Check floor numbers to be correct
        if lower_floor_number >= upper_floor_number:
            raise InputError('Lower floor number must be lower than upper floor number in a stair')
        # Set internal variables
        self._lower_room = None
        self._upper_room = None
        # Forced scenario:
        if polygon:
            # Check there are no redunadancies and, if so, warn the user
            # If the polygon has been forced then all arguments for the random setup of the stairs make no sense
            if slope:
                print('WARNING: Redundant input "slope" for stairs with already forced polygon')
            if width:
                print('WARNING: Redundant input "width" for stairs with already forced polygon')
            if length:
                print('WARNING: Redundant input "length" for stairs with already forced polygon')
            if configuration:
                print('WARNING: Redundant input "configuration" for stairs with already forced polygon')
            # Check doors to be over the polygon
            if self._lower_door and self._lower_door.point not in polygon:
                raise InputError('Lower door is not over the polygon')
            if self._upper_door and self._upper_door.point not in polygon:
                raise InputError('Upper door is not over the polygon')
        # Random scenario:
        else:
            # Doors can not be forced if the polygon is not forced
            if self._lower_door or self._upper_door:
                raise InputError('You can not force stairs doors if the polygon is not forced as well')
            # Set the defaults
            # Slope and length cannot be both passed
            if slope and length:
                raise InputError('Length and slope are dependent and thus you can not pass both arguments')
            # If one of the two is defined (slope or length) then to know the other we need the height
            # We may not have the height yet so we must wait
            # If any of the two parameters is passed (slope nor length) then we set the default value for the slop
            if not slope and not length:
                self._slope = 45
            # Set the default configuration
            if configuration == None:
                self.configuration = 1

    # Stairs will always have two floors: one upper and one lower

    # Get the upper floor
    def get_upper_floor (self) -> Room:
        if not self.parent_building:
            raise RuntimeError('Trying to get upper floor of a stair with no parent building')
        if self.upper_floor_number == None:
            raise RuntimeError('Trying to get upper floor of a stair with no upper floor number')
        return self.parent_building.floors[self.upper_floor_number]

    # The upper floor
    upper_floor = property(get_upper_floor, None, None, "The upper floor (read only)")

    # Get the lower floor
    def get_lower_floor (self) -> Room:
        if not self.parent_building:
            raise RuntimeError('Trying to get lower floor of a stair with no parent building')
        if self.lower_floor_number == None:
            raise RuntimeError('Trying to get lower floor of a stair with no lower floor number')
        return self.parent_building.floors[self.lower_floor_number]

    # The lower floor
    lower_floor = property(get_lower_floor, None, None, "The lower floor (read only)")

    # Stairs will always have two rooms: one lower and one upper
    # These rooms will always overlap in the boundaries
    # These rooms set the place for the actual stairs, which take place in both floors
    # However, the stairs between 2 fllors are always defined in the lower floor stairs
    # Note that stairs may stack, thus sharing the same room along different couples of floors

    # Generate both the upper and lower rooms
    def setup_rooms (self):
        # If polygon is forced:
        if self.polygon:
            # Set the lower door
            if not self._lower_door:
                self._lower_door = Door(rigid=True)
            # Set the upper door
            if not self._upper_door:
                self._upper_door = Door(rigid=True)
        # If polygon is random:
        else:
            # This function sets self polygon and doors
            self.set_place()
        # Set the rooms
        boundary = Boundary(self.polygon)
        self._lower_room = Room(boundary=boundary, doors=[ self._lower_door ], parent=self.lower_floor, rigid=True, min_size=self.width, name='Lower')
        self._upper_room = Room(boundary=boundary, doors=[ self._upper_door ], parent=self.upper_floor, rigid=True, min_size=self.width, name='Upper')

    # Get the lower door
    def get_lower_door (self) -> Door:
        # Return the internal value if it exists already
        if self._lower_door:
            return self._lower_door
        # Otherwise we must set the lower door
        self.setup_rooms()
        return self._lower_door

    # The lower door
    lower_door = property(get_lower_door, None, None, "The lower door (read only)")

    # Get the lower room
    def get_lower_room (self) -> Room:
        # Return the internal value if it exists already
        if self._lower_room:
            return self._lower_room
        # Otherwise we must set the lower room
        self.setup_rooms()
        return self._lower_room

    # The lower room
    lower_room = property(get_lower_room, None, None, "The lower room (read only)")

    # Get the upper door
    def get_upper_door (self) -> Door:
        # Return the internal value if it exists already
        if self._upper_door:
            return self._upper_door
        # Otherwise we must set the upper door
        self.setup_rooms()
        return self._upper_door

    # The upper door
    upper_door = property(get_upper_door, None, None, "The upper door (read only)")

    # Get the upper room
    def get_upper_room (self) -> Room:
        # Return the internal value if it exists already
        if self._upper_room:
            return self._upper_room
        # Otherwise we must set the upper room
        self.setup_rooms()
        return self._upper_room

    # The upper room
    upper_room = property(get_upper_room, None, None, "The upper room (read only)")

    # Get the height
    def get_height (self) -> number:
        if not self.lower_floor:
            raise ValueError('Trying to get height of a stair with no lower floor')
        return self.lower_floor.height

    # The height
    height = property(get_height, None, None, "The height (read only)")

    # Get the width
    def get_width (self) -> number:
        # Return the stored value, if any
        # This means the width has been forced from the arguments
        if self._width != None:
            return self._width
        # Otherwise we must use the lower floor corridor size
        if not self.lower_floor:
            raise ValueError('You are requesting the width of stairs with no lower floor when this values was not passed')
        width = self.lower_floor.corridor_size
        #print('WIDTH: ' + str(width))
        if width == None:
            raise ValueError('Lower floor has no corridor size')
        return width

    # The width
    width = property(get_width, None, None, "The width (read only)")

    # Get the length
    def get_length (self) -> number:
        # Return the stored value, if any
        # This means the length has been forced from the arguments or previously calculated
        if self._length != None:
            return self._length
        # Otherwise we must calculate the length
        # Note that the height is required for this calculation and thus the lower floor must be set already
        self._length = resolute( self.height / tan(radians(self.slope)) )
        return self._length

    # Set the length
    # Modify the slope to make it coherent
    def set_length (self, new_length : number):
        self._length = new_length
        # Note that the height is required for this calculation and thus the lower floor must be set already
        new_slope = degrees( atan( self.height / new_length ) )
        self._slope = resolute(new_slope)

    # The length
    length = property(get_length, set_length, None, "The length")

    # Get the slope
    def get_slope (self) -> number:
        # Return the stored value, if any
        # This means the slope has been forced from the arguments, set by default, or previously calculated
        if self._slope != None:
            return self._slope
        # Otherwise we must calculate the slope
        # Note that this happens when length is passed as argument
        # Note that the height is required for this calculation and thus the lower floor must be set already
        self._slope = degrees( atan( self.height / self.length ) )
        return self._slope

    # Set the slope
    # Modify the length to make it coherent
    def set_slope (self, new_slope : number):
        self._slope = new_slope
        # Note that the height is required for this calculation and thus the lower floor must be set already
        new_length = self.height / tan(radians(new_slope))
        self._length = new_length

    # The slope
    slope = property(get_slope, set_slope, None, "The slope")

    # Set the polygon, lower door and upper door according to the stair parameters and floor avilable space
    # Here is decided both shape and position of the stair rooms and their doors also considering extra space for the corridor
    # Note that the extra space may not be the same in upper and lower floors
    # Note that this additional space is unserstood as simply an extension of the stair length which keeps the stair width
    def set_place (self):
        # Set the lower and upper floor grids
        # If there is no upper floor then asume the available space will be the lower floor
        lower_floor_grid = self.lower_floor.grid
        upper_floor_grid = self.upper_floor.grid
        if not upper_floor_grid:
            upper_floor_grid = lower_floor_grid
        # Set the additional space length required
        # Note that the additional space width will be the stairs width to avoid difficulties, but length may change
        required_free_space_length = self.lower_floor.door_args['length']
        # Rectangle with both doors in the same place (e.g. elevators)
        if self.configuration == 0:
            raise SystemExit('Configuration 0 is not yet programmed :(')
        # Rectangle with a door on each extreme (e.g. regular linear stairs)
        elif self.configuration == 1:
            # Stair may be oriented in any direction and it may tell the difference between fitting or not
            for direction in [ UP, RIGHT, DOWN, LEFT ]:
                is_vertical = direction.is_vertical()
                x_size = self.width if is_vertical else self.length + required_free_space_length
                y_size = self.length + required_free_space_length if is_vertical else self.width
                x_offset = 0 if is_vertical else required_free_space_length * direction.x # Direction x is to be -1 or 1
                y_offset = required_free_space_length * direction.y if is_vertical else 0 # Direction y is to be -1 or 1
                # Get both available spaces, lower and upper, taking in count each other
                lower_available_space, upper_available_space = get_tandem_fitting_grid(
                    grid_1 = lower_floor_grid, x_fit_size_1 = x_size, y_fit_size_1 = y_size,
                    grid_2 = upper_floor_grid, x_fit_size_2 = x_size, y_fit_size_2 = y_size,
                    x_position_offset = x_offset, y_position_offset = y_offset )
                # Now we must define a random spot in the available lower space and set its equivalent in the upper room
                lower_spot = next(lower_available_space.generate_fitting_spots(x_size, y_size))
                upper_spot = lower_spot.get_offset_rect(x_position_offset = x_offset, y_position_offset = y_offset)
                # Now we tell apart the free space from the actual room space
                x_size_offset = 0 if is_vertical else -self.length
                y_size_offset = -self.length if is_vertical else 0
                lower_x_position_offset = self.length if direction == LEFT else 0
                lower_y_position_offset = self.length if direction == DOWN else 0
                lower_free_space = lower_spot.get_offset_rect( lower_x_position_offset, lower_y_position_offset, x_size_offset, y_size_offset )
                lower_free_grid = Grid([lower_free_space])
                upper_x_position_offset = self.length if direction == RIGHT else 0
                upper_y_position_offset = self.length if direction == UP else 0
                upper_free_space = upper_spot.get_offset_rect( upper_x_position_offset, upper_y_position_offset, x_size_offset, y_size_offset )
                upper_free_grid = Grid([upper_free_space])
                # Check boths spots match after substracting free spaces
                lower_grid = Grid([lower_spot]) - lower_free_grid
                upper_grid = Grid([upper_spot]) - upper_free_grid
                if lower_grid != upper_grid:
                    raise Exception('Lower and upper grids must be identical')
                if len(lower_grid.boundaries) != 1:
                    raise Exception('Grid must have only one boundary')
                # Set the polygon
                self.polygon = lower_grid.boundaries[0].exterior_polygon
                # Set the door positions
                if is_vertical:
                    lower_door_point_x_position = upper_door_point_x_position = lower_free_space.x_min + lower_free_space.x_size / 2
                    lower_door_point_y_position = lower_free_space.y_min if direction == DOWN else lower_free_space.y_max
                    upper_door_point_y_position = upper_free_space.y_max if direction == DOWN else upper_free_space.y_min
                else:
                    lower_door_point_x_position = lower_free_space.x_min if direction == LEFT else lower_free_space.x_max
                    upper_door_point_x_position = upper_free_space.x_max if direction == LEFT else upper_free_space.x_min
                    lower_door_point_y_position = upper_door_point_y_position = lower_free_space.y_min + lower_free_space.y_size / 2
                lower_door_point = Point(lower_door_point_x_position, lower_door_point_y_position)
                self._lower_door = Door(point = lower_door_point, rigid=True)
                upper_door_point = Point(upper_door_point_x_position, upper_door_point_y_position)
                self._upper_door = Door(point = upper_door_point, rigid=True)
        # If the configuration is not recognized then raise an input error
        else:
            raise InputError('Stairs configuration ' + str(self.configuration) + ' is not defined')


# The building which may contain several floors
class Building:
    def __init__ (self,
        # A dict containing all floors in the building
        # Each floor is a room which may contains several rooms
        # Keys are the floor number. The 0 is the base. Negative numbers stand for the basements
        floors : dict,
        # Set the default room arguments
        # These values are applied to all rooms unless other values are assigned explicitly
        room_args : Optional[dict] = None,
        # Set the stairs
        # Note that each floor may have several stairs
        # If not stairs are provided then the default is one stairs per floor
        stairs : Optional[ List['Stairs'] ] = None,
    ):
        # Set the floors
        self.floors = floors
        # Check input floors to make sense according to floor indices
        floor_indices = sorted(list(floors.keys()))
        self.floor_indices = floor_indices
        lowest_floor_number = min(floor_indices)
        self.lowest_floor_number = lowest_floor_number
        highest_floor_number = max(floor_indices)
        self.highest_floor_number = highest_floor_number
        # There must always be a floor 0
        if not 0 in floor_indices:
            raise InputError('A building must have a floor 0 (i.e. the first floor)')
        # Floor indices must not have gaps
        # i.e. if there is a floor 1 and a floor 3 then there must be a floor 2
        for index in range(lowest_floor_number +1, highest_floor_number):
            if index not in floor_indices:
                raise InputError('Missing floor ' + str(index))
        # Set the parent building in all floors
        for floor in self.floors.values():
            floor.parent_building = self
        # Set the stairs
        self.stairs = stairs
        # If stairs are set by the user then make a few checks
        if stairs:
            for stair in stairs:
                if stair.lower_floor_number >= highest_floor_number:
                    raise InputError('Stair lower floor number must be below the building highest floor number')
                if stair.upper_floor_number <= lowest_floor_number:
                    raise InputError('Stair upper floor number must be above the building lowest floor number')
        # If stairs are missing then set the default values
        # By default all floors are connected by one stairs
        else:
            self.stairs = []
            for floor_number_a, floor_number_b in pairwise(floor_indices):
                default_stairs = Stairs(parent_building=self, lower_floor_number=floor_number_a, upper_floor_number=floor_number_b)
                self.stairs.append(default_stairs)
        # Set the room args
        self.room_args = room_args if room_args else {}
        # Complete the missing values with some default values
        # If the height is not passed try to guess a reasonable height
        if room_args.get('height', None) == None:
            # DANI: Ya pensaré algo un poco más elaborado
            self.room_args['height'] = 30
        # If minimum size is not passed then guess a reasonable value as well
        if room_args.get('min_size', None) == None:
            # DANI: Lo mismo
            self.room_args['min_size'] = 10
        # Calcualte the overall min size which may be useful to set some default values
        floor_min_sizes = [ floor.get_min_size_recursive() for floor in self.floors.values() ]
        overall_min_size = min(floor_min_sizes)
        if not overall_min_size:
            raise InputError('Cannot guess the overall minimum size. Please set a minimum size somewhere')
        # If the corridor size is missing guess a resonable size from the overall minimum size
        # Note that making the corridor slightly thiner than the min size makes things easier specially in upper floors
        # This avoids the inherited corridor regions between rigid boundaries to be filled by a room
        # This is not a "problem", but depending on the configuration the result may be "not elegant"
        if room_args.get('corridor_size', None) == None:
            # DANI: No me, lo suyo es que mida lo mismo el pasillo
            self.room_args['corridor_size'] = overall_min_size
        # If the door args are missing guess resonable values from the corridor size
        if room_args.get('door_args', None) == None:
            # Make the margined width of all doors equal to the corridor size
            # Make the width of all doors the 80% of the margined width
            margined_width = room_args['corridor_size']
            width = margined_width * 0.8
            margin = margined_width * 0.1
            length = margined_width
            self.room_args['door_args'] = {
                'width': width,
                'margin': margin,
                'length' : length
            }

    # DANI: Al final esto no lo he implementado
    # Fix corridor sized regions in a non-base floor which are produced by actual corridors in previous floors
    # These regions are common and they appear between rigid rooms (e.g. stairs) and the floor perimeter
    # They are problematic since they are too thin to place rooms (even if possible the result is not realistic)
    # Note that if stairs are stackable and this is not the last roof then this space must become corridor inmediatelly
    # These regions are handled here to improve the result. For now they are just discarded.
    # In future implementations these regions may be recycled as:
    # - Additional corridor (useless but estetical)
    # And if they are in contact with the perimeter (which is very common) they may be also recycled as:
    # - Exterior balconies
    # - Little roofs (i.e. perimeter truncation)
    # *** Note that there may be no problem in reducing the perimeter since there may be no upper floors (also very commom)
    def _fix_inherited_ghost_corridor_regions (self, floor : Room):
        # Find these inherited ghost corridor regions
        # Get all free regions not respecting the parent free limit (i.e. the highest minimum size among its children)
        corridor_free_regions = floor.free_grid
        correct_regions = corridor_free_regions.keep_minimum(floor.parent_free_limit)
        wrong_regions = corridor_free_regions - correct_regions
        if not wrong_regions:
            return
        # Keep only those which respect the corridor size
        ghost_corridor_regions = Grid()
        for wrong_region in wrong_regions.find_connected_grids():
            if wrong_region.check_minimum(floor.corridor_size):
                ghost_corridor_regions += wrong_region
            

    # Solve all floors
    def solve (self):
        # Solve each floor starting by the floor 0 (the first floor), then solving the superior floors and finally the basements
        sorted_floor_indices = list(range(self.highest_floor_number +1)) + list(range(self.lowest_floor_number, 0))
        for floor_number in sorted_floor_indices:
            floor = self.floors[floor_number]
            upward_stairs = [ stair for stair in self.stairs if stair.lower_floor_number == floor_number ]
            downward_stairs = [ stair for stair in self.stairs if stair.upper_floor_number == floor_number ]
            # Get the lower floor, it may be useful to aset a few parameter of the current one
            lower_floor_number = floor_number - 1
            # Get the upper floor, it may be useful to aset a few parameter of the current one
            upper_floor_number = floor_number + 1
            # In case this floor has not a forced boundary,
            if not floor.input_boundary:
                # Basement floors
                if floor_number < 0:
                    # We set same boundary as the base
                    base_boundary = self.floors[0].boundary
                    floor.boundary = base_boundary
                    floor._child_adaptable_boundary = False
                # Upper floors
                elif floor_number > 0:
                    # We set same boundary as the lower floor
                    lower_floor_boundary = self.floors[lower_floor_number].boundary
                    floor.boundary = lower_floor_boundary
                    floor._child_adaptable_boundary = False
                # Base floor
                else:
                    # If we have an area range in the inputs then generate a random polygon
                    if floor.min_area != None and floor.max_area != None:
                        # Generate a random polygon and set it as the floor boundary
                        random_polygon = generate_random_polygon(
                            min_area=floor.min_area,
                            max_area=floor.max_area,
                            min_size=floor.min_size
                        )
                        floor.boundary = Boundary(random_polygon)
                    # Otherwise, we have a child adaptable boundary scenario
                    else:
                        raise Exception('This scenario is not yet supported')
                        # # We may need to set a bit of extra free space next to the stairs
                        # # This is for the perimeter in the floor above to have space for the corridor to reach the stairs door
                        # # Note that this has to be done now that the stair rooms are children of the floor
                        # # Otherwise the door is not able to find its arguments (width and margin) at this point since it has not root
                        # if stairs:
                        #     # Now calculate the space required by the door
                        #     for stair in stairs:
                        #         extra_space = Grid([ stair.upper_room.doors[0].get_required_space(inside=False) ])
                        #         floor.forced_grid += extra_space
            # In case this floor has not forced doors there will be not doors incase it is not the base
            if floor.input_doors == None and floor_number != 0:
                floor.doors = []
            # Set the stairs
            for stair in upward_stairs + downward_stairs:
                # At this moment we set both lower and upper room positions and shapes
                # If stairs already have a polygon then skip this part
                if not stair.polygon:
                    stair.set_place()
            # Set the lower and upper rooms to the floor children doors
            # Note that by declaring them they are parented automatically
            for stair in upward_stairs:
                stair.lower_room
            for stair in downward_stairs:
                stair.upper_room
            # Start the whole solving process
            floor.solve()

# -----------------------------------------------------