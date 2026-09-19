import random
from math import sqrt, inf, tan, atan, degrees, radians
from typing import List, Tuple, Dict, Union, Optional

from utils.auxiliar import *
from vectorial_base import *

# A door is a segment in a boundary
# When boundaries are transformed to walls with tickness, doors become holes in the wall
# Doors connect children rooms to the parent room corridor
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
        room : Optional['Room'] = None,
        # Set a name for the room
        # This is a representation parameters and it has no effect in the logic
        name : Optional[str] = None,
    ):  
        # Save input values as internal values
        # These values are usually None at this point
        # They are usually set further from the door room 'door_args' value
        self.name = name
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
            raise InputError('A point must be defined if the door is to be rigid')
        self.reverse = reverse
        # The room this door belongs to
        self.room = room

    def __repr__ (self):
        name = self.name if self.name else 'Unnamed'
        point = f'placed in {self._point}' if self._point else '(not placed)'
        width = f'with a width of {self._width}' if self._width else '(widthless)'
        margin = f'and with a margin of {self._margin}' if self._margin else '(marginless)'
        return f'<Door "{name}" {point} {width} {margin}>'

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

    # Make a copy of this door
    # Make sure we can mutate the copy without having any effect on the original
    def copy(self) -> 'Door':
        return Door(
            point = self.point,
            width = self.width,
            margin = self.margin,
            length = self.length,
            direction = self.direction,
            pivot = self.pivot,
            rigid = self.rigid,
            reverse = self.reverse,
            room = self.room,
        )
