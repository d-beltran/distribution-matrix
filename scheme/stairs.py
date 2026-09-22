import random
from math import sqrt, inf, tan, atan, degrees, radians
from typing import List, Tuple, Dict, Union, Optional

from utils.auxiliar import *
from vectorial_base import *
from utils.display import add_frame

from scheme.door import Door
from scheme.room import Room

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
        self._lower_room = Room(boundary=boundary, doors=[ self._lower_door ], rigid=True, min_size=self.width, name='Lower staris')
        self.lower_floor.add_child(self._lower_room)
        self._upper_room = Room(boundary=boundary, doors=[ self._upper_door ], rigid=True, min_size=self.width, name='Upper staris')
        self.upper_floor.add_child(self._upper_room)

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
    def set_place (self, debug : bool = True) -> bool:
        # Set the lower and upper floor grids
        # If there is no upper floor then asume the available space will be the lower floor
        lower_floor_grid = self.lower_floor.grid
        upper_floor_grid = self.upper_floor.grid
        if not upper_floor_grid:
            upper_floor_grid = lower_floor_grid
        # Set the additional space length required
        # Note that the additional space width will be the stairs width to avoid difficulties, but length may change
        required_lower_space = max(self.lower_floor.door_args['length'], self.lower_floor.min_size)
        required_upper_space = max(self.upper_floor.door_args['length'], self.upper_floor.min_size)
        # Rectangle with both doors in the same place (e.g. elevators)
        if self.configuration == 0:
            raise SystemExit('Configuration 0 is not yet programmed :(')
        # Rectangle with a door on each extreme (e.g. regular linear stairs)
        elif self.configuration == 1:
            # Stair may be oriented in any direction and it may tell the difference between fitting or not
            for direction in [ UP, RIGHT, DOWN, LEFT ]:
                is_vertical = direction.is_vertical()
                lower_x_size = self.width if is_vertical else self.length + required_lower_space
                lower_y_size = self.length + required_lower_space if is_vertical else self.width
                upper_x_size = self.width if is_vertical else self.length + required_upper_space
                upper_y_size = self.length + required_upper_space if is_vertical else self.width
                x_offset = 0 if is_vertical else required_lower_space * direction.x # Direction x is to be -1 or 1
                y_offset = required_lower_space * direction.y if is_vertical else 0 # Direction y is to be -1 or 1
                # Get both available spaces, lower and upper, taking in count each other
                lower_available_space, upper_available_space = get_tandem_fitting_grid(
                    grid_1 = lower_floor_grid, x_fit_size_1 = lower_x_size, y_fit_size_1 = lower_y_size,
                    grid_2 = upper_floor_grid, x_fit_size_2 = upper_x_size, y_fit_size_2 = upper_y_size,
                    x_position_offset = x_offset, y_position_offset = y_offset )
                # Now we must define a random spot in the available lower space and set its equivalent in the upper room
                # DANI: Aquí la lógica es limitada y no resolverá espacios complejos
                # DANI: Lo suyo sería usar lower_available_space.generate_fitting_regions_with_margin()
                # DANI: Sin embargo esta función aún no soporta el "espacio libre ya asignado" de la escalera
                lower_spot = next(lower_available_space.generate_fitting_spots(lower_x_size, lower_y_size), None)
                if lower_spot is None: continue
                # Make sure the lower spot respects the minimum size
                if not lower_available_space.does_rect_fit(lower_spot, self.lower_floor.min_size): continue
                # Get the corresponding upper spot
                upper_spot = lower_spot.get_offset_rect(x_position_offset = x_offset, y_position_offset = y_offset)
                # Make sure the upper spot fits as well
                if not upper_available_space.does_rect_fit(upper_spot, self.upper_floor.min_size): continue
                # Show the current position of the spots if we are to debug
                if debug:
                    lower_available_space.color, upper_available_space.color = 'red', 'blue'
                    lower_spot.color, upper_spot.color = 'orange', 'purple'
                    elements_to_display = [ lower_available_space, upper_available_space, lower_spot, upper_spot ]
                    add_frame(elements_to_display, title='Available space to set the stairs')
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
                if lower_grid and not lower_grid.is_unified():
                    raise Exception('Grid must not be splitted, but unified')
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
                # If we made it this far then it means we successfully placed the stairs
                return True
            # If we could not find a single spot in any direction then we surrender
            return False
        # If the configuration is not recognized then raise an input error
        else:
            raise InputError(f'Stairs configuration {self.configuration} is not defined')
