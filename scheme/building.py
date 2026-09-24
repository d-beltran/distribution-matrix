import random
from math import sqrt, inf, tan, atan, degrees, radians
from typing import List, Tuple, Dict, Union, Optional

from utils.auxiliar import *
from vectorial_base import *
from utils.display import add_frame

from scheme.door import Door
from scheme.window import Window
from scheme.room import Room
from scheme.stairs import Stairs

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
        # Get the floor's parent free limit to be respected
        free_limit = floor.get_parent_free_limit()
        # Find these inherited ghost corridor regions
        # Get all free regions not respecting the parent free limit (i.e. the highest minimum size among its children)
        corridor_free_regions = floor.free_grid
        correct_regions = corridor_free_regions.keep_minimum(free_limit)
        wrong_regions = corridor_free_regions - correct_regions
        if not wrong_regions:
            return
        # Keep only those which respect the corridor size
        ghost_corridor_regions = Grid()
        for wrong_region in wrong_regions.find_connected_grids():
            if wrong_region.check_minimum(floor.corridor_size):
                ghost_corridor_regions += wrong_region
            

    # Solve all floors
    def solve (self) -> bool:
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
            # Set the facade windows now that the floor boundary is set, before the stairs are placed
            # Floors with a boundary inherited from the lower floor inherit its windows as well, so all floors share the same axes
            # Basements have no windows unless they are forced
            if floor.input_windows == None and not floor.input_boundary and floor_number != 0:
                if floor_number > 0:
                    lower_floor_windows = self.floors[lower_floor_number].windows
                    floor.windows = [ Window(point=window.point, width=window.width, margin=window.margin) for window in lower_floor_windows ]
                floor.set_facade_windows(fill=False)
            else:
                floor.set_facade_windows()
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
                # Now that the staris have been placed, update the view from the lower floor point of view
                stair.lower_room.update_display(title=f'Stairs have been placed, where room "{stair.lower_room.name}" is the lower floor')
            for stair in downward_stairs:
                stair.upper_room
            # Start the whole solving process
            if not floor.solve():
                return False
        return True