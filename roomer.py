from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some predefined test polygons
from tests import *

# This is for windows to dont loop
if __name__ == '__main__':

    # Set if we want to display the solving process
    display = True

    # Represent current rooms tagged as display = True
    if display:
        setup_display()

    test = test_room_2
    #test = test_building_1

    test.solve(display=display)

    print('Done!')