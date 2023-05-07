from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some python libraries to trace and benchmark
from traceback import print_exc
from time import time

# Import some predefined test polygons
from tests import *

# This is for windows to dont loop
if __name__ == '__main__':

    # Set if we want to display the solving process
    display = True

    # Represent current rooms tagged as display = True
    if display:
        setup_display()

    #test = test_room_1
    #test = test_building_1
    test = test_building_2

    # Record the current time so we can then calculate how much time it took to run all the process
    start_time = time()

    # Start the solving process
    try:
        test.solve(display=display)
    except Exception as e:
        print_exc()

    # Calculate how much it took to run the whole process and output the result
    end_time = time()
    total_time = round((end_time - start_time) * 100) / 100
    print(' -- The process took ' + str(total_time) + ' seconds to run --')