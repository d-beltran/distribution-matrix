from scheme_display import setup_display
from vectorial_base import *
from scheme import *
from auxiliar import round_to_hundredths

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
        setup_display(frames_limit=10000)


    #test = test_room_1
    #test = test_building_1
    test = test_building_3

    # Record the current time so we can then calculate how much time it took to run all the process
    start_time = time()

    #test = generate_random_polygon()

    # Start the solving process
    try:
        if test.solve():
            print('Done :)')
        else:
            print('Failed :(')
    # An input error means the problem itself can not be solved as it was defined
    # Note that it is caught apart so it is clearly reported and not mistaken for a crash
    # Note that it inherits from SystemExit, so otherwise it would kill the script silently
    # WARNING: Its message would be printed by python through stderr, while the whole log goes to stdout
    # WARNING: Thus it would show up misplaced among thousands of lines and it would be missed
    except InputError as error:
        print(f'Input error: {error}')
    except Exception as e:
        print_exc()

    # Calculate how much it took to run the whole process and output the result
    end_time = time()
    total_time = round_to_hundredths(end_time - start_time)
    print(f' -- The process took {total_time} seconds to run --')