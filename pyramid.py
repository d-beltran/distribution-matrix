from multiprocessing import Manager, Pool, cpu_count
from multiprocessing.managers import BaseManager
from threading import Lock
from queue import Empty
from traceback import format_exception

from typing import List, Union, Optional

# Get the number of available CPUs
AVAILABLE_CPUS : int = cpu_count()
# Set the number of parallel process to launch within reason
N_PARALEL = AVAILABLE_CPUS - 2 if AVAILABLE_CPUS > 0 else 1
print(f'Allowing {N_PARALEL} paralelizable processes from {AVAILABLE_CPUS} available CPUs')

# In order to solve children configurations faster we will paralelize several processes looking for a solution.
# To do so, we will use the following strategy:
#
# Given a fixed order of rooms to be solved, think of a pyramid of possible room configurations.
# We can place the first room in different positions and this would be the first layer of the pyramid.
# Then, for every initial placement, we can place the following room in different positions and so forth.
# This is an iterative process with as many layers as rooms to solve and the possibilites may be a lot.
# The last layer of the pyramid is made of configurations including all the rooms, which may be tested.
# Testing a configuration means trying to fit the corridor. Note that a valid configuration may fail here.
# If we test a final configuration and it works then our job is done here.
# However there may be a scenario where there are few or no solutions.
# And the only way to know this is to explore the whole pyramid, so we must be prepared to do it.
# 
# In order to explore this pyramid of possible configurations we will use a bottom-top strategy.
# Every process will go straight to the last layer of the pyramid by placing one room after another.
# Then it will try if the configuration works and, if not, it will just change the last room and try again.
# This means every process will potentially find a solution in the first try
# In the other hand, the paralelization will throw new processes in a top-bottom strategy.
# Thus every independent process will explore a very different region of the pyramid.
# This is important since sometimes the starting configurations may be "doomed" by just luck.
# A single process will take a lot of time to check all possibilities before surrendering.
#
# We must coordinate these two opposite strategies which may eventually find each other:
# The independent processes may not find any solutions and climb up in the pyramid
# If they climb up enought they will reach the level where the "paralelizer" is working
# We must detect this scenario so the processes stop at this point to not waste the work
# To do so, we will keep track of every region of the pyramid which has been already explored
#
# The pyramid object stores the progress by every layer
# Layer keys are the number of rooms still left to solve at the corresponding layer
# Then for every layer we keep the generator and a list with every generated spot
#
# This script includes auxiliar logic for this coordinated paralelization strategy
# Most of this logic was reimplemented, improved and tested by Claude

# This function starts searches in the pyramid from the top to the bottom
# When it finds a generator which is still alive, it throws new digger processes from it which run in paralel
# However this function is not paralelized but it runs in the main process only
def paralelizer (
    starting_configuration : 'Room',
    rooms_to_solve : List['Room'],
    verbose : bool = False,
) -> bool:
    # Get the first room as the starting point
    current_room = rooms_to_solve[0]
    # Setup the manager which holds the pyramid of possible room configurations
    # Note that the pyramid lives inside this manager process and is never sent to other processes
    # Instead every process holds a proxy to it, so all of them claim spots from the very same pyramid
    pyramid_manager = PyramidManager()
    pyramid_manager.start()
    pyramid = pyramid_manager.Pyramid()
    # Set the root level of the pyramid, which is the only level the main process creates
    # Note that all further levels are created by the diggers as they go down the pyramid
    target_parent_grid, fitting_grid = starting_configuration._get_child_fitting_space(current_room, verbose)
    pyramid.add_level(starting_configuration, rooms_to_solve, 0, target_parent_grid, fitting_grid)

    # Setup the queue to handle results from multiple processes at once
    manager = Manager()
    pool_results = manager.Queue()

    # Set a flag to be returned when there is no result in the queue yet
    # Note that this can not be None, since None is the result reported by a process which failed
    no_result = 'no result'

    def get_result_if_exists() -> Union[str, None, 'Room']:
        try: return pool_results.get_nowait()
        except Empty: return no_result

    def wait_for_result() -> Optional['Room']:
        return pool_results.get()

    def show_error(error):
        # Get the whole traceback and not just the error message
        # Note that the traceback of an error raised in a different process comes as the error cause
        # This is why we must format the whole exception chain and not only the exception itself
        traceback_text = ''.join(format_exception(type(error), error, error.__traceback__))
        print(f'There was an error: {error}')
        print(traceback_text)
        # Report the crash to the main loop so it aborts the whole search
        # Note that this callback is run in the main process, but in a different thread
        # Thus the queue is used to make the main loop stop waiting for results which are not coming
        pool_results.put(DiggerError(traceback_text))

    # Set the search which will be run by every process
    # Note that the results queue is passed along, so every process can report on its own
    # Note that the pyramid proxy is passed along as well, so every process claims spots from it
    digger_task = DiggerTask(verbose=verbose, results_queue=pool_results, pyramid=pyramid)

    # Instantiate the pool to run multiple processes in paralel
    # Note that the pyramid manager is included here so it is shut down whatever happens
    # Otherwise its process would stay alive forever, since it is a process of its own
    with pyramid_manager, Pool(processes=N_PARALEL) as pool:

        # Handle a result which has been reported by any of the running processes
        # If the result is a successful configuration then save it and stop every other process
        # Return True when a configuration was saved and False when there is nothing to do yet
        # Note that this function is defined here since it requires the pool to exist already
        def handle_result (result) -> bool:
            # If any process crashed then we must abort the whole search
            # Note that we can not simply ignore it and carry on with the other processes
            # This is because the crashed process leaves a whole region of the pyramid unexplored
            # Thus we would think the pyramid was fully explored when it was not
            if type(result) == DiggerError:
                pool.terminate()
                pool.join()
                raise RuntimeError(f'A digger process crashed:\n{result.traceback_text}')
            # There is nothing to do if no process reported yet or if the reporting process failed
            if result is no_result or result is None:
                return False
            # Save the successful configuration in the room which is being solved
            # Note that the configuration was solved in a different process, so it must be pasted
            starting_configuration.paste(result)
            starting_configuration.update_display(title='Found successful children distribution')
            # Now that we have a solution there is no point in keeping the other processes alive
            pool.terminate()
            pool.join()
            return True

        # Launch as many processes as possible and wait for a positive result
        # If a process fails then launch a new process to fill the available CPU
        while True:
            # Before anything, check if we already had a result from any already running process
            if handle_result(get_result_if_exists()):
                return True
            # Check if there are CPUs available and, if so, launch a new process
            # Note that the pool cache holds the tasks which have not been completed yet
            running_processes_count = len(pool._cache)
            print(f'RUNNING PROCESSES {running_processes_count}')
            has_available_cpus = running_processes_count < N_PARALEL
            if has_available_cpus:
                # Claim a spot as high in the pyramid as possible
                # This way every new process explores a configuration very different from the others
                next_pyramid_claim = pyramid.claim_shallowest_spot()
                # If there are no more spots then we can not launch further processes
                if next_pyramid_claim is None:
                    # If there are no processes left then the whole pyramid was explored and we failed
                    # Note that both conditions are required, since a running process may add new levels
                    if running_processes_count == 0:
                        return False
                    # Otherwise we must wait for the already running processes to report
                    if handle_result(wait_for_result()):
                        return True
                    continue
                # Unpack some values
                level_id, distribution, remaining_rooms, depth, level_parent_grid, spot = next_pyramid_claim
                # Build the arguments to try this spot out of the claimed distribution
                # Note that the tried grids are not shared with the process which owns this level
                # Thus some work may be repeated, but only within this very spot
                next_room = remaining_rooms[0]
                following_rooms = remaining_rooms[1:]
                room_copy, child_spot_args = distribution._get_child_spot_args(
                    next_room, spot, level_parent_grid, set(), verbose)
                # Launch a new process
                # Note that the arguments must be picklable, since they are sent to another process
                # This is the reason why the pyramid is not sent but a proxy to it
                # Note that the depth tells how high in the pyramid this new search starts
                print(f'LAUNCHING PROCESS (pyramid depth {depth})')
                digger_args = (room_copy, child_spot_args, following_rooms, depth)
                pool.apply_async(digger_task.run, digger_args, error_callback=show_error)
                continue
            # Await for the next result
            if handle_result(wait_for_result()):
                return True

# This is what is reported to the main process when one of the paralel processes crashes
# Note that a crash is not a failure to find a solution, but a bug which must stop the whole process
# WARNING: The traceback must be sent as text, since traceback objects can not be pickled
class DiggerError:
    def __init__ (self, traceback_text : str):
        self.traceback_text = traceback_text

# A single level of the pyramid of possible room configurations
# It holds the distribution it starts from and a generator of raw spots to place the next room
# Note that every spot in the generator leads to a different configuration, thus to a different branch
class PyramidLevel:
    def __init__ (self,
        # The room configuration this level starts from
        distribution : 'Room',
        # The rooms which are still to be placed, being the first one the room to be placed here
        remaining_rooms : List['Room'],
        # How deep in the pyramid this level is, being 0 the root
        depth : int,
        # The parent grid where the next room is to be placed
        target_parent_grid : 'Grid',
        # The space in the parent grid where the next room actually fits
        fitting_grid : 'Grid',
    ):
        self.distribution = distribution
        self.remaining_rooms = remaining_rooms
        self.depth = depth
        self.target_parent_grid = target_parent_grid
        # Set the generator of raw spots for the next room to be placed
        # Note that this generator is created here, thus inside the manager process
        # This is essential since generators can not be pickled and thus they can not be sent between processes
        next_room = remaining_rooms[0]
        self.spots_generator = fitting_grid.generate_fitting_spots(next_room.min_size, next_room.min_size)

# The pyramid of possible room configurations, which is shared between every process
# Both the main process and the paralel processes claim spots from it to explore new branches
# Note that a spot is a destructive read, so every branch of the pyramid has one single explorer
# WARNING: This object is never sent between processes, since generators can not be pickled
# WARNING: Instead it lives inside the manager process and every other process holds a proxy to it
class Pyramid:
    def __init__ (self):
        # Keep every level by its id
        self.levels = {}
        # Set a counter to give a different id to every new level
        self.level_count = 0
        # WARNING: The manager runs a different thread per connected process
        # WARNING: Thus every access to the levels must be protected or two processes may claim a same spot
        self.lock = Lock()

    # Add a new level to the pyramid and return its id
    # Note that the fitting grid is calculated by the caller, which already has the distribution at hand
    def add_level (self,
        distribution : 'Room',
        remaining_rooms : List['Room'],
        depth : int,
        target_parent_grid : 'Grid',
        fitting_grid : 'Grid',
    ) -> int:
        with self.lock:
            level_id = self.level_count
            self.level_count += 1
            self.levels[level_id] = PyramidLevel(
                distribution, remaining_rooms, depth, target_parent_grid, fitting_grid)
            return level_id

    # Claim the next spot from a specific level
    # Return None when the level has no more spots to be claimed, thus it is fully explored
    # Note that this is used by the process which owns the level, since it already has the distribution
    def claim_spot (self, level_id : int) -> Optional['Rect']:
        with self.lock:
            level = self.levels.get(level_id, None)
            # The level may have been retired already by whoever consumed its last spot
            if level is None:
                return None
            spot = next(level.spots_generator, None)
            if spot is None:
                self._retire_level(level_id)
                return None
            return spot

    # Claim the next spot from the shallowest level which is still alive
    # Return None when there is no level left with spots to be claimed
    # Note that this is used by the main process to launch new processes as high in the pyramid as possible
    # This way every new process explores a very different configuration from the already running ones
    def claim_shallowest_spot (self) -> Optional[tuple]:
        with self.lock:
            # Sort levels by depth, so we always dig as high in the pyramid as possible
            for level_id in sorted(self.levels, key=lambda level_id: self.levels[level_id].depth):
                level = self.levels[level_id]
                spot = next(level.spots_generator, None)
                # If this level has no more spots then retire it and try with the next level
                if spot is None:
                    self._retire_level(level_id)
                    continue
                return level_id, level.distribution, level.remaining_rooms, level.depth, level.target_parent_grid, spot
            return None

    # Remove a level whose spots are over so we do not keep its distribution in memory forever
    # Otherwise the pyramid would hold the whole explored region, which may be huge
    # Note that the process working under this level keeps its own copy of the distribution
    # WARNING: This must be called with the lock already taken
    def _retire_level (self, level_id : int):
        del self.levels[level_id]

# Set a manager to hold the pyramid, so it is shared between every process
# Note that the registered class is instantiated inside the manager process and never sent
# WARNING: The registration must happen at the module level, since the manager process imports this module
class PyramidManager(BaseManager):
    pass
PyramidManager.register('Pyramid', Pyramid)

# This is the search which is run by every paralel process while solving a room children
# It searches the pyramid of possible room configurations from the bottom to the top
# Note that this must be a class at the module level and not a function inside the room solving logic
# This is because processes communicate by pickling their arguments and pickle can not handle local objects
# However pickle does handle instances of module level classes, and thus their methods as well
# Then every value which the search requires and is not a room may be stored here as an attribute
class DiggerTask:
    def __init__ (self,
        verbose : bool = True,
        # The queue where results are sent back to the main process
        # Note that this must be a queue from a multiprocessing Manager, since a regular queue is not picklable
        results_queue = None,
        # The proxy to the pyramid, which is shared between every process
        pyramid = None,
    ):
        self.verbose = verbose
        self.results_queue = results_queue
        self.pyramid = pyramid

    # This is the function which is actually sent to the paralel processes
    # It runs the search and reports the result back to the main process through the queue
    def run (self,
        rooms_distribution : 'Room',
        child_spot_args : tuple,
        remaining_rooms : List['Room'],
        depth : int,
    ):
        # Note that errors are not caught here on purpose
        # A search which finds no solution must report it as a regular result, thus returning None
        # However a search which crashes must abort the whole solving process and not be mistaken for a failure
        # Otherwise we would consider a whole region of the pyramid as explored when it was not
        # Note that a crash makes the pool run its error callback, which reports the crash to the main process
        result = self.dig(rooms_distribution, child_spot_args, remaining_rooms, depth)
        # Note that a single result is reported per process, when the whole search is over
        self.results_queue.put(result)

    # Try to place all rooms starting from a distribution which may be at any layer of the pyramid
    # The child spot arguments are the new spot for the next child to be placed in the distribution
    # The depth is the level of the pyramid this distribution comes from, being 0 the root
    # Return the solved configuration when it succeeds and None when it fails
    # WARNING: We can not return just a boolean since this may be running in a different process
    # WARNING: Thus the solved configuration must travel back to the main process to be useful
    def dig (self,
        rooms_distribution : 'Room',
        child_spot_args : tuple,
        remaining_rooms : List['Room'],
        depth : int,
    ) -> Optional['Room']:
        # Fit the child out of its initial grid
        if not rooms_distribution._try_initial_grid(*child_spot_args):
            return None
        # If there are no remaining rooms then we have a complete configuration here
        if len(remaining_rooms) == 0:
            # Try to use this configuration by adjusting the corners and adding the corridor
            # If the corridor does not fit then this branch of the pyramid is over
            if rooms_distribution._try_rooms_configuration():
                return rooms_distribution
            return None
        # If there are still rooms to solve then we run this same logic recursively
        # Get the first remaining room as the next room to be placed
        next_room = remaining_rooms[0]
        # Get the remaining rooms
        following_rooms = remaining_rooms[1:]
        # Find the parent space where the next room fits
        # Note that this is calculated here since we already have the distribution at hand
        target_parent_grid, fitting_grid = rooms_distribution._get_child_fitting_space(next_room, self.verbose)
        # Add the new level of spots to the shared pyramid
        # Note that the spots generator is created inside the manager process, where the pyramid lives
        # This way the main process may claim spots from this very level to launch further processes
        next_depth = depth + 1
        level_id = self.pyramid.add_level(
            rooms_distribution, remaining_rooms, next_depth, target_parent_grid, fitting_grid)
        # Keep track already tried initial grids so we do not loose time repeating the same
        tried_fitted_initial_grids = set()
        # Claim spots from our own level until there are no more spots to be claimed
        # Note that a claimed spot is never claimed again, not even by a different process
        # Thus every branch of the pyramid is explored by one single process and never twice
        while True:
            spot = self.pyramid.claim_spot(level_id)
            # If there are no more spots then this whole region of the pyramid is over
            if spot is None:
                return None
            # Build the arguments to try this spot from our own copy of the distribution
            # Note that only the raw spot travels between processes, which is way lighter
            next_rooms_distribution, next_child_spot_args = rooms_distribution._get_child_spot_args(
                next_room, spot, target_parent_grid, tried_fitted_initial_grids, self.verbose)
            # Call this function recursively for this spot
            solved_configuration = self.dig(
                next_rooms_distribution, next_child_spot_args, following_rooms, next_depth)
            if solved_configuration:
                return solved_configuration