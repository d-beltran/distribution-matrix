from typing import Union, Generator, Optional
from math import isclose

# Auxiliar functions widely used along the whole code

# Exception for when user input is wrong
class InputError(SystemExit):
    pass

# Float numbers are not allowed internally
number = Union[int, float]
def is_number(var):
    return isinstance(var, int) or isinstance(var, float)

# Set a resoltion limit
# Resolution is the number of decimals to take in count when comparing different coordinates
# When working with substractions (e.g. distances) the resolution should be decreased in 1 order
# When working with squared substractions (e.g. areas) the resolution should be decreased in 2 orders
base_resolution = 4
def resolute(num, base_offset = 0):
    base = base_resolution + base_offset
    resolution = 10**base
    return round(num * resolution) / resolution

# Calculate the minimum resolution
# i.e. changes below this resolution make no sense
minimum_resolution = 1 / 10 ** base_resolution
# Set functions to compare float according to the minimum resolution
def equal (a : float, b : float) -> bool:
    return a < b + minimum_resolution and a > b - minimum_resolution
def greater (a : float, b : float) -> bool:
    return a > b + minimum_resolution
def lower (a : float, b : float) -> bool:
    return a < b - minimum_resolution
def equal_or_greater (a : float, b : float) -> bool:
    return equal(a, b) or greater(a, b)
def equal_or_lower (a : float, b : float) -> bool:
    return equal(a, b) or lower(a, b)

# Set another function to check if two numbers are very close, out of resolution matters
# We do the double check because the 'isclose' function may fail for values close to 0
# see https://stackoverflow.com/questions/35324893/using-math-isclose-function-with-values-close-to-0
def same_number(a : number, b : number) -> bool:
    # DANI: Esto debería funcionar siempre pero son 4 operaciones
    #return isclose(a,b) or isclose(a+1,b+1)
    # Esto debería funcionar bien
    return isclose(a,b, abs_tol=minimum_resolution)

# Just for better display
def round_to_hundredths (number : number) -> number:
    return round(number * 100) / 100

# Return unique values of a list
def unique (values : list) -> list:
    return list(set(values))

# Set a special iteration system
# Return each value of the array and the next value
# By default, the final array value is skipped, since it has no next value
# However, if the 'retro' argument is True, the final array value is returned with the first array value as the next value
# By default, values are returned as follows: A with B, B with C, C with D ...
# However, if the 'loyals' argument is True, values are returned as follows: A with B, C with D, E with F ...
def pairwise (values : list, retro : bool = False, loyals : bool = False) -> Generator[tuple, None, None]:
    last = len(values) - 1
    step = 2 if loyals else 1
    for i in range(0, last, step):
        yield values[i], values[i+1]
    if retro:
        yield values[last], values[0]

# Set a special iteration system
# Return one value of the array and a new array with all other values for each value
def otherwise (values : list) -> Generator[tuple, None, None]:
    for v, value in enumerate(values):
        others = values[0:v] + values[v+1:]
        yield value, others

# Set a special iteration system
# Return one value of the array and a new array with the rest of values after for each value
def afterwise (values : list) -> Generator[tuple, None, None]:
    for v, value in enumerate(values):
        others = values[v+1:]
        yield value, others

# Get the first match index of a value in a list
# If there is no match then return None instead of crashing as the list.index() would
def get_index (values : list, target_value) -> Optional[int]:
    for v, value in enumerate(values):
        if value == target_value: return v
    return None