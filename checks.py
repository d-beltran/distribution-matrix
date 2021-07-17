from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some predefined test perimeters
from test_perimeters import *

def test (actual_value, expected_value) -> str:
    if actual_value == expected_value:
        return print(str(actual_value) + ' -> ' + 'OK')
    return print(str(actual_value) + ' -> ' + 'ERROR')


# Should be false
# However, with rounded float to the 4th deciaml is true (experimentally proved)
point = Point(-22.4997, -19.999)
segment = Segment(Point(-20.0,-20.0),Point(-82.8538,-20.0))
test(point in segment, False)

# Should be false
# However, with the 'isclose' float comparision is true (experimentally proved)
point = Point(-44.9348, -11.8746)
segment = Segment(Point(-22.112,-11.8742),Point(-148.6031,-11.8742))
test(point in segment, False)

# Should be true
# However, with natural float comparision is false (experimentally proved)
point = Point(49.3939, -20.0)
segment = Segment(Point(49.3939, -1.2122), Point(49.3939, -60.0))
test(point in segment, True)

#setup_display()