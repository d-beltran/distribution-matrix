from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some predefined test perimeters
from test_perimeters import *


room_1 = Room(
    perimeter=Perimeter.from_corners(test_perimeter_7),
    forced_area=3600,
    min_size=15,
    name='Comedor',
    fill_color='blue')
room_2 = Room(
    perimeter=Perimeter.from_corners(test_perimeter_8),
    forced_area=1600,
    min_size=15,
    name='Cocina',
    fill_color='yellow')

#print(room_1.get_frontiers(room_2))

segment_1 = Segment(Point(0,1), Point(0,5))
segment_2 = Segment(Point(0,2), Point(0,4))

#print(segment_1.substract_segments([segment_2]))

segment_1 = Segment(Point(-40,-40), Point(-40,40))
segment_2 = Segment(Point(-40,40), Point(-40,-40))

#print(segment_1.substract_segments([segment_2]))

segment_1 = Segment(Point(-40,-40), Point(-40,40))
segment_2 = Segment(Point(-40,40), Point(-40,-40))
segment_3 = Segment(Point(-40,-40), Point(-40,-60)) # Este segmento hace que falle

#print(segment_1.substract_segments([segment_2, segment_3]))

segment_1 = Segment(Point(19.8201,-20.1856), Point(19.8201,-20.0))
segment_2 = Segment(Point(-20,-20), Point(20,-20))

print(segment_1.get_intersection_point(segment_2))