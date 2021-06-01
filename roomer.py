from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some predefined test perimeters
from test_perimeters import *

# This segment is for windows to dont loop
if __name__ == '__main__':

    # Represent current rooms tagged as display = True
    setup_display()

    distribution = Room(
        perimeter=Perimeter.from_corners(test_perimeter_4),
        display=True,
        name='Planta',
        children=[
            Room(
                #perimeter=Perimeter.from_corners(test_perimeter_6),
                forced_area=3600,
                min_size=15,
                name='Comedor',
                fill_color='blue'),
            Room(
                forced_area=1600,
                min_size=15,
                name='Cocina',
                fill_color='yellow'),
            Room(
                forced_area=1800,
                min_size=15,
                name='Habitación',
                fill_color='red'),
            Room(
                forced_area=1800,
                min_size=15,
                name='Habitación',
                fill_color='orange'),
            Room(
                forced_area=800,
                min_size=10,
                name='Lavabo',
                fill_color='green')
        ]
    )

    # room_1 = Room(
    #     perimeter=Perimeter.from_corners(test_perimeter_7),
    #     forced_area=3600,
    #     min_size=15,
    #     name='Comedor',
    #     fill_color='blue')
    # room_2 = Room(
    #     perimeter=Perimeter.from_corners(test_perimeter_8),
    #     forced_area=1600,
    #     min_size=15,
    #     name='Cocina',
    #     fill_color='yellow')

    # print(room_1.get_frontiers(room_2))