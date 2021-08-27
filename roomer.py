from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some predefined test polygons
from test_boundaries import *

# This segment is for windows to dont loop
if __name__ == '__main__':

    # Represent current rooms tagged as display = True
    setup_display()

    distribution = Room(
        boundary=test_boundary_4,
        display=True,
        name='Planta',
        children=[
            Room(
                #polygon=Polygon.from_corners(test_polygon_6),
                #forced_area=3600,
                forced_area='37.5%',
                min_size=15,
                name='Comedor',
                fill_color='blue'),
            Room(
                #forced_area=1600,
                forced_area='18%',
                min_size=15,
                name='Cocina',
                fill_color='yellow'),
            Room(
                #forced_area=1800,
                forced_area='18.75%',
                min_size=15,
                name='Habitación 1',
                fill_color='red'),
            Room(
                #forced_area=1800,
                forced_area='18.75%',
                min_size=15,
                name='Habitación 2',
                fill_color='orange'),
            Room(
                #forced_area=800,
                forced_area='7%',
                min_size=15,
                name='Lavabo',
                fill_color='green')
        ]
    )

    print('Done!')