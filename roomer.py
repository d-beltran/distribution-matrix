from vectorial_base import *
from scheme import *
from scheme_display import setup_display

# Import some predefined test polygons
from test_boundaries import *

# This is for windows to dont loop
if __name__ == '__main__':

    # Set if we want to display the solving process
    display = True

    # Represent current rooms tagged as display = True
    if display:
        setup_display()

    test = Room(
        #boundary=test_boundary_3,
        max_corners = 6,
        display=True,
        name='Planta',
        doors=[ Door() ],
        #doors=[ Door(), Door() ],
        children=[
            Room(
                #polygon=Polygon.from_corners(test_polygon_6),
                forced_area=3600,
                #forced_area='37.5%',
                min_size=15,
                name='Comedor',
                fill_color='blue'),
            Room(
                forced_area=1600,
                #forced_area='18%',
                min_size=15,
                name='Cocina',
                fill_color='yellow'),
            Room(
                forced_area=1800,
                #forced_area='21%',
                display=True,
                min_size=15,
                name='Habitación 1',
                fill_color='red',
                children=[
                    # Room(
                    #     forced_area='25%',
                    #     min_size=10,
                    #     name='Lavabo de habitación',
                    #     fill_color='purple',
                    # )
                ]),
            Room(
                forced_area=1800,
                #forced_area='16.50%',
                min_size=15,
                name='Habitación 2',
                fill_color='orange'),
            Room(
                forced_area=800,
                #forced_area='7%',
                min_size=10,
                name='Lavabo',
                fill_color='green')
        ]
    )

    solve(test, display=display)

    print('Done!')