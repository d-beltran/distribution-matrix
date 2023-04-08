from scheme import *
from vectorial_base import *

# Perímetro en L
test_boundary_1 = Boundary(Polygon.from_corners([
    Point(-60,-60),
    Point(-60,+60),
    Point(+60,+60),
    Point(+60,-40),
    Point(+20,-40),
    Point(+20,-60)
]))

# Perímetro en Z
test_boundary_2 = Boundary(Polygon.from_corners([
    Point(-60,-60),
    Point(-60,+40),
    Point(-40,+40),
    Point(-40,+60),
    Point(+60,+60),
    Point(+60,-40),
    Point(+20,-40),
    Point(+20,-60)
]))

# Perímetro en peineta
test_boundary_3 = Boundary(Polygon.from_corners([
    Point(-60,-60),
    Point(-60,+20),
    Point(-20,+20),
    Point(-20,+60),
    Point(+20,+60),
    Point(+20,-20),
    Point(+60,-20),
    Point(+60,-60),
]))

# perímetro en +
test_boundary_4 = Boundary(Polygon.from_corners([
    Point(-60,-20),
    Point(-60,+20),
    Point(-40,+20),
    Point(-40,+40),
    Point(-20,+40),
    Point(-20,+60),
    Point(+20,+60),
    Point(+20,+40),
    Point(+40,+40),
    Point(+40,+20),
    Point(+60,+20),
    Point(+60,-20),
    Point(+40,-20),
    Point(+40,-40),
    Point(+20,-40),
    Point(+20,-60),
    Point(-20,-60),
    Point(-20,-40),
    Point(-40,-40),
    Point(-40,-20),
]))

# Perímetro en U
test_boundary_5 = Boundary(Polygon.from_corners([
    Point(-60,-60),
    Point(-60,+60),
    Point(-20,+60),
    Point(-20,-20),
    Point(+20,-20),
    Point(+20,+20),
    Point(+60,+20),
    Point(+60,-60),
]))

# Rectángulo 1
test_boundary_6 = Boundary(Polygon.from_corners([
    Point(-40,-40),
    Point(-40,+0),
    Point(0,0),
    Point(0,-40),
]))

# Rectángulo 2
test_boundary_7 = Boundary(Polygon.from_corners([
    Point(-80,-80),
    Point(-80,-40),
    Point(-40,-40),
    Point(-40,-80),
]))

# Rectángulo 3
test_boundary_8 = Boundary(Polygon.from_corners([
    Point(-40,-80),
    Point(-40, 0),
    Point(0,0),
    Point(0,-80),
]))

# Test distribución de habitaciones para pasillo complejo
# Espacio libre en el centro y necesidad de establecer pasillo a ambos lados del espacio libre
test_room_1 = Room(
    name='Planta',
    boundary=Boundary(Polygon.from_corners([
        Point(0,0),
        Point(0,60),
        Point(60,60),
        Point(60,40),
        Point(100,40),
        Point(100,0),
    ])),
    corridor_size=7,
    door_args={
        'width': 5,
        'margin': 1
    },
    doors=[ Door(point=Point(0,30)) ],
    children=[
        Room(
            boundary=Boundary(Polygon.from_corners([
                Point(0,0),
                Point(0,30),
                Point(20,30),
                Point(20,0),
            ])),
            name='Hab 1',
            fill_color='blue'),
        Room(
            boundary=Boundary(Polygon.from_corners([
                Point(0,30),
                Point(0,60),
                Point(20,60),
                Point(20,30),
            ])),
            name='Hab 2',
            fill_color='yellow'),
        Room(
            boundary=Boundary(Polygon.from_corners([
                Point(20,40),
                Point(20,60),
                Point(60,60),
                Point(60,40),
            ])),
            name='Hab 3',
            fill_color='red'),
        Room(
            boundary=Boundary(Polygon.from_corners([
                Point(60,0),
                Point(60,20),
                Point(80,20),
                Point(80,0),
            ])),
            name='Hab 4',
            fill_color='green'),
        Room(
            boundary=Boundary(Polygon.from_corners([
                Point(60,20),
                Point(60,40),
                Point(80,40),
                Point(80,20),
            ])),
            name='Hab 5',
            fill_color='cyan'),
        Room(
            boundary=Boundary(Polygon.from_corners([
                Point(80,0),
                Point(80,40),
                Point(100,40),
                Point(100,0),
            ])),
            name='Hab 6',
            fill_color='orange',
            children=[
                Room(
                    boundary=Boundary(Polygon.from_corners([
                        Point(80,0),
                        Point(80,10),
                        Point(100,10),
                        Point(100,0),
                    ])),
                    name='Hab 7',
                    fill_color='grey'),
            ]),
    ]
)

# Nada establecido
test_room_2 = Room(
        #boundary=test_boundary_3,
        max_corners = 6,
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
                min_size=15,
                name='Habitación 1',
                fill_color='red',
                children=[
                    Room(
                        forced_area='25%',
                        min_size=10,
                        name='Lavabo de habitación',
                        fill_color='purple',
                    )
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

# Polygono exterior preestablecido
test_room_3 = Room(
        boundary=test_boundary_5,
        name='Planta',
        doors=[ Door() ],
        children=[
            Room(
                forced_area='37.5%',
                min_size=15,
                name='Comedor',
                fill_color='blue'),
            Room(
                forced_area='18%',
                min_size=15,
                name='Cocina',
                fill_color='yellow'),
            Room(
                forced_area='21%',
                min_size=15,
                name='Habitación 1',
                fill_color='red',
                children=[
                    Room(
                        forced_area=200,
                        min_size=10,
                        name='Lavabo de habitación',
                        fill_color='purple',
                    )
                ]),
            Room(
                forced_area='16.50%',
                min_size=15,
                name='Habitación 2',
                fill_color='orange'),
            Room(
                forced_area='7%',
                min_size=10,
                name='Lavabo',
                fill_color='green')
        ]
    )

# Familar house with 2 floors
test_building_1 = Building(
    floors = {
        0: Room(
            name='Primera planta',
            max_corners = 6,
            children = [
                Room(
                    name='Comedor',
                    forced_area=3600,
                    min_size=15,
                    fill_color='blue'
                ),
                Room(
                    name='Cocina',
                    forced_area=2000,
                    min_size=15,
                    fill_color='yellow'
                ),
            ]
        ),
        1: Room(
            name='Segunda planta',
            children = [
                Room(
                    name='Habitación 1',
                    forced_area=1800,
                    min_size=15,
                    fill_color='red',
                    children=[
                        Room(
                            name='Lavabo de habitación',
                            forced_area='20%',
                            min_size=10,
                            fill_color='purple',
                        )
                    ]),
                Room(
                    name='Habitación 2',
                    forced_area=1800,
                    min_size=15,
                    fill_color='orange'),
                Room(
                    name='Lavabo',
                    forced_area=800,
                    min_size=10,
                    fill_color='green')
            ]
        )
    },
    floor_heights = 20
)