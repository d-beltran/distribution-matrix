from vectorial_base import *
from scheme import *
from scheme_display import plot_lines, plot_everything
        
# INPUTS ------------------------------------------------------------

corners0 = [
    Point(-60,-60),
    Point(-60,+60),
    Point(+60,+60),
    Point(+60,-40),
    Point(+20,-40),
    Point(+20,-60)
]

corners1 = [
    Point(-60,-60),
    Point(-60,+40),
    Point(-40,+40),
    Point(-40,+60),
    Point(+60,+60),
    Point(+60,-40),
    Point(+20,-40),
    Point(+20,-60)
]

corners2 = [
    Point(-60,-60),
    Point(-60,+20),
    Point(-20,+20),
    Point(-20,+60),
    Point(+20,+60),
    Point(+20,-20),
    Point(+60,-20),
    Point(+60,-60),
]

corners3 = [
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
]

corners4 = [
    Point(-60,-60),
    Point(-60,+60),
    Point(-20,+60),
    Point(-20,-20),
    Point(+20,-20),
    Point(+20,+20),
    Point(+60,+20),
    Point(+60,-60),
]

corners5 = [
    Point(-40,-40),
    Point(-40,+0),
    Point(0,0),
    Point(0,-40),
]

corners6 = [
    Point(-80,-80),
    Point(-80,-40),
    Point(-40,-40),
    Point(-40,-80),
]

distribution = Room(
    perimeter=Perimeter.from_corners(corners3),
    display=True,
    name='Planta',
    children=[
        Room(
            #perimeter=Perimeter.from_corners(corners5),
            forced_area=600,
            min_size=15,
            name='Habitación'),
        Room(
            forced_area=400,
            min_size=10,
            name='Lavabo')
    ]
)

#door = Point(-60,0),
#parentRoom = Room('pasillo', 0.1, 15, maxWidth=15, priorizeBorder=-1, childRooms=[
#    room('comedor', 0.25, 30),
#    room('cocina', 0.15, 30),
#    room('habitacion1', 0.2, 30),
#    room('habitacion2', 0.2, 30),
#    room('baño', 0.1, 20),


# This line is for windows to dont loop
if __name__ == '__main__':

    # Represent current rooms tagged as display = True
    setup_display()

    test = distribution.free_mrects