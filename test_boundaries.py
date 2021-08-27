from vectorial_base import Point, Polygon, Boundary

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