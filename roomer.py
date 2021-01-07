from grid import *
from grid_display import *

from scheme import *
from scheme_display import plot_lines

# CLASSES ----------------------------------------------------------

class floor:
    # Set initial attributes
    def __init__(self, name, corners, door, parentRoom):
        self.name = name
        self.corners = corners
        self.door = door
        self.parentRoom = parentRoom

class room:
    # Set initial attributes
    def __init__(self, name, size, minWidth, maxWidth=math.inf, priorizeBorder=0, childRooms=[]):
        self.name = name
        self.size = size
        self.minWidth = minWidth
        self.maxWidth = maxWidth
        self.childRooms = childRooms
        self.priorizeBorder = priorizeBorder
    # Get the greatest common divisor of all minWidths
    def getGCD(self):
        mins = [self.minWidth]
        for c in self.childRooms:
            mins.append(c.getGCD())
        if (len(mins) == 1):
            return mins[0]
        else:
            return functools.reduce(math.gcd, mins)
    

# MAIN CODE ------------------------------------------------------------

def grid_distribution (distribution):

    # PRUEBAS
    matrix = Matrix(distribution.corners, distribution.parentRoom.getGCD(), updater = addFrame)
    #matrix = Matrix(distribution.corners, distribution.parentRoom.getGCD(), tracked=1)

    matrix.setArea(Cell(0,10),1)
    #matrix.setCluster(Cluster(2,0.1,3,3), 1)

    clusters = []

    # El pasillo
    parent = distribution.parentRoom

    corridor_cluster =  Cluster (2, parent.size, parent.minWidth, parent.maxWidth, parent.priorizeBorder, 2)

    matrix.setCluster(corridor_cluster, 1, origin = distribution.door)
    #matrix.fixCluster(corridor_cluster, 1)

    setupRepresentProcess()

    room_clusters = []
    count = 3
    for room in parent.childRooms:
        room_cluster = Cluster(count, room.size, room.minWidth, room.maxWidth, room.priorizeBorder)
        matrix.setCluster(room_cluster, 1)
        room_clusters.append(room_cluster)
        count += 1

    clusters = room_clusters
    #clusters = [corridor_cluster] + room_clusters
    print(clusters)
    for cluster in clusters:
        matrix.fixCluster(cluster, 1)
        
    for i in range(8):
        if matrix.getCluster(i):
            print(str(matrix.getCluster(i).size) + " -> " + str(matrix.countCells(i)));
        else:
            print('free -> ' + str(matrix.countCells(i)))

    print('yasta')


def vectors_distribution (distribution):

    # Set the limits of the whole scheme
    limits = Perimeter.from_corners(distribution.corners)

    # Set the base scheme
    #scheme = Scheme(limits, display = True)
    scheme = Scheme(limits)

    test1 = Rect(Point(-10,-10),Point(10,10))
    test2 = Rect(Point(5,-15),Point(25,5))

    line1 = Line(Point(-20, -20), Point(-20, 10))
    line2 = Line(Point(-20, -20), Point(-10, -20))

    #test3 = subtract_rects(test1, test2)
    test3 = limits.split_in_rectangles()

    lines = []
    for rect in test3:
        for line in rect.get_lines():
            lines.append(line)

    test4 = limits.area
    print(test4)

    #plot_lines([ *limits.lines, *test1.get_lines(), *test2.get_lines() ])
    plot_lines([ *limits.lines, *lines ])

    print('yasta')
    
        
# INPUTS ------------------------------------------------------------
             
distribution = floor('floor0',
    corners = [
        Point(-60,-60),
        Point(-60,+60),
        Point(+60,+60),
        Point(+60,-40),
        Point(+20,-40),
        Point(+20,-60)],
    door = Cell(0,10),
    parentRoom = room('pasillo', 0.1, 15, maxWidth=15, priorizeBorder=-1, childRooms=[
        room('comedor', 0.25, 30),
        room('cocina', 0.15, 30),
        room('habitacion1', 0.2, 30),
        room('habitacion2', 0.2, 30),
        room('baño', 0.1, 20),
]))


# This line is for windows to dont loop
if __name__ == '__main__':

    #grid_distribution(distribution)
    vectors_distribution(distribution)
