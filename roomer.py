from distributionMatrix import *
from matrixDisplay import *

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
        
    
# INPUTS ------------------------------------------------------------
             
distribution = floor('floor0',
    corners = [
        (-60,-60,0),
        (-60,+60,0),
        (+60,+60,0),
        (+60,-40,0),
        (+20,-40,0),
        (+20,-60,0)],
    door = Cell(0,10),
    parentRoom = room('pasillo', 0.1, 15, maxWidth=15, priorizeBorder=-1, childRooms=[
        room('comedor', 0.25, 30),
        room('cocina', 0.15, 30),
        room('habitacion1', 0.2, 30),
        room('habitacion2', 0.2, 30),
        room('baño', 0.1, 20),
]))
    
# MAIN CODE ------------------------------------------------------------

# This line is for windows to dont loop
if __name__ == '__main__':

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

