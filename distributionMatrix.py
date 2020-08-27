# Room distributor

# MODULE IMPORTS ------------------------------------------------------------
import math
import functools
import random  

# CLASS DEFINITIONS ------------------------------------------------------------
    
# An x,y coordinate
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __str__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'
    def __repr__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'
        
# A segment defined by 2 coordinates
class Line:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    def __str__(self):
        return 'A: ' + str(self.a) + ', B: ' + str(self.b)

# An x,y coordinate
class Cell:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        return False
    def __hash__(self):
        return hash((self.x, self.y))
    def __str__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'
    def __repr__(self):
        return '(x: ' + str(self.x) + ', y: ' + str(self.y) + ')'

# An inteligent group of cells
# value: The int value which cells will be set to
# size: The int number of cells or the float percent of the matrix size
# min and max widths: The minimum and maximum wide of the cluster in both dimension
# priorize border: Set if the cluster is meant to be in the matrix border or not
#   1 - Must be in the border
#   0 - Indiferent
#   -1 - Must be away from the border
# behaviour: Set how cluster must expand
#   0 - Caotic
#   1 - As a compact square
#   2 - As a long rectange
class Cluster:
    def __init__(self, value, size, minWidth, maxWidth = None, priorizeBorder = 0, behaviour = 0):
        self.value = value
        self.size = size
        self.minWidth = minWidth
        self.maxWidth = maxWidth
        self.priorizeBorder = priorizeBorder
        self.behaviour = behaviour
        # Define if the cluster has been set successfully
        self.set = False
    def __str__(self):
        return 'c-' + str(self.value)
    def __repr__(self):
        return 'c-' + str(self.value)

# Set a matrix with as cells with the specified cell size
# Some cells of the matrix may be out of the corners perimeter
# This will happend in not perfectly squared/rectangular perimeters
class Matrix:
    def __init__(self, corners, cellSize, tracked = 0, updater = None):
        self.corners = corners
        self.cellSize = cellSize
        self.tracked = tracked
        self.updater = updater
        # Get maximum and minimum x and y coordinates, then calculate the x and y ranges
        xcoords = []
        ycoords = []
        # Get all xs and ys in separated arrays
        for c in corners:
            xcoords.append(c[0])
            ycoords.append(c[1])
        # Get maximum and minimums
        self.maxx = max(xcoords)
        self.minx = min(xcoords)
        self.maxy = max(ycoords)
        self.miny = min(ycoords)
        # Set the cells ranges
        self.xrange = int((self.maxx - self.minx) / cellSize)
        self.yrange = int((self.maxy - self.miny) / cellSize)
        self.size = self.xrange * self.yrange
        # Set all cells as 0 at first
        self.values = [0] * self.size
        # Set the perimeter lines
        self.lines = self.setLines(corners)
        
    # Set the lines array between these coordinates
    def setLines (self, corners):
        lines = [ Line(
                Point(corners[len(corners)-1][0], corners[len(corners)-1][1]),
                Point(corners[0][0], corners[0][1])) ]
        for c in range(0, len(corners)-1):
            lines.append( Line(
                Point(corners[c][0], corners[c][1]),
                Point(corners[c+1][0], corners[c+1][1])))
        return lines
    
    # Set the values from the origin cell and all "connected" cells
    # Connected cells are cells from the same value not separated by lines
    # In addition, return the cells array
    def setArea (self, origin, value):
        area = [self.getCell(origin)]
        originValue = self.getValue(area[0])
        self.setValue(area[0], value)
        i = 0
        while(i < len(area)):
            cell = area[i]
            for c in self.getColliderCells(cell, value=originValue, diagonals=True, connected=True):
                area.append(c)
                self.setValue(c, value)
            i += 1
        return area
    
    # Get available cells in the matrix for the new cluster
    # Consider available cells whose value is equal to 'sourceValue'
    # The first cell to be placed is the 'origin', and this cell is never lost
    clusters = []
    def setCluster (self, cluster, sourceValue, origin = None):
        self.track('Setting new cluster ' + str(cluster), deepen = 1)
        # Check that there is not another cluster with the same value
        if (self.getCluster(cluster.value)):
            self.track('FAIL: Cluster ' + str(cluster) + ' already exists', deepen = -1)
            return False
        # Save the new cluster in the clusters list
        self.clusters.append(cluster)
        # If the cluster size is in float format it must be a percent
        # Reset it as number of cells
        if isinstance(cluster.size, float):
            # Solución chapucera: como numero total de source cells cojo las que no son 0
            # Si coges las que son = al sourceValue el numero se reduce con cada iteración
            sourceNumber = len(self.filterCells(list(range(self.size)), 0, False))
            #print(str(cluster.size) + " -> " + str(round(sourceNumber * cluster.size)))
            cluster.size = round(sourceNumber * cluster.size)
        # Adapt minimum and maximum sizes
        # WARNING: These values should be integers. There is no control about it
        if cluster.minWidth:
            cluster.minWidth = round(cluster.minWidth / self.cellSize)
        if cluster.maxWidth and not math.isinf(cluster.maxWidth):
            #print(str(cluster.maxWidth) + ' / ' + str(self.cellSize))
            cluster.maxWidth = round(cluster.maxWidth / self.cellSize)
        # Check that there is enough space
        if self.countCells(sourceValue) < cluster.size:
            self.track("FAIL: Not enought avilable cells to allocate the cluster ("+
            str(self.countCells(sourceValue)) + "/" + str(cluster.size) + ")", deepen = -1)
            return False
        # Set the cluster inital minimum area
        first = False
        for area in self.priorizeArea(sourceValue, cluster.minWidth, cluster.minWidth, origin):
            first = self.claimCellGroup(area, cluster.value, sourceValue)
            if first:
                break
        if first == False:
            self.track("FAIL: New cluster failed to get the origin cell", deepen = -1)
            return False
        # Try to expand the cluster cells to reach the cluster size
        # We substract the origin cell from the size count (cluster.size - 1)
        if self.expandCluster(cluster, sourceValue, expansion = cluster.size - cluster.minWidth * cluster.minWidth) == False:
            self.track("FAIL: New cluster failed to expand", deepen = -1)
            return False
        # Finally set the cluster as set and return True
        cluster.set = True
        self.track("SUCCESS: New cluster " + str(cluster) + " was set successfully", deepen = -1)
        return True
    
    # Get the number of cells with the specified value
    def countCells (self, value):
        #print(sum(v == value for v in self.values))
        return sum(v == value for v in self.values)
    
    # Get a random cell with the specified value (if so) or any value
    # If input is empty return None
    def getRandomCell (self, cells = None, value = None):
        if cells == None or len(cells) == 0:
            return None
        random.shuffle(cells)
        if value == None:
            return self.getCell(cells[0])
        else:
            for c in cells:
                if self.getValue(cells[c]) == value:
                    #print (self.getCell(cells[c]))
                    return self.getCell(cells[c])
        return None        

    # Get a random element from the specified input array
    def getRandom (self, inputs = None):
        if inputs == None or len(inputs) == 0:
            return None
        random.shuffle(inputs)
        return inputs[0]

    # Return the left cell. If there is no left cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getLeftCell (self, cell):
        x = cell.x - 1
        if x >= 0:
            return Cell(x, cell.y)
        return None

    # Return the right cell. If there is no right cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getRightCell (self, cell):
        x = cell.x + 1
        if x < self.xrange:
            return Cell(x, cell.y)
        return None

    # Return the down cell. If there is no down cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getDownCell (self, cell):
        y = cell.y - 1
        if y >= 0:
            return Cell(cell.x, y)
        return None

    # Return the up cell. If there is no up cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getUpCell (self, cell):
        y = cell.y + 1
        if y < self.yrange:
            return Cell(cell.x, y)
        return None

    # Return the left down cell. If there is no left down cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getLeftDownCell (self, cell):
        x = cell.x - 1
        y = cell.y - 1
        if x >= 0 and y >= 0:
            return Cell(x, y)
        return None

    # Return the right down cell. If there is no right down cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getRightDownCell (self, cell):
        x = cell.x + 1
        y = cell.y - 1
        if x < self.xrange and y >= 0:
            return Cell(x, y)
        return None

    # Return the left up cell. If there is no left up cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getLeftUpCell (self, cell):
        x = cell.x - 1
        y = cell.y + 1
        if x >= 0 and y < self.yrange:
            return Cell(x, y)
        return None

    # Return the right up cell. If there is no right up cell, return None.
    # The input cell must exist. Otherwise the returned cell may not exists.
    def getRightUpCell (self, cell):
        x = cell.x + 1
        y = cell.y + 1
        if x < self.xrange and y < self.yrange:
            return Cell(x, y)
        return None
    
    # Return cells which collide the specified cell
    # When 'diagonals' is true also de the 'corners' are taken
    # They can be filtered to be connected or from a specific value
    # Impossible cells (i.e. out of the matrix limits) are not returned
    def getColliderCells (self, cell, value = None, diagonals = True, connected = False):
        self.track('Getting cells colliding ' + str(cell), level = 1, deepen = 1)
        # Get all cells around
        colliders = [self.getLeftCell(cell),
                     self.getRightCell(cell),
                     self.getDownCell(cell),
                     self.getUpCell(cell)]
        if diagonals:
            colliders += [self.getLeftDownCell(cell),
                          self.getRightDownCell(cell),
                          self.getLeftUpCell(cell),
                          self.getRightUpCell(cell)]
        result = []
        for c in colliders:
            if( c
            and (value == None or self.getValue(c) == value)
            and (connected == False or self.connected(cell,c))):
                result.append(c)
        self.track('Cell colliding Results: ' + str(result), level = 1, deepen = -1)
        return result
    
    # Check if the specified cell is in the matrix limit or next to any edge
    def inBorder (self, cell):
        colliders = self.getColliderCells(cell, diagonals = True, connected = True)
        return len(colliders) < 8
    
    # Check if any of the specified cells is in the matrix limit or next to any edge
    def inBorderGroup (self, cells):
        for cell in cells:
            if self.inBorder(cell):
                return True
        return False
    
    # Get all cells colliding any cell of the specified cluster
    # By default corners are not taken
    def getClusterColliderCells (self, clusterValue, value = None, diagonals = False, connected = True):
        self.track('Getting cells colliding cluster ' + str(clusterValue), deepen = 1)
        result = []
        # Get the colliding cells of each cluster cells
        for c in self.filterCells(list(range(self.size)), clusterValue):
            for col in self.getColliderCells(c, diagonals = diagonals, connected = connected):
                result.append(col)
        # Select only cells which are not already from the cluster
        # Remove deplicates with the "-> set -> list" workaround
        result = list(set(self.filterCells(result, clusterValue, False)))
        self.track('Cluster collinding results: ' + str(result), deepen = -1)
        return result
    
    # When include = True, return only cells with the specified value
    # When include = False, return all cells but the ones with the specified value
    # A list of cells may be provided. Otherwise, return None
    # Input cells may be any format, but output is always cells
    def filterCells (self, cells, value, include = True):
        self.track('Filtering cells from: ' + str(cells), level = 1, deepen = 1)
        if cells == None:
            self.track('FILTER: here are no results', level = -1, deepen = -1)
            return None
            #cells = list(range(self.size))
        result = []
        if include == True:
            self.track('FILTER: Cell values must be ' + str(value), level = 1)
            for c in cells:
                if self.getValue(c) == value:
                    result.append(self.getCell(c))
        if include == False:
            self.track('FILTER: Cell values must not be ' + str(value), level = 1)
            for c in cells:
                if self.getValue(c) != value:
                    result.append(self.getCell(c))
        if len(result) == 0:
            self.track('FILTER: There are no results', level = 1, deepen = -1)
            return None
        self.track('Filter results: ' + str(result), level = 1, deepen = -1)
        return result
    
    # Get the suitable cells from a cells array
    # Suitable cells are thouse which dont belog to a 'banned' cluster
    def getSuitableCells (self, cells, source, banned = []):
        self.track('Getting suitable cells from: ' + str(cells), deepen = 1)
        suitables = cells
        # Filter those cells excluding each 'previous' values
        for p in banned:
            suitables = self.filterCells(suitables, p, False)
        # Return filtered cells. They may be 'None'
        self.track('Suitable results: ' + str(suitables), deepen = -1)
        return suitables
    
    # Count cells in each direction from the specified cell with the specified value
    # Also return if the specified cell belongs or not (true/false) to this value
    # The result is structured as follows:
    # (belongs, -x, +x, -y, +y, -x-y, +x-y, -x+y, +x+y)
    def getCellRelatives (self, cell, value):
        # Get the cell value
        cellValue = self.getValue(cell)
        # Count cells in the '-x' direction
        nx = 0
        for i in range(self.maxx):
            if self.getValue(Cell(cell.x - (nx + 1), cell.y)) == value:
                nx += 1
            else:
                break
        # Count cells in the '+x' direction
        px = 0
        for i in range(self.maxx):
            if self.getValue(Cell(cell.x + (px + 1), cell.y)) == value:
                px += 1
            else:
                break
        # Count cells in the '-y' direction
        ny = 0
        for i in range(self.maxy):
            if self.getValue(Cell(cell.x, cell.y - (ny + 1))) == value:
                ny += 1
            else:
                break
        # Count cells in the '+y' direction
        py = 0
        for i in range(self.maxy):
            if self.getValue(Cell(cell.x, cell.y + (py + 1))) == value:
                py += 1
            else:
                break
        # Count cells in the '-x-y' direction
        wider = max(self.maxx, self.maxy)
        nxny = 0
        for i in range(wider):
            if self.getValue(Cell(cell.x - (nxny + 1), cell.y - (nxny + 1))) == value:
                nxny += 1
            else:
                break
        # Count cells in the '+x-y' direction
        pxny = 0
        for i in range(wider):
            if self.getValue(Cell(cell.x + (pxny + 1), cell.y - (pxny + 1))) == value:
                pxny += 1
            else:
                break
        # Count cells in the '-x+y' direction
        nxpy = 0
        for i in range(wider):
            if self.getValue(Cell(cell.x - (nxpy + 1), cell.y + (nxpy + 1))) == value:
                nxpy += 1
            else:
                break
        # Count cells in the '+x+y' direction
        pxpy = 0
        for i in range(wider):
            if self.getValue(Cell(cell.x + (pxpy + 1), cell.y + (pxpy + 1))) == value:
                pxpy += 1
            else:
                break
        return (cellValue, nx, px, ny, py, nxny, pxny, nxpy, pxpy)

    # Evaluate the specified cell relatives to determine how suitable is to be claimed next
    def getCellPrioirty (self, cell, cluster, source):
        # Set the value, minimum and maximum values
        # The maximum may be null
        value = cluster.value
        minWidth = 0
        if cluster.minWidth:
            minWidth = cluster.minWidth
        maxWidth = 0
        if cluster.maxWidth:
            maxWidth = cluster.maxWidth
        # (belongs, -x, +x, -y, +y, -x-y, +x-y, -x+y, +x+y)
        relatives = self.getCellRelatives(cell, value)
        # If the cell already belongs to this cluster then return 0
        if relatives[0] == value:
            return 0
        # Set the priority, which measure how suitable is a cell to be claimed
        priority = 0
        # From 1 to 4
        for i in range(1,5):
            # 0 points when this side is 0
            if relatives[i] == 0:
                priority += 0 
            # 1 point when this side is exactly in the maximum limit
            elif relatives[i] == maxWidth:
                priority += 1
            # 10 points when this side is exactly in the minimum limit
            elif relatives[i] == minWidth:
                priority += 10
            # 100 points when this side has already crossed the maximum
            elif relatives[i] > maxWidth:
                priority += 100
            # 100000 points when this side has not reached the minimum yet
            elif relatives[i] < minWidth:
                priority += 100000
            # 1000 points if it is between the minimum and the maximum
            else:
                priority += 1000
        # Next, evaluate cells which join disconnected parts
        # If -x-y < minWidth and -x >= minWidth and -y >= minWidth
        if relatives[5] < minWidth and relatives[1] >= minWidth and relatives[3] >= minWidth:
            priority += 100000
        # If +x-y < minWidth and +x >= minWidth and -y >= minWidth
        if relatives[6] < minWidth and relatives[2] >= minWidth and relatives[3] >= minWidth:
            priority += 100000
        # If -x+y < minWidth and -x >= minWidth and +y >= minWidth
        if relatives[7] < minWidth and relatives[1] >= minWidth and relatives[4] >= minWidth:
            priority += 100000
        # If +x+y < minWidth and +x >= minWidth and +y >= minWidth
        if relatives[8] < minWidth and relatives[2] >= minWidth and relatives[4] >= minWidth:
            priority += 100000
        # Independently from the previous clasification, if the priority is not 0...
        # Add 10000 if the cell belongs to the source
        if relatives[0] == source and priority > 0:
            priority += 10000
        # Finally, return priority
        return priority          
    
    # Get an array with all rows in a cell array
    # The 'connected' argument may be set as true to get only connected rows
    # If so, multiple arrays from the same row may be returned if it is interrupted
    # Cells in every row are sortered
    # WARNING: A column of 'x' cells would be understood as 'x' rows of 1 cell
    def getRows (self, cells, connected = False, minimumLength = None):
        rows = []
        # Keep track of the currently calculated cells to avoid repeating searches
        searchedRows = []
        # Sorting function
        def byX (elem):
            return elem.x
        for c in cells:
            if c in searchedRows:
                continue
            newRow = [c]
            for i in range(self.maxx):
                newCell = Cell(c.x + i + 1, c.y)
                if newCell in cells:
                    newRow.append(newCell)
                    searchedRows.append(newCell)
                elif connected:
                    break
            for i in range(self.maxx):
                newCell = Cell(c.x - i - 1, c.y)
                if newCell in cells:
                    newRow.append(newCell)
                    searchedRows.append(newCell)
                elif connected:
                    break
            newRow.sort(key=byX)
            rows.append(newRow)
        if minimumLength:
            def minLength (row):
                return len(row) >= minimumLength
            rows = filter(minLength, rows)
        return rows
    
    # Get an array with all columns in a cell array
    # The 'connected' argument may be set as true to get only connected columns
    # If so, multiple arrays from the same column may be returned if it is interrupted
    # Cells in every column are sortered
    # WARNING: A row of 'x' cells would be understood as 'x' columns of 1 cell
    def getColumns (self, cells, connected = False, minimumLength = None):
        columns = []
        # Keep track of the currently calculated cells to avoid repeating searches
        searchedColumns = []
        # Sorting function
        def byY (elem):
            return elem.y
        for c in cells:
            if c in searchedColumns:
                continue
            newColumn = [c]
            for i in range(self.maxx):
                newCell = Cell(c.x, c.y + i + 1)
                if newCell in cells:
                    newColumn.append(newCell)
                    searchedColumns.append(newCell)
                elif connected:
                    break
            for i in range(self.maxx):
                newCell = Cell(c.x, c.y - i - 1)
                if newCell in cells:
                    newColumn.append(newCell)
                    searchedColumns.append(newCell)
                elif connected:
                    break
            newColumn.sort(key=byY)
            columns.append(newColumn)
        if minimumLength:
            def minLength (column):
                return len(column) >= minimumLength
            columns = filter(minLength, columns)
        return columns
    
    # Get an array with all possible cell groups in a cell array
    # Cell groups are multiple cells colliding in the same row/column
    # The minimum of cells in each group must be defined
    # The 'collapse' option is used to skip small groups when a bigger group contains its cells
    def getCellGroups (self, cells, minimum, maximum = None, collapse = False):
        # First of all, try to get complete rows/columns
        rows = self.getRows(cells, connected = True, minimumLength = minimum)
        columns = self.getColumns(cells, connected = True, minimumLength = minimum)
        # Define a function to convert an array of seried cells in groups of a given size
        def setGroups (serie, size):
            print('setting group with size ' + str(size) + ' in ' + str(serie))
            grps = []
            # Set the expected number of posible groups
            grpsNumber = len(serie) - size + 1
            if(grpsNumber > 0):
                # Find all posible groups
                for g in range(grpsNumber):
                    newGroup = []
                    for i in range(size):
                        newGroup.append(serie[g+i])
                    grps.append(newGroup)
            return grps
        # Get all possible groups in those rows and columns
        groups = []
        for r in rows:
            maxCells = len(r)
            if maximum == None or maximum > maxCells:
                maximum = maxCells
            if collapse:
                groups += setGroups(r, maximum)
            else:
                for n in range(minimum, maximum + 1):
                    groups += setGroups(r,n)
        for c in columns:
            maxCells = len(c)
            if maximum == None or maximum > maxCells:
                maximum = maxCells
            if collapse:
                groups += setGroups(c, maximum)
            else:
                for n in range(minimum, maximum + 1):
                    groups += setGroups(c,n)
        return groups
            
    # Find if the cluster is respecting each of the cluster rules
    # 1. Find if it has the minimum number of cells to fill the minimum widths
    # NOTA: Esto no se está usando de momento
    def getClusterStatus (self, cluster):
        print('getting status...')
        clusterCells = self.filterCells(list(range(self.size)), cluster.value)
        return len(clusterCells) >= cluser.minWidth * cluser.minWidth
    
    # Check that all cluster parameters are correct
    # 1 - Check that the cluster has as many cells as it is meant to
    # 2 - Check that all minimums are filled
    # 3 - Check that all parts are connected with the requiered minimum of cells
    def isClusterCorrect (self, cluster):
        print('checking cluster...')
        value = cluster.value
        clusterCells = self.filterCells(list(range(self.size)), value)
        # Check that the cluster has as many cells as it is meant to
        if len(clusterCells) != cluster.size:
            return False
        # Check that all minimums are filled
        rows = self.getRows(clusterCells, connected = True)
        for r in rows:
            if len(r) < cluster.minWidth:
                return False
        columns = self.getColumns(clusterCells, connected = True)
        for c in columns:
            if len(c) < cluster.minWidth:
                return False
        # Check that all parts are connected with the requiered minimum of cells
        def isRowConnected (row):
            count = 0
            for cell in row:
                leftCell = self.getLeftCell(cell)
                if leftCell and self.getValue(leftCell) == value:
                    count += 1
                else:
                    count = 0
                if count == minWidth:
                    return True
            count = 0
            for cell in row:
                rightCell = self.getRightCell(cell)
                if rightCell and self.getValue(rightCell) == value:
                    count += 1
                else:
                    count = 0
                if count == minWidth:
                    return True
            return False
        def isColumnConnected (column):
            count = 0
            for cell in column:
                downCell = self.getDownCell(cell)
                if downCell and self.getValue(downCell) == value:
                    count += 1
                else:
                    count = 0
                if count == minWidth:
                    return True
            count = 0
            for cell in column:
                upCell = self.getUpCell(cell)
                if upCell and self.getValue(upCell) == value:
                    count += 1
                else:
                    count = 0
                if count == minWidth:
                    return True
            return False
        for row in rows:
            if not isRowConnected(row):
                return False
        for column in columns:
            if not isRowConnected(row):
                return False
        # If everything was fine return true
        return True

    # Try to change a cell 'value' respecting all clusters
    # The value of the finally lost cell is the source
    # If the claimed cell has the source type it is converted and we are done
    # If not, try to claim other random source cell for the displaced cluster
    # The 'previous' variable tracks all previously displaced clusters to avoid loops
    def claimCell (self, cell, value, source, previous = []):
        self.track('Claiming cell ' + str(cell))
        targetValue = self.getValue(cell)
        if targetValue == value:
            self.track('FAIL: The cell already belongs to the specified cluster')
            return False
        if targetValue == source:
            self.track('SUCCESS: The cell has the source value')
            self.setValue(cell, value)
            print('easy claim: ' + str(cell))
            if (self.updater and len(previous) == 0):
                self.updater(self.format())
            return True
        self.track('The cell belongs to cluster ' + str(targetValue))
        # Add the value to the current list of previous clusters if not there 
        if (value not in previous):
            previous.append(value)
        # Expand the cluster whose cell we just claimed, to compensate
        if self.expandCluster(self.getCluster(targetValue), source, banned = previous):
            self.track('SUCCESS: Cell ' + str(cell) + ' was claimed')
            self.setValue(cell, value)
            print('claim: ' + str(cell))
            if (self.updater and len(previous) == 0):
                self.updater(self.format())
            return True
        else:
            self.track('FAIL: Cell ' + str(cell) + ' was not claimed')
            return False

    # Claim all cells in a group at the same time. If any of the claims fails, none is claimed.
    # All matrix values are saved at the begining, so we can back up when a claim fails
    def claimCellGroup (self, cells, value, source, previous = []):
        print('claiming group...')
        targetValues = [self.getValue(cell) for cell in cells]
        # If any of them belongs to the clamed value return false
        if any([v == value for v in targetValues]):
            return False
        # If all of them are free just get them
        if all([v == source for v in targetValues]):
            for cell in cells:
                self.setValue(cell, value)
                print('easy claim group: ' + str(cell))
            if (self.updater and len(previous) == 0):
                self.updater(self.format())
            return True
        # If not, save curent values and proceed to claim them one by one
        backup = [v for v in self.values]
        for cell in cells:
            if self.claimCell(cell, value, source, previous) == False:
                self.values = backup
                return False
        if (self.updater and len(previous) == 0):
            self.updater(self.format())
        return True

    # Get a random area with cells from the specified value with the specified size
    # If there is not any area with all cells with the specified value return an area with as much as possible
    # You can specify a cell which must be included in all yield areas
    def priorizeArea (self, value, xwide, ywide, includeCell = None, connected = True):
        results = []
        # Iterate over all matrix cells
        for i in range(0, self.size):
            origin = self.getCell(i)
            originX = origin.x
            originY = origin.y
            # If the expected area excels the matrix we skip this area
            if originX + xwide >= self.xrange or originY + ywide >= self.yrange:
                continue
            # If there is an 'includeCell' and it is not included we skip this area
            if (includeCell and not
            (includeCell.x >= originX and includeCell.x < originX + xwide
            and includeCell.y >= originY and includeCell.y < originY + ywide)):
                continue
            # Count how many cells have the source value in each area
            cells = []
            matchingCells = 0
            for x in range(originX, originX + xwide):
                for y in range(originY, originY + ywide):
                    cell = Cell(x, y)
                    cells.append(cell)
                    v = self.getValue(cell)
                    if v == value:
                        matchingCells += 1
            # Check all cells are connected when it is required
            if connected:
                allConnected = True
                for c in range(1, len(cells)):
                    if not self.connected(cells[c - 1], cells[c]):
                        allConnected = False
                        break
                if not allConnected:
                    continue
            # Save the origin cell and the number of matching cells
            results.append((cells, matchingCells))
        # Order areas randomly and by matching cells
        random.shuffle(results)
        def orderer(i):
            return i[1]
        results.sort(key = orderer, reverse=True)
        # Return results 1 by 1
        # Return a getter function
        for r in results:
            yield r[0]

    # Expected input cells are a cluster collider cells
    # Find all possible groups in input cells and priorize the most suitable groups to be claimed
    # (e.g. groups with less non-source cells)
    # Yield the result groups in order 1 by 1
    def priorizeGroups (self, cells, cluster, source, forceMaximum = None):
        self.track('Priorizing groups from: ' + str(cells))
        # Set the value, minimum and maximum values
        # The maximum may be null
        value = cluster.value
        minWidth = 0
        if cluster.minWidth:
            minWidth = cluster.minWidth
        maxWidth = 0
        if cluster.maxWidth:
            maxWidth = cluster.maxWidth
        if forceMaximum:
            maxWidth = min(maxWidth, forceMaximum)
        # Find all possible groups
        groups = self.getCellGroups(cells, minWidth, maxWidth, collapse = True)
        # Reorder them randomly
        random.shuffle(groups)
        # Now priorize groups in border according to cluster settings
        priorizeBorder = cluster.priorizeBorder
        if priorizeBorder:
            def borderer(g):
                # This is a boolean multiplied by an integer
                order = self.inBorderGroup(g) * priorizeBorder
                #print(str(cluster.priorizeBorder) + ' -> ' + str(order))
                return order
            groups.sort(key=borderer, reverse=True)
        # Now priorize groups by length
        behaviour = cluster.behaviour
        if behaviour:
            def behaviorder(gp1):
                count = 0
                for gp2 in groups:
                    if gp1 == gp2:
                        continue
                    for cell in gp1:
                        if cell in gp2:
                            count += 1
                            break
                if behaviour == 1:
                    order = count
                elif behaviour == 2:
                    order = -count
                return order
            groups = sorted(groups, key = behaviorder, reverse=True)
        # Now priorize groups with as less non-source cells as posible
        def sourcer(g):
            order = 0
            for c in g:
                if self.getValue(c) == source:
                    order += 1
            return order
        groups.sort(key = sourcer, reverse=True)
        # Yield groups 1 by 1
        for g in groups:
            # Yield cell groups 1 by 1
            yield g
        # Finally, if all previous yields have failed we stop here
        # This also happens when the only available cells are not together to form groups
        self.track('OVER: No more posible groups')
        #yield None

    # Expected input cells are a cluster collider and suitable cells
    # Return an array of tuples with this format: (Cells, Priority)
    # Cells are classified and pointed according how they conserve the min/max limits
    # First, cells which contribute to fill the minimums
    # Second, free cells
    # Third, cells which make the cluster stay as a rectangle
    # Fourth, cells which make the cluster stay as an L
    # Fifth, cells which keep the minimums/maximums somehow
    # Finally, the rest
    def priorizeCells (self, cells, cluster, source):
        self.track('Priorizing cells from: ' + str(cells))
        # Set the minimum and maximum values
        # The maximum may be null
        value = cluster.value
        minWidth = 0
        if cluster.minWidth:
            minWidth = cluster.minWidth
        maxWidth = 0
        if cluster.maxWidth:
            maxWidth = cluster.maxWidth
        # Set how 'advisable' is each cell to be the next expanded cell
        results = []
        for c in cells:
            priority = self.getCellPrioirty(c, cluster, source)
            results.append((c, priority))
        # Group cells by priority
        priorities = list(set([pc[1] for pc in results]))
        groupedResults = []
        for p in priorities:
            pcells = []
            for r in results:
                if r[1] == p:
                    pcells.append(r[0])
            groupedResults.append((pcells, p))
        # Reorder the cell groups by priority (from highest to lowest)
        def orderer(i):
            return i[1]
        groupedResults.sort(key = orderer, reverse=True)
        # Define functions to sort in order cells and group cells
        # Order functions according to if cells/groups are in border
        def borderer(c):
            # This is a boolean multiplied by an integer
            order = self.inBorder(c) * cluster.priorizeBorder
            #print(str(cluster.priorizeBorder) + ' -> ' + str(order))
            return order
        # Start to generate and yield cells or cell groups in priority order
        # Try all possibile cells or cell groups in first priority group
        # If none works in the first group jump to the next priority group and so on
        self.track('Priority order: ' + str(groupedResults))
        for (priorityGroup, priority) in groupedResults:
            # In border priorization
            priorityGroup.sort(key = borderer, reverse=True)
            for c in priorityGroup:
                # Yield cells 1 by 1
                yield c
        # Finally, if all previous yields have failed we return all cells
        # This also happens when the only available cells are not together to form groups
        self.track('OVER: Returning all cells')
        #yield None
        
    # Expand a the specified cluster as many cells as the 'expansion' number
    # New cells are claimed from cluster colliding cells
    # Cells with values in the 'banned' list are discarded during the expansion
    # This prevents an expansion chain to claim cells from the original cluster
    # WARNING: The 'banned' array seems to be conserved through different calls
    # WARNING: The last value must be 'poped' before returning or it remains
    def expandCluster (self, cluster, source, expansion = 1, banned = []):
        self.track('Expanding cluster ' + str(cluster) + ' for ' + str(expansion) + ' cells', deepen = 1)
        self.track('The following clusters are banned for this expansion:' + str(banned))
        while(expansion > 0):
            # Get all cells colliding the target cluster
            colliders = self.getClusterColliderCells(cluster.value)
            # Filter out the banned cluster cells
            self.track('Getting suitable cells from: ' + str(colliders), deepen = 1)
            suitables = colliders
            # Filter those cells excluding each 'previous' values
            for p in banned:
                suitables = self.filterCells(suitables, p, False)
            # Return filtered cells. They may be 'None'
            self.track('Suitable results: ' + str(suitables), deepen = -1)
            # If there are no suitable cells to expand we are trapped
            if suitables == None:
                self.track('FAIL: Expansion is trapped', deepen = -1)
                if(len(banned) > 0):
                    banned.pop()
                return False
            # Get priorized cells
            # This must always return results
            PREpreviousExpansion = expansion
            # When the number of cells to expand is equal or greater than the minimum, we claim in groups
            if expansion >= cluster.minWidth:
                forcedMax = None
                if expansion < cluster.maxWidth:
                    forcedMax = expansion
                priorizedGroups = self.priorizeGroups(suitables, cluster, source, forcedMax)
                for expanders in priorizedGroups:
                    self.track('Next group to expand: ' + str(expanders))
                    # Save the expansion number before searching suitable cells
                    previousExpansion = expansion
                    # Claim the whole group in 1 single step
                    if self.claimCellGroup(expanders, cluster.value, source, previous = banned):
                        expansion -= len(expanders)
                        # If expansion is 0 we are done
                        if expansion == 0:
                            self.track('SUCCESS: Expansion completed', deepen = -1)
                            if(len(banned) > 0):
                                banned.pop()
                            return True
                    # If there is no successful expansion
                    if expansion == previousExpansion:
                        self.track('NEXT: The maximum priority group is not available')
                        continue
                    # If we expanded with the last priorized cells we must stop
                    # It is mandatory to find new piorized cells to keep expanding
                    else:
                        self.track('OK: The maximum priority group was available')
                        break
            # When the number of cells to expand is lower than the minimum, we claim cells 1 by 1
            else:
                priorizedCells = self.priorizeCells(suitables, cluster, source)

                for cell in priorizedCells:
                    self.track('Next cell to expand: ' + str(cell))
                    # Save the expansion number before searching suitable cells
                    previousExpansion = expansion
                    # Claim each cell individually
                    if self.claimCell(cell, cluster.value, source, previous = banned):
                        expansion -= 1
                        # If expansion is 0 we are done
                        if expansion == 0:
                            self.track('SUCCESS: Expansion completed', deepen = -1)
                            if(len(banned) > 0):
                                banned.pop()
                            return True
                        # If we expanded with the last priorized cells we must stop
                        # It is mandatory to find new piorized cells to keep expandinge
                        self.track('OK: The maximum priority group was available')
                        break
                    # If there is no successful expansion
                    else:
                        self.track('NEXT: The maximum priority group is not available')
                        continue

            # If there is no successful expansion it means we are trapped
            if expansion == PREpreviousExpansion:
                print('Error: Cluster ' + str(cluster.value) + ' expansion has failed')
                self.track('FAIL: Last iteration failed to claim any cell', deepen = -1)
                if(len(banned) > 0):
                    banned.pop()
                return False
                
    # Return the cluster with the specified value
    # Return False if there is no cluster with the specified value
    def getCluster (self, value):
        for cluster in self.clusters:
            if cluster.value == value:
                return cluster
        return None
    
    # Functions to manage the matrix
    # Get the cell array index from the x and y matrix coordenates
    def cell2index (self, cell):
        return int(cell.y * self.xrange + cell.x)
    
    # Get the x and y matrix coordenates from a cell array index
    def index2cell (self, index):
        y = math.floor(index / self.xrange)
        x = index - y * self.xrange
        return Cell(x,y)
    
    # Get the real position equivalent to the center of a matrix cell
    def cell2point (self, cell):
        x = self.minx + cell.x * self.cellSize + self.cellSize/2
        y = self.miny + cell.y * self.cellSize + self.cellSize/2
        return Point(x,y)
    
    # Set the value of an index (int) transformed from any supported format
    # Return True if anything was fine and False if were problems
    # WARNING, there is no way to check if the input is already a value, so it will be transformed anyway
    def setValue (self, input, value):
        self.values[self.getIndex(input)] = value
        
    # Transform any supported format to value (int)
    # WARNING, there is no way to check if the input is already a value, so it will be transformed anyway
    def getValue (self, input):
        index = self.getIndex(input)
        if index >= self.size or index < 0:
            return None
        return self.values[self.getIndex(input)]
    # Transform any supported format to index (int)
    def getIndex (self, input):
        if(type(input) == Cell):
            return self.cell2index(input)
        elif(type(input) == int):
            return input
        else:
            return print('getIndex Error: not supported type -> ' + str(type(input)))
        
    # Transform any supported format to cell (int)
    def getCell (self, input):
        if(type(input) == int):
            return self.index2cell(input)
        elif(type(input) == Cell):
            return input
        else:
            return print('getCell Error: not supported type -> ' + str(type(input)))
        
    # Transform any supported format to Point
    def getPoint (self, input):
        if(type(input) == int):
            return self.cell2point(self.index2cell(input))
        elif(type(input) == Cell):
            return self.cell2point(input)
        elif(type(input) == Point):
            return input
        else:
            return print('getPoint Error: not supported type -> ' + str(type(input)))
        
    # Used for interesection calcules
    # This function was copied from stackoverflow
    def ccw (self,A,B,C):
        return (C.y-A.y) * (B.x-A.x) > (B.y-A.y) * (C.x-A.x)
    
    # Check if two cells are connected (true) or they are separated by a line (false)
    # Accept all formats: indexes, cells and points
    def connected (self, p1, p2):
        c = self.getPoint(p1)
        d = self.getPoint(p2)
        for l in self.lines:
            #print(str(Line(c,d)) + ' -> ' + str(l))
            if self.ccw(c,l.a,l.b) != self.ccw(d,l.a,l.b) and self.ccw(c,d,l.a) != self.ccw(c,d,l.b):
                return False
        return True
    
    # Display in console the logic progress
    # The level defines how important a message is
    # Only messages with levels lower than the 'tracked' value will be displayed
    # The 'deepen' index makes firther messages display after more 'tabulators'
    # A message can be 'forced' to be displayed even when the matrix is not tracked
    deep = 0
    def track (self, message, level = 0, deepen = 0, forced = False):
        previous = ''
        for d in range(self.deep):
            previous = previous + '  '
        if self.tracked > level or forced:
            print(previous + message)
        self.deep += deepen

    # Format the values list as a matrix (i.e. an array of arrays) 
    def format (self):
        formated = []
        rowSize = self.xrange
        rowCount = int(self.size/rowSize)
        for row in range(rowCount):
            firstIndex = row * rowSize
            lastIndex = (row + 1) * rowSize
            rowIndexes = self.values[ firstIndex : lastIndex ]
            formated.append(rowIndexes)
        return formated
        
    # Create a new mesh for debugging
    def represent (self):
        print('hi :)')