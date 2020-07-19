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
# min and max widths: ...
class Cluster:
    def __init__(self, value, size, minWidth, maxWidth = None, priorizeBorder = 0):
        self.value = value
        self.size = size
        self.minWidth = minWidth
        self.maxWidth = maxWidth
        self.priorizeBorder = priorizeBorder
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
        self.value = [0] * self.size
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
            for c in self.getColliderCells(cell, value=originValue, connected=True):
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
        # Save the new cluster in the clusters list
        self.clusters.append(cluster)
        # Set the first cluster cell
        if origin == None:
            first = self.claimCell(self.getRandomCell(list(range(self.size)), value = sourceValue), cluster.value, sourceValue)
        else:
            first = self.claimCell(origin, cluster.value, sourceValue)
        # If we didn't get the origin cell return error
        if first == False:
            self.track("FAIL: New cluster failed to get the origin cell", deepen = -1)
            return False
        # Try to expand the cluster cells to reach the cluster size
        # We substract the origin cell from the size count (cluster.size - 1)
        if self.expandCluster(cluster, sourceValue, expansion = cluster.size - 1) == False:
            self.track("FAIL: New cluster failed to expand", deepen = -1)
            return False
        # Finally set the cluster as set and return True
        cluster.set = True
        self.track("SUCCESS: New cluster " + str(cluster) + " was set successfully", deepen = -1)
        return True
    
    # Get the number of cells with the specified value
    def countCells (self, value):
        #print(sum(v == value for v in self.value))
        return sum(v == value for v in self.value)
    
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
    
    # Return cells which collide the specified cell
    # They can be filtered to be connected or from a specific value
    # Impossible cells (i.e. out of the matrix limits) are not returned
    def getColliderCells (self, cell, value = None, connected = False):
        self.track('Getting cells colliding ' + str(cell), level = 1, deepen = 1)
        # Get all cells around
        colliders = [Cell(cell.x - 1, cell.y),
                     Cell(cell.x + 1, cell.y),
                     Cell(cell.x, cell.y - 1),
                     Cell(cell.x, cell.y + 1),
                     Cell(cell.x - 1, cell.y - 1),
                     Cell(cell.x + 1, cell.y - 1),
                     Cell(cell.x - 1, cell.y + 1),
                     Cell(cell.x + 1, cell.y + 1)]
        result = []
        for c in colliders:
            if( c.x >= 0 and c.x < self.xrange
            and c.y >= 0 and c.y < self.yrange
            and (value == None or self.getValue(c) == value)
            and (connected == False or self.connected(cell,c))):
                result.append(c)
        self.track('Cell colliding Results: ' + str(result), level = 1, deepen = -1)
        return result
    
    # Check if the specified cell is in the matrix limit or next to any edge
    def inBorder (self, cell):
        colliders = self.getColliderCells(cell, connected = True)
        return len(colliders) < 8
    
    # Check if any of the specified cells is in the matrix limit or next to any edge
    def inBorderGroup (self, cells):
        for cell in cells:
            if self.inBorder(cell):
                return True
        return False
    
    # Get all cells colliding any cell of the specified cluster
    def getClusterColliderCells (self, clusterValue, value = None, connected = False):
        self.track('Getting cells colliding cluster ' + str(clusterValue), deepen = 1)
        result = []
        # Get the colliding cells of each cluster cells
        for c in self.filterCells(list(range(self.size)), clusterValue):
            for col in self.getColliderCells(c, connected = connected):
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
    # (belongs, -x. +x, -y, +y, -x-y, +x-y, -x+y, +x+y)
    def getCellRelatives (self, cell, value):
        # Check if the cell belongs to the specified value
        belongs = self.getValue(cell) == value
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
        return (belongs, nx, px, ny, py, nxny, pxny, nxpy, pxpy)
    
    # Get an array with all rows in a cell array
    # The 'connected' argument may be set as true to get only connected rows
    # If so, multiple arrays from the same row may be returned if it is interrupted
    # Cells in every row are sortered
    def getRows (self, cells, connected = False):
        rows = []
        # Keep track of the currently calculated cells to avoid repeating searches
        searchedRows = []
        # Sorting function
        def byX (elem):
            return elem.x
        for c in cells:
            if c not in searchedRows:
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
        return rows;
    
    # Get an array with all columns in a cell array
    # The 'connected' argument may be set as true to get only connected columns
    # If so, multiple arrays from the same column may be returned if it is interrupted
    # Cells in every column are sortered
    def getColumns (self, cells, connected = False):
        columns = []
        # Keep track of the currently calculated cells to avoid repeating searches
        searchedColumns = []
        # Sorting function
        def byY (elem):
            return elem.y
        for c in cells:
            if c not in searchedColumns:
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
        return columns;
    
    # Get an array with all possible cell groups in a cell array
    # Cell groups are multiple cells colliding in the same row/column
    # The number of cells in each group is defined by 'number'
    def getCellGroups (self, cells, number):
        # First of all, try to get complete rows/columns
        rows = self.getRows(cells, connected = True)
        columns = self.getColumns(cells, connected = True)
        # Get all possible groups in those rows and columns
        groups = []
        for r in rows:
            groupsNumber = len(r) - number + 1
            if(groupsNumber > 0):
                for g in range(groupsNumber):
                    newGroup = []
                    for i in range(number):
                        newGroup.append(r[g+i])
                    groups.append(newGroup)
        for c in columns:
            groupsNumber = len(c) - number + 1
            if(groupsNumber > 0):
                for g in range(groupsNumber):
                    newGroup = []
                    for i in range(number):
                        newGroup.append(c[g+i])
                    groups.append(newGroup)
        return(groups)
            
    # Find if the cluster is respecting each of the cluster rules
    # 1. Find if it has the minimum number of cells to fill the minimum widths
    def getClusterStatus (self, cluster):
        print('getting status...')
        clusterCells = self.filterCells(list(range(self.size)), cluster.value)
        return len(clusterCells) >= cluser.minWidth * cluser.minWidth
    
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
            if (self.updater):
                self.updater(self.format())
            self.setValue(cell, value)
            return True
        self.track('The cell belongs to cluster ' + str(targetValue))
        # Add the value to the current list of previous clusters if not there 
        if (value not in previous):
            previous.append(value)
        # Expand the cluster whose cell we just claimed, to compensate
        if self.expandCluster(self.getCluster(targetValue), source, banned = previous):
            self.track('SUCCESS: Cell ' + str(cell) + ' was claimed')
            if (self.updater):
                self.updater(self.format())
            self.setValue(cell, value)
            return True
        else:
            self.track('FAIL: Cell ' + str(cell) + ' was not claimed')
            return False
        
    # Expected input cells are a cluster collider and suitable cells
    # Return an array of tuples with this format: (Cells, Priority)
    # Cells are classified and pointed according how they conserve the min/max limits
    # First, cells which contribute to fill the minimums
    # Second, cells which make the cluster stay as a rectangle
    # Third, cells which make the cluster stay as an L
    # Fourth, cells which keep the minimums/maximums somehow
    # Fifth, the rest
    def priorizeFitCells (self, cells, cluster, source):
        self.track('Priorizing cells from: ' + str(cells))
        results = []
        # Set the minimum and maximum values
        # The maximum may be null
        minWidth = cluster.minWidth
        maxWidth = cluster.maxWidth
        # Set how 'advisable' is each cell to be the next expanded cell
        for c in cells:
            relatives = self.getCellRelatives(c, cluster.value)
            if relatives[0] != False:
                continue
            priority = 0
            # From 1 to 4
            for i in range(1,5):
                # -100 points when this side is exactly in the minimum limit
                if relatives[i] == minWidth:
                    priority -= 100
                # -10 points when this side is exactly in the maximum limit
                elif relatives[i] == maxWidth:
                    priority -= 10
                # 0 points when this side is 0
                # If the final points are 0 this cell will be discarded
                elif relatives[i] == 0:
                    continue
                # +1 points when this side has already crossed the maximum
                elif relatives[i] > maxWidth:
                    priority += 1
                # +1000 points when this side has not reached the minimum yet
                elif relatives[i] < minWidth:
                    priority += 1000
            # Degrade the priority of cells with priority 0
            if priority == 0:
                priority = -1000
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
        # Once inside the loop, reorder the groups/cells inside each priority group
        # They must reorder according to if cells/groups belong to the source
        def sourcer(c):
            order = (self.getValue(c) == source) * 1
            return order
        def groupSourcer(g):
            order = 0
            for c in g:
                if self.getValue(c) == source:
                    order += 1
            return order
        # They must reorder according to if cells/groups are in border
        def borderer(c):
            # This is a boolean multiplied by an integer
            order = self.inBorder(c) * cluster.priorizeBorder
            #print(str(cluster.priorizeBorder) + ' -> ' + str(order))
            return order
        def groupBorderer(g):
            # This is a boolean multiplied by an integer
            order = self.inBorderGroup(g) * cluster.priorizeBorder
            #print(str(cluster.priorizeBorder) + ' -> ' + str(order))
            return order
        # Start to generate and yield cells or cell groups in priority order
        # Try all possibile cells or cell groups in first priority group
        # If none works in the first group jump to the next priority group and so on
        self.track('Priority order: ' + str(groupedResults))
        for (group, priority) in groupedResults:
            # If the priorize function finds cells to complete the minimums
            random.shuffle(group)
            if(priority > 4):
                # In border priorization
                group.sort(key = sourcer, reverse=True)
                group.sort(key = borderer, reverse=True)
                for c in group:
                    # Yield cells 1 by 1
                    yield ([c], False)
            # Otherwise
            else:
                # They must be claimed in groups
                groups = self.getCellGroups(group, minWidth)
                random.shuffle(groups)
                groups.sort(key = groupSourcer, reverse=True)
                groups.sort(key = groupBorderer, reverse=True)
                for g in groups:
                    # Yield cell groups 1 by 1
                    yield (g, True);
            # Finally, if all previous yields have failed we return all cells
            # This also happens when the only available cells are not together to form groups
            self.track('OVER: Returning all cells')
            yield (cells, False)
        
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
            colliders = self.getClusterColliderCells(cluster.value, connected = True)
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
            priorizedGroups = self.priorizeFitCells(suitables, cluster, source)
            for (expanders, isGroup) in priorizedGroups:
                self.track('Next cells to expand: ' + str(expanders))
                print(expanders)
                print(isGroup)
                # Save the expansion number before searching suitable cells
                previousExpansion = expansion
                # Claim as many cells as needed
                for cell in expanders:
                    if self.claimCell(cell, cluster.value, source, previous = banned):
                        expansion -= 1
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
            # If there is no successful expansion it means we are trapped
            if expansion == PREpreviousExpansion:
                print('Error: Cluster '+str(cluster.value)+' expansion has failed')
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
        self.value[self.getIndex(input)] = value
        
    # Transform any supported format to value (int)
    # WARNING, there is no way to check if the input is already a value, so it will be transformed anyway
    def getValue (self, input):
        index = self.getIndex(input)
        if index >= self.size or index < 0:
            return None
        return self.value[self.getIndex(input)]
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
            rowIndexes = self.value[ firstIndex : lastIndex ]
            formated.append(rowIndexes)
        return formated
        
    # Create a new mesh for debugging
    def represent (self):
        print('hi :)')