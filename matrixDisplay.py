import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import animation

from multiprocessing import Process, Array, Value

# Set the heatmap colors
colors = ListedColormap([
    'black',
    'white',
    'blue',
    'green',
    'yellow',
    'red',
    'pink',
])

frames = []
#frames = Array('i',0)
#frames = Value('list',0)
def addFrame (data):
    frames.append(data)
    #if matrice:
    #    matrice.set_array(data)

# Show the heatmap
def represent ():
    fig = plt.figure()
    matrice = plt.imshow(frames[0], vmax=8)
    plt.gca().invert_yaxis()
    print('Representing ' + str(len(frames)) + ' frames')

    def updateFrame (i):
        matrice.set_data(frames[i])

    anim = animation.FuncAnimation(fig, updateFrame, repeat = False)
    plt.show()

# Set the heatmap representation in a paralel process so the matrix can keep beeing calculated
def setupRepresentProcess ():
    p = Process(target=represent)
    p.start()