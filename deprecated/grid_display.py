#from matplotlib import use 
#use('AGG')
#import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import animation
from matplotlib.widgets import Slider 

import numpy as np
from multiprocessing import Process, Queue

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

# Set a list with all matrix value at each recorded step
frames = []
# Set a queue for the frames, since they are passed to a process
queue = Queue()

# Updater called from the matrix
def addFrame (data):
    print(' [ frame ' + str(len(frames)) + ' ] ')
    frames.append(data)
    queue.put(frames)

# Show the heatmap
def represent (queue):
    # Setup
    fig = plt.figure()
    frames = queue.get()
    matrice = plt.imshow(frames[0], vmax=8)
    ref = plt.gca()
    ref.invert_yaxis()
    # Get data x and y ranges
    x_range = len(frames[0][0])
    y_range = len(frames[0])
    # Gridlines
    ref.set_xticks(np.arange(-.5, x_range, 1), minor=True)
    ref.set_yticks(np.arange(-.5, y_range, 1), minor=True)
    ref.grid(which='minor', color='w', linestyle='-', linewidth=1)

    # Slider
    axslider = plt.axes([0.25, .03, 0.50, 0.02])
    slider = Slider(axslider, 'Frame', 0, len(frames), valinit=len(frames), valstep=1)

    # Animation updater
    def updateFrame (i):
        global frames
        # Update frames when the queue is not empty
        if queue.qsize() > 0:
            frames = queue.get()
        maximum = len(frames) - 1
        updated = slider.val == slider.valmax
        # Update the maximum slider value
        slider.valmax = maximum
        slider.ax.set_xlim(slider.valmin,slider.valmax) # This is necessary to make the valmax stable
        # If the slider is in the maximum value we keep it updated
        if updated:
            slider.set_val(maximum)
        matrice.set_data(frames[slider.val])
        
    # Run the animation and show the plot
    anim = animation.FuncAnimation(fig, updateFrame)
    plt.show()

# Set the heatmap representation in a paralel process so the matrix can keep beeing calculated
def setupRepresentProcess ():
    queue.put(frames)
    p = Process(target=represent, args=(queue, ))
    p.start()