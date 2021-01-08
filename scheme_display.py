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

# Set a list with all system values at each recorded step
frames = []
# Set a queue for the frames, since they are passed to a process
queue = Queue()

# Updater called from the system
def add_frame (data):
    print(' [ frame ' + str(len(frames)) + ' ] ')
    frames.append(data)
    queue.put(frames)

# Show the heatmap
def represent (queue):

    # Setup
    #fig = plt.figure()
    fig, ax = plt.subplots()
    frames = queue.get()

    # Remove top and right box lines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Slider
    axslider = plt.axes([0.25, .03, 0.50, 0.02])
    slider = Slider(axslider, 'Frame', 0, len(frames), valinit=len(frames), valstep=1)

    # Animation updater
    def update_frame (i):
        global frames
        # Update frames when the queue is not empty
        if queue.qsize() > 0:
            frames = queue.get()
        maximum = len(frames) - 1
        updated = slider.val == slider.valmax
        # Update the maximum slider value
        slider.valmax = maximum
        if maximum != 0:
            slider.ax.set_xlim(slider.valmin,slider.valmax) # This is necessary to make the valmax stable
        # If the slider is in the maximum value we keep it updated
        if updated:
            slider.set_val(maximum)

        # Clear previous lines
        ax.lines = []

        # Draw all lines
        lines = frames[slider.val]
        for line in lines:
            xs = [line.a.x, line.b.x]
            ys = [line.a.y, line.b.y]
            ploted_lines = ax.plot(xs,ys,color='black')
        
    # Run the animation and show the plot
    anim = animation.FuncAnimation(fig, update_frame)
    plt.show()

# Set the heatmap representation in a paralel process so the matrix can keep beeing calculated
def setup_display ():
    queue.put(frames)
    p = Process(target=represent, args=(queue, ))
    p.start()

# --------------------------------------------------------------------------------------------------

# Plot lines manually

def plot_lines (lines : list):
    fig, ax = plt.subplots()
    # Remove top and right box lines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for line in lines:
        xs = [line.a.x, line.b.x]
        ys = [line.a.y, line.b.y]
        ax.plot(xs,ys,color='black')
    plt.show()