import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider 

import numpy as np
from multiprocessing import Process, Queue

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
            slider.ax.set_xlim(slider.valmin, slider.valmax) # This is necessary to make the valmax stable
        # If the slider is in the maximum value we keep it updated
        if updated:
            slider.set_val(maximum)

        # Clear previous lines
        ax.lines = []
        ax.texts = []

        # Get everything to be displayed in the current frame
        data = frames[int(slider.val)]

        # Draw all lines
        lines = get_lines_from_anything(data)
        for line in lines:
            xs = [line.a.x, line.b.x]
            ys = [line.a.y, line.b.y]
            ploted_lines = ax.plot(xs,ys,color=line.color)

        # Draw all room names
        rooms = [ room for room in data if hasattr(room, 'name') ]
        for room in rooms:
            # Find the most wide maximum free rectangle
            # We focus in the x axis since text is horizontal
            # DANI: Esta es la linea que está generando tanta mierda
            place = room.text_place
            ploted_text = ax.text(place.x, place.y, room.name, color='black')
            #ploted_text = ax.annotate(room.name, (pmin.x, pmin.y), color='black')
        
    # Run the animation and show the plot
    anim = animation.FuncAnimation(fig, update_frame)
    plt.show()

# Set the heatmap representation in a paralel process so the matrix can keep beeing calculated
def setup_display ():
    queue.put(frames)
    p = Process(target=represent, args=(queue, ))
    p.start()

# --------------------------------------------------------------------------------------------------

# Plot lines, rectangles and perimeters manually
def get_lines_from_anything (things : list):
    lines = []
    for thing in things:
        # If it is a line or something with a and b (i.e. something "linealizable")
        if hasattr(thing, 'a') and hasattr(thing, 'b'):
            lines.append(thing)
        # If it is a rectangle, perimeter, or something with lines
        if hasattr(thing, 'lines'):
            lines += thing.lines
        # If it is a rectangle or something with a "crossing line" getter
        if hasattr(thing, 'get_crossing_line'):
            lines.append(thing.get_crossing_line())
        # If it is a room or something with a perimeter
        if hasattr(thing, 'perimeter'):
            if thing.perimeter:
                lines += thing.perimeter.lines
    return lines