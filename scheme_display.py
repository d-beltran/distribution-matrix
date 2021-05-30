import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.widgets import Slider 
import matplotlib.patches as mpatches

import math
import numpy as np
from multiprocessing import Process, Queue

# Set a list with all system values at each recorded step
frames = []
# Set a queue for the frames, since they are passed to a process
queue = Queue()

# Updater called from the system
def add_frame (data):
    print(' [ frame ' + str(len(frames)) + ' ] ')
    segments = get_segments_from_anything(data)
    rects = get_rects_from_anything(data)
    traced = [ element for element in data if hasattr(element, 'name') ]
    frames.append((segments, rects, traced))
    queue.put(frames)

# Show the heatmap
def represent (queue):

    # Setup
    fig, ax = plt.subplots()
    frames = queue.get()

    # Remove top and right box segments
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

        # Clear previous segments and rects
        #ax.segments = []
        ax.clear()

        # Get everything to be displayed in the current frame
        segments, rects, traced = frames[int(slider.val)]

        # Draw all segments
        for segment in segments:
            xs = [segment.a.x, segment.b.x]
            ys = [segment.a.y, segment.b.y]
            ploted_segments = ax.plot(xs,ys,color=segment.color)

        # Draw all rect areas
        for rect in rects:
            pulc = rect.get_upper_left_corner()
            pbrc = rect.get_bottom_right_corner()
            xs = [rect.pmin.x, pulc.x, rect.pmax.x, pbrc.x]
            ys = [rect.pmin.y, pulc.y, rect.pmax.y, pbrc.y]
            # WARNING: Use 'facecolor' instead of 'color' to hide separation segments between fills
            ploted_rects = ax.fill(xs, ys, facecolor=rect.fill_color or 'C0', alpha=0.2)

        # Set the legend with all room names
        handles = []
        for track in traced:
            patch = mpatches.Patch(color=track.fill_color, label=track.name)
            handles.append(patch)
        columns_number = math.ceil( len(handles) / 2 )
        ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=columns_number)
        #legend.handles = handles
        
    # Run the animation and show the plot
    anim = animation.FuncAnimation(fig, update_frame)
    plt.show()

# Set the heatmap representation in a paralel process so the matrix can keep beeing calculated
def setup_display ():
    queue.put(frames)
    p = Process(target=represent, args=(queue, ))
    p.start()

# --------------------------------------------------------------------------------------------------

# Mine all possible segments from a list of different vectorial_base elements
def get_segments_from_anything (things : list):
    segments = []
    for thing in things:
        # If it is a segment or something with a and b (i.e. something "segmentalizable")
        if hasattr(thing, 'a') and hasattr(thing, 'b'):
            segments.append(thing)
        # If it is a rectangle, perimeter, or something with segments
        if hasattr(thing, 'segments'):
            segments += thing.segments
        # If it is a rectangle or something with a "crossing segment" getter
        if hasattr(thing, 'get_crossing_segment'):
            segments.append(thing.get_crossing_segment())
        # If it is a room or something with a perimeter
        if hasattr(thing, 'perimeter'):
            if thing.perimeter:
                segments += thing.perimeter.segments
    return segments

# Mine all possible rectangles from a list of different vectorial_base elements
def get_rects_from_anything (things : list):
    rects = []
    for thing in things:
        # If it is a rectangle or something with pmin and pmax (i.e. something "rectanglizable")
        if hasattr(thing, 'pmin') and hasattr(thing, 'pmax'):
            rects.append(thing)
        # If it is a perimeter
        if hasattr(thing, 'rects'):
            rects += thing.rects
        # If it is a room
        if hasattr(thing, 'free_rects'):
            if thing.free_rects:
                rects += thing.free_rects
    return rects