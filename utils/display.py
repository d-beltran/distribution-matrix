import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import to_rgba 
from matplotlib.widgets import Slider, Button
import matplotlib.patches as mpatches

import math
from multiprocessing import Process
from threading import Thread

from typing import Optional, List

# Remove the slider silly warning
import warnings
warnings.filterwarnings("ignore")

# Import global variables
from utils.auxiliar import GLOBAL

# Set a list with all system values at each recorded step
global_frames = []

# Updater called from the system
def add_frame (data : list, title : Optional[str] = None):
    # If the display is not enabled then there is nothing to do here
    # Note that this check is also done in the room update_display function
    # However there are many debug calls to this function which do not go through it
    # WARNING: Otherwise we would fill the queue with frames which nobody is ever going to read
    # WARNING: This would also keep the process alive forever, since the queue is never emptied
    if not GLOBAL['enabled_display']: return
    # First check if we reached the limkit
    # If so stop here
    if GLOBAL['frame_count'] > GLOBAL['frames_limit']: raise SystemExit('Reached displayed frames limit')
    display_message = title if title else 'No title'
    print(f' [ frame {len(global_frames)} ] - {display_message}')
    if type(data) != list:
        data = [data]
    # Remove duplicates
    data = list(set(data))
    # Mine displayable segments from data
    segments = get_segments_from_anything(data)
    # This is for the fillings only, to add color
    rects = get_rects_from_anything(data)
    traced = [ element for element in data if hasattr(element, 'name') ]
    frame = (segments, rects, traced, display_message)
    global_frames.append(frame)
    GLOBAL['frame_count'] += 1
    # Send only the new frame to the display process, which keeps its own list of frames
    # Note that a frame is never modified after this point, so it does not matter when the queue pickles it
    frames_queue = GLOBAL['frames_queue']
    frames_queue.put(frame)

# Show the heatmap
# The frames queue must be passed as argument
# WARNING: Do not get it from GLOBAL, since the child process may not inherit the parent memory (e.g. forkserver)
# WARNING: In such case the child would import a new and disconnected queue
def represent (frames_queue):

    # Setup
    # Use the constrained layout so the legend (above) and the title (below) never overlap the image
    # Reserve the bottom of the figure for the slider and buttons, which are not managed by the layout
    # Make the window bigger than the default so tall or wide drawings are not too small
    fig, ax = plt.subplots(layout='constrained', figsize=(12, 8.5))
    fig.get_layout_engine().set(rect=(0, 0.05, 1, 0.95))
    # Set the list of frames received so far
    frames = []
    # Keep reading frames from the queue in a background thread
    # Thus the queue is always emptied, even while the display is busy drawing
    # WARNING: Frames are heavy and the queue pipe is small, so reading only between draws is too slow
    def read_frames ():
        while True:
            frames.append(frames_queue.get())
    Thread(target=read_frames, daemon=True).start()

    # Remove top and right box segments
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Slider
    axslider = plt.axes([0.25, .02, 0.50, 0.02])
    slider = Slider(axslider, label='Frame', valmin=0, valmax=len(frames), valinit=len(frames), valstep=1, valfmt='%0.0f')

    # Functions for when the user clicks in the previous/next arrows
    def previous_frame (event):
        slider.set_val(slider.val - 1)
    def next_frame (event):
        slider.set_val(slider.val + 1)

    axprev = plt.axes([0.86, 0.01, 0.06, 0.04])
    axnext = plt.axes([0.92, 0.01, 0.06, 0.04])
    bprev = Button(axprev, '<-')
    bprev.on_clicked(previous_frame)
    bnext = Button(axnext, '->')
    bnext.on_clicked(next_frame)

    # Add the title
    # Use the x label as title, so it is always placed below the tick labels
    ax.set_xlabel('(empty)')

    # Track any time the previous slider value
    previous_slider_value = None

    # Animation updater
    def update_frame (i):
        nonlocal previous_slider_value
        # Set the slider range
        minimum = 0
        maximum = len(frames) - 1
        # Update the slider val in case it is out of minimum/maximum range
        # This may happen when using the buttons
        if slider.val < minimum:
            slider.set_val(minimum)
        if slider.val > maximum:
            slider.set_val(maximum)
        # Find out if the slider value is the last (current) value
        updated = slider.val == slider.valmax
        # Update the maximum slider value
        slider.valmax = maximum
        if maximum != 0:
            slider.ax.set_xlim(slider.valmin, slider.valmax) # This is necessary to make the valmax stable
        # If the slider is in the maximum value we keep it updated
        if updated:
            slider.set_val(maximum)

        # Check the current slider value
        # If it has no changed since last update then do nothing
        slider_value = int(slider.val)
        if slider_value == previous_slider_value:
            return
        # Update the tracked slider value
        previous_slider_value = slider_value

        # Clear previous segments and rects
        ax.clear()

        # Make axes respect the 1:1 ration
        # Otherwise the whole image is deformed in case there is a dmension larger than the other
        ax.axes.set_aspect('equal')
        # Stick the image to the bottom so the title is close to the frames slider
        ax.set_anchor('S')

        # In case something went wrong in the solving process and there are no frames from the begining
        # Prevent error logs from the display to be shown in the console
        # But do not close the display window (i.e. do not call raise SystemExit)
        # It is better for the developer to see an empty display to understand that it failed before returning any frame
        if len(frames) == 0:
            return

        # Get everything to be displayed in the current frame
        segments, rects, traced, display_message = frames[slider_value]

        # Draw all segments
        for segment in segments:
            xs = [segment.a.x, segment.b.x]
            ys = [segment.a.y, segment.b.y]
            zs = segment.z if hasattr(segment, 'z') else None
            ploted_segments = ax.plot(xs, ys, color=segment.color, zorder=zs)

        # Draw all rect areas
        for rect in rects:
            xs = [rect.x_min, rect.x_max, rect.x_max, rect.x_min]
            ys = [rect.y_min, rect.y_min, rect.y_max, rect.y_max]
            # WARNING: Use 'facecolor' instead of 'color' to hide separation segments between fills
            ploted_rects = ax.fill(xs, ys, facecolor=rect.fill_color or 'C0', alpha=rect.opacity)
            # Apply textures for those rects which have one
            if rect.texture:
                ploted_rects[0].set_hatch(rect.texture)

        # Set the legend with all room names
        handles = []
        for track in traced:
            facecolor = to_rgba(track.fill_color)
            facecolor = (facecolor[0], facecolor[1], facecolor[2], 0.2) # Reduce the opacity
            patch = mpatches.Patch(facecolor=facecolor, edgecolor=track.segments_color, label=track.name)
            handles.append(patch)
        columns_number = math.ceil( len(handles) / 2 )
        if columns_number > 0:
            # Anchor the legend bottom just above the axes top, so it grows upwards and never overlaps the image
            ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=columns_number)
        #legend.handles = handles

        # Add the title
        ax.set_xlabel(display_message)
        
    # Run the animation and show the plot
    anim = animation.FuncAnimation(fig, update_frame)
    plt.show()

# Set the heatmap representation in a paralel process so the matrix can keep beeing calculated
def setup_display (frames_limit : Optional[int] = None):
    # Set the global display option as true
    global GLOBAL
    GLOBAL['enabled_display'] = True
    if frames_limit != None:
        GLOBAL['frames_limit'] = frames_limit
    # Start the display logic
    # Note that the queue is passed as argument so both processes share the same queue
    p = Process(target=represent, args=(GLOBAL['frames_queue'], ))
    p.start()

# --------------------------------------------------------------------------------------------------

# Mine all possible segments from a list of different vectorial_base elements
def get_segments_from_anything (things : list) -> List['Segment']:
    total_segments = []
    for thing in things:
        thing_segments = []
        # If it is a segment or something with a and b (i.e. something "segmentalizable")
        if hasattr(thing, 'a') and hasattr(thing, 'b'):
            thing_segments.append(thing)
        # If it is a rectangle, polygon, or something with segments
        if hasattr(thing, 'segments'):
            thing_segments += thing.segments
        # If it is a rectangle or something with a "crossing segment" getter
        if hasattr(thing, 'get_crossing_segment'):
            thing_segments.append(thing.get_crossing_segment())
        # If it is a grid or something with rects
        if hasattr(thing, 'rects'):
            for rect in thing.rects:
                thing_segments += rect.segments
                thing_segments.append(rect.get_crossing_segment())
        # If it is a boundary or something with a polygon
        if hasattr(thing, 'polygon'):
            if thing.polygon:
                thing_segments += thing.polygon.segments
        # If it is a room or something with a boundary
        if hasattr(thing, 'boundary'):
            boundary = thing.boundary
            if boundary:
                thing_segments += boundary.exterior_polygon.segments
                for polygon in boundary.interior_polygons:
                    thing_segments += polygon.segments
        # If it has a corridor grid (i.e. it is a room)
        if hasattr(thing, 'corridor_grid'):
            if thing.corridor_grid != None:
                corridor_segments = sum([ boundary.segments for boundary in thing.corridor_grid.boundaries ], [])
                thing_segments += corridor_segments
        # If it is a room or something with doors
        if hasattr(thing, 'doors'):
            doors = thing.doors or []
            for door in doors:
                segment = door.segment
                if not segment:
                    continue
                segment.color = 'white'
                segment.z = 17 # Make this segment display in the top layer
                thing_segments.append(segment)
                # Create a new segment to represent the door open
                open_door_segment = door.get_open_door()
                thing_segments.append(open_door_segment)
        # If it is a room or something with windows
        if hasattr(thing, 'windows'):
            windows = thing.windows or []
            for window in windows:
                segment = window.segment
                if not segment:
                    continue
                segment.color = 'grey'
                segment.z = 17 # Make this segment display in the top layer
                thing_segments.append(segment)
        # Make sure the new segments inherit the z property, if they have it
        if hasattr(thing, 'z'):
            for segment in thing_segments:
                # If a segment already has a z attribute then do not overwrite it
                if hasattr(segment, 'z'): continue
                segment.z = thing.z
        # Make sure the new segments inherit the color, if they have it
        if hasattr(thing, 'color'):
            for segment in thing_segments:
                segment.color = thing.color
        total_segments += thing_segments
    return total_segments

# Mine all possible rectangles from a list of different vectorial_base elements
def get_rects_from_anything (things : list) -> List['Rect']:
    rects = []
    for thing in things:
        # If it is a rectangle or something with x_min and y_min
        if hasattr(thing, 'x_min') and hasattr(thing, 'y_min'):
            rects.append(thing)
        # If it is a grid
        if hasattr(thing, 'rects'):
            rects += thing.rects
        # If it is a boundary
        # if hasattr(thing, 'grid'):
        #     if thing.grid:
        #         rects += thing.grid.rects
        # If it has a discarded grid (i.e. it is a room)
        if hasattr(thing, 'discarded_grid'):
            if thing.discarded_grid != None:
                discareded_rects = thing.discarded_grid.rects
                for rect in discareded_rects:
                    rect.fill_color = thing.fill_color
                    rect.texture = '///'
                rects += discareded_rects
        # If it has a corridor grid (i.e. it is a room)
        if hasattr(thing, 'corridor_grid'):
            if thing.corridor_grid != None:
                corridor_rects = thing.corridor_grid.rects
                for rect in corridor_rects:
                    rect.texture = '--'
                rects += corridor_rects
        # If it has a free grid (i.e. it is a room)
        # This may fail for a parent free grid in steps where children overlap
        try:
            if hasattr(thing, 'free_grid'):
                # WARNING: Do not put any code below this part or it may not be run in some frames
                # This part is prone to fail
                if thing.boundary:
                    free_grid_rects = thing.free_grid.rects
                    for rect in free_grid_rects:
                        rect.segments_color = thing.segments_color
                        rect.fill_color = thing.fill_color
                    rects += free_grid_rects
                # WARNING: Do not write code here!!
        except:
            pass
    # Remove duplicates
    return rects