import seaborn as sb
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


colors = ListedColormap([
    'black',
    'white',
    'blue',
    'green',
    'yellow',
    'red',
    'pink',
])

def represent (data):
    mr = sb.heatmap(data, linewidths=.5, cmap=colors)
    mr.invert_yaxis()
    plt.show()