import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from collections import deque
import sys
import random

def plot_grid_time(dir, T):
    fig = plt.figure()

    while True:
        plt.ion()

        for t in T:
            file = dir + '/t' + str(t) + '.png'
            plt.clf()
            image = img.imread(file)
            plt.imshow(image)
            plt.axis('off')
            plt.title('t=%i' % (t), fontsize=15)
            plt.pause(0.005)
            plt.draw()


'''
# Task 2:
T = np.arange(0,2300,20)
#dir = 'plots4' #eta 2pi
#dir = 'plots3' #eta 0.2pi
dir = 'plots1' #eta 0.002pi
plot_grid_time(dir,T)
'''


# Task 3:
T = np.arange(0,3400,20)
#dir = 'plots_path_6' #eta 2pi
#dir = 'plots_path_5' #eta .2pi
dir = 'plots_path_4' #eta .002pi
plot_grid_time(dir,T)
