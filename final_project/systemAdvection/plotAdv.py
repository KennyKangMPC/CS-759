#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 13:55:34 2020

@author: Kangqi Fu
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os

plt.ioff()

fileNames = glob("./output/Solution*.dat")
fileNames.sort()
fig = plt.figure()

for fileName in fileNames:    
    #ax = fig.add_subplot(111)
    ax = fig.gca(projection = '3d')
    f = open(fileName, "r")
    xCells = int(f.readline().split(":")[1])
    yCells = int(f.readline().split(":")[1])
    numGhostCells = int(f.readline().split(":")[1])
    time = float(f.readline().split(":")[1])
    cfl = float(f.readline().split(":")[1])
    f.close()
    
    x, y, h, hu, hv = np.loadtxt(fileName, skiprows = 5, unpack=True)
    x = np.reshape(x, (xCells + 2 * numGhostCells, yCells + 2 * numGhostCells))
    y = np.reshape(y, (xCells + 2 * numGhostCells, yCells + 2 * numGhostCells))
    h = np.reshape(h, (xCells + 2 * numGhostCells, yCells + 2 * numGhostCells))
    #ax.set_xlim(-1.5, 1.5)
    #ax.set_ylim(-0.1, 1.1)
    # contour plot
        
    plt.cla()
    plt.clf()
    plt.xticks([])
    plt.yticks([])
    plt.contourf(x, y, h, 100, cmap='jet') # yeah jet is much better
    plt.axes().set_aspect("equal")
    #plt.contourf(x, y, u,100, cmap='ocean_r')
    
    plt.title("CFL = %5.2f"%cfl + ", Times = %5.3f"%time)
    fig.savefig(fileName.replace(".dat", ".png"))
    plt.colorbar()
    
# =============================================================================
#     # Surface plot
#     ax.plot_surface(x[2:-2, 2:-2], y[2:-2, 2:-2], h[2:-2, 2:-2], rstride=2, cstride=2, linewidth=0.1)
#     tray = np.zeros_like(x)
#     tray[0, :] = 0.5; tray[:, 0] = 0.5; tray[-1,:] = 0.5; tray[:,- 1] = 0.5;
#     ax.plot_surface(x, y, tray, rstride=2, cstride=2, linewidth = 0, color=(0.5, 0.5, 0.5, 0.1))
#     ax.set_zlim(-0.01, 1.01)
#     fig.savefig(fileName.replace(".dat", ".png"))
#     plt.clf()
# =============================================================================
    
os.system("eog " + fileNames[0].replace(".dat",".png"))
