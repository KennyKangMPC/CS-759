#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 13:55:34 2020

@author: Kangqi Fu
"""

from numpy import loadtxt, reshape
from pylab import ioff
import matplotlib.pyplot as plt
from glob import glob
import os

ioff()

fileNames = glob("./output/Solution*.dat")
fileNames.sort()

for fileName in fileNames:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    f = open(fileName, "r")
    xCells = int(f.readline().split(":")[1])
    yCells = int(f.readline().split(":")[1])
    numGhostCells = int(f.readline().split(":")[1])
    time = float(f.readline().split(":")[1])
    cfl = float(f.readline().split(":")[1])
    f.close()
    
    x, y, u = loadtxt(fileName, skiprows = 5, unpack=True)
    x = reshape(x, (xCells + 2 * numGhostCells, yCells + 2 * numGhostCells))
    y = reshape(y, (xCells + 2 * numGhostCells, yCells + 2 * numGhostCells))
    u = reshape(u, (xCells + 2 * numGhostCells, yCells + 2 * numGhostCells))
    #ax.set_xlim(-1.5, 1.5)
    #ax.set_ylim(-0.1, 1.1)
    plt.contourf(x, y, u, 100, cmap='jet')
    #plt.contourf(x, y, u,100, cmap='ocean_r')
    
    plt.colorbar()
    ax.set_title("CFL = %5.2f"%cfl + ", Times = %5.3f"%time)
    fig.savefig(fileName.replace(".dat", ".png"))

os.system("eog " + fileNames[0].replace(".dat",".png"))
