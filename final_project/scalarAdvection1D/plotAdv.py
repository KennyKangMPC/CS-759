#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 13:55:34 2020

@author: kenny
"""

from numpy import loadtxt
from pylab import figure, ioff
from glob import glob
import os

ioff()

fileNames = glob("./output/Solution*.dat")
fileNames.sort()

for fileName in fileNames:
    fig = figure()
    ax = fig.add_subplot(111)
    f = open(fileName, "r")
    xCells = int(f.readline().split(":")[1])
    numGhostCells = int(f.readline().split(":")[1])
    time = float(f.readline().split(":")[1])
    cfl = float(f.readline().split(":")[1])
    f.close()
    
    x, u = loadtxt(fileName, skiprows = 4, unpack=True)
    #x, u = loadtxt(fileName, skiprows = 3, unpack=True)
    p, = ax.plot(x, u, lw=5)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-0.1, 1.1)
    #ax.set_title("CFL = %5.3f, Times = %5.3f"%cfl, %time)
    ax.set_title("")
    ax.set_title("CFL = %5.2f"%cfl + ", Times = %5.3f"%time);
    fig.savefig(fileName.replace(".dat", ".png"))

os.system("eog " + fileNames[0].replace(".dat",".png"))
