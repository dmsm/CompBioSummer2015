# costscapeScore.py
# Ran Libeskind-Hadas, Jessica Yi-Chieh Wu, Mukul Bansal, November 2013
# Updated by Juliet Forman and Srinidhi Srinivasan, July 2015

# This file contains functions that find the centers of the regions in 
# costscape for a given .newick file. The main function is findCenters, and 
# the rest are helper functions.

# python libraries
import time
import math


from newickFormatReader import *

# xscape libraries
try:
    import xscape
except ImportError:
    import sys
    from os.path import realpath, dirname, join
    sys.path.append(join(realpath(dirname(dirname(__file__))), "python"))
    import xscape
from xscape.commonAnalytic import *
from xscape.CostVector import *
from xscape import reconcile
from xscape import plotcostsAnalyticNew as plotcosts


import shapely
from shapely.wkt import loads as load_wkt


def getNewCoordList(newickFile, switchLo, switchHi, lossLo, lossHi):
    """Takes as input a newick file in the form <filename>.newick, and low 
    and high values for costscape for both switches and losses. Returns a 
    list of strings, where each string contains all the verteces of one 
    region from costscape."""

    hostTree, parasiteTree, phi = newickFormatReader(newickFile)
    CVlist = reconcile.reconcile(parasiteTree, hostTree, phi, switchLo, \
        switchHi, lossLo, lossHi)
    coordList = plotcosts.plotcosts(CVlist, lossLo, lossHi, switchLo, \
        switchHi, "", False, False)
    print coordList
    newCoordList = []
    for vertexList in coordList:
        string = "POLYGON(("
        for vertex in vertexList:
            string = string + str(vertex[0]) + ' ' + str(vertex[1]) + ','
        string = string[:-1] + '))'
        newCoordList.append(string)
    return newCoordList


def findCenters(newickFile, switchLo, switchHi, lossLo, lossHi):
    """This function takes as input a .newick file in the form 
    <filename>.newick, and low and high values for costscape for both 
    switches and losses. It returns a list of the centroids of each region in 
    the costscape associated with the given .newick file."""

    hostTree, parasiteTree, phi = newickFormatReader(newickFile)
    CVlist = reconcile.reconcile(parasiteTree, hostTree, phi, switchLo, \
        switchHi, lossLo, lossHi)
    coordList = plotcosts.plotcosts(CVlist, lossLo, lossHi, switchLo, \
        switchHi, "", False, False)
    polygonList = getNewCoordList(newickFile, switchLo, switchHi, lossLo, \
        lossHi)
    pointList = []
    for i in range(len(polygonList)):
        point = polygonList[i]
        numCommas = 0
        for j in range(len(point)):
            if point[j] == ",":
                numCommas = numCommas + 1
        if numCommas > 1:
            #polygon case
            region = load_wkt(point)
            pointList.append((region.centroid.wkt))
        elif numCommas == 1:
            #line case
            x1 = coordList[i][0][0]
            y1 = coordList[i][0][1]
            x2 = coordList[i][1][0]
            y2 = coordList[i][1][1]
            midx = (x1 + x2)*1.0/2
            midy = (y1 + y2)*1.0/2
            pointList.append("POINT (" + str(midx) + " " + str(midy) + ")")
        else:
            #point case
            pointList.append("POINT " + str(coordList[i][0]))

    return pointList
