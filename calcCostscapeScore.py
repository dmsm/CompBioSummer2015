# calcCostscapeScore.py
# Juliet Forman and Srinidhi Srinivasan
# July 2015

# This file contains functions for calculating scores for events in the 
# DTLReconGraph based on their frequency among the different regions of 
# costscape. The score for each event is the fraction of regions in costscape 
# in which the event appears. The main function is newScoreWrapper, and the 
# other functions are helper functions for newScoreWrapper.

from dp import DP
from costscapeScore import findCenters
import newickFormatReader


def getDTLReconGraphVals(pointList):
    """This function takes as input a list of centroids of the costscape
    regions, and returns a list of tuples containing the T, L costs for each
    region."""
    pointList = [point.replace(',', '') for point in pointList]
    DTLReconGraphPairs = []
    for point in pointList:
        coordList = point[7:-1].split()
        pair = [float(coord) for coord in coordList]
        DTLReconGraphPairs.append(tuple(pair))
    return DTLReconGraphPairs


def getCostscapeDTLReconGraphs(DTLReconGraphPairs, hostTree, parasiteTree,
                               phi):
    """This function takes as input DTLReconGraphPairs, a list of tuples with
    transfer and loss costs, and the hostTree, parasiteTree, and phi. It
    returns a list of DTLReconGraphs whose scores are computed with the
    transfer and loss values from each element in DTLReconGraphPairs."""

    DTLReconGraphList = []
    for cost in DTLReconGraphPairs:
        # assign those associated costs to the newDTLReconGraph
        newDTLReconGraph = DP(hostTree, parasiteTree, phi, 1, cost[0],
                              cost[1])[0]
        DTLReconGraphList.append(newDTLReconGraph)
    return DTLReconGraphList


def changeDTLReconGraphScores(originalDTLReconGraph, DTLReconGraphList):
    """This function takes as input the originalDTLReconGraph and a list
    DTLReconGraphList of the DTLReconGraphs from each region in costscape.
    This function calculates a new score for each event in
    originalDTLReconGraph, and returns a newDTLReconGraph with these
    scores."""

    newDTLReconGraph = {}
    numDTLReconGraph = len(DTLReconGraphList)
    for event in originalDTLReconGraph:
        counter = 0
        for DTLReconGraph in DTLReconGraphList:
            if event in DTLReconGraph:
                counter += 1
        newScore = 1.0 * counter / numDTLReconGraph
        oldVal = originalDTLReconGraph[event]
        oldVal[0][-1] = newScore  # assign new score to DTLReconGraph
        newDTLReconGraph[event] = oldVal
    return newDTLReconGraph


def newScoreWrapper(newickFile, switchLo, switchHi, lossLo, lossHi, D, T, L):
    """This function takes as input hostTree, parasiteTree, phi, duplication
    cost D, transfer cost T, and loss cost L, and returns the
    newDTLReconGraph whose scores were calculated from costscape."""

    H, P, phi = newickFormatReader.getInput(newickFile)
    originalDTLReconGraph, _ = DP(H, P, phi, D, T, L)
    pointList = findCenters(newickFile, switchLo, switchHi, lossLo, lossHi)
    DTLReconGraphPairs = getDTLReconGraphVals(pointList)
    DTLReconGraphList = getCostscapeDTLReconGraphs(DTLReconGraphPairs, H, P,
                                                   phi)
    newDTLReconGraph = changeDTLReconGraphScores(originalDTLReconGraph,
                                                 DTLReconGraphList)
    return newDTLReconGraph
