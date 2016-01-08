# detectCycles.py
# Juliet Forman and Srinidhi Srinivasan
# July 2015

# This file contains functions for detecting and marking temporal 
# inconsistencies in reconciliattions. It uses depth first search to find 
# cycles in graphs that represent temporal relationships between nodes in the 
# parasite and host trees. When temporal inconsistencies are found, the 
# reconciliation is updated to mark which transfers are guilty of causing
# temporal inconsistencies. The main function in this file is 
# detectCyclesWrapper, and the rest of the functions are helper functions 
# that are used by detectCyclesWrapper.

import copy
from cycleCheckingGraph import createParentsDict, treeFormat, findRoot


def buildReconciliation(HostTree, ParasiteTree, reconciliation):
    """Takes as input a host tree, a parasite tree, and a reconciliation, and
    returns a graph where the keys are host or parasite nodes, and the values
    are a list of the children of a particular node. The graph represents
    temporal relationships between events. The function also returns a list
    transferList containing all the transfers in the reconciliation in the
    form """

    # create a dictionary with a list of parents of each host and parasite node
    parents = createParentsDict(HostTree, ParasiteTree)
    H = treeFormat(HostTree)
    P = treeFormat(ParasiteTree)
    cycleCheckingGraph = H
    cycleCheckingGraph.update(P)
    transferList = []
    for key in reconciliation:
        # deal with transfer case:
        if reconciliation[key][0] == 'T':
            # add the children of the parasite node to the list of children
            # of the host node in cycleCheckingGraph
            cycleCheckingGraph[key[0]] = P[key[0]] + \
                                         [reconciliation[key][1][1], reconciliation[key][2][1]]
            # find the parents of the take-off and landing host nodes
            parent1 = parents[reconciliation[key][1][1]]
            parent2 = parents[reconciliation[key][2][1]]
            # add the parasite node as a child of parent1 and parent2
            cycleCheckingGraph[parent1] = cycleCheckingGraph[parent1] + \
                                          [key[0]]
            cycleCheckingGraph[parent2] = cycleCheckingGraph[parent2] + \
                                          [key[0]]
            transferEdge1 = reconciliation[key][1][1]
            transferEdge2 = reconciliation[key][2][1]
            transferList.append([key[0], parent1, transferEdge1, parent2, \
                                 transferEdge2])

        # deal with speciation case:
        elif reconciliation[key][0] == 'S':
            parent = parents[key[0]]
            if parent != 'Top':
                cycleCheckingGraph[parent] = cycleCheckingGraph[parent] + \
                                             [key[1]]
            cycleCheckingGraph[key[1]] = cycleCheckingGraph[key[1]] + \
                                         cycleCheckingGraph[key[0]]

        # deal with duplication case:
        elif reconciliation[key][0] == 'D':
            parent = parents[key[1]]
            if parent != 'Top':
                cycleCheckingGraph[parent] = cycleCheckingGraph[parent] + \
                                             [key[0]]
            cycleCheckingGraph[key[0]] = cycleCheckingGraph[key[0]] + [key[1]]

        # deal with contemporary case:
        elif reconciliation[key][0] == 'C':
            cycleCheckingGraph[key[1]] = [None]
            cycleCheckingGraph[key[0]] = [None]

    for key in cycleCheckingGraph:
        cycleCheckingGraph[key] = list(set(cycleCheckingGraph[key]))

    return cycleCheckingGraph, transferList


def detectCycles(HostTree, ParasiteTree, reconciliation):
    """This function takes as input the cycle checking graph,
    cycleCheckingGraph. It returns a new version of cycleCheckingGraph,
    newCycleCheckingGraph, from which the transfer events responsible for the
    cycles have been removed. It also returns a list, guiltyTransferList, of
    the guilty transfers."""

    guiltyTransferList = []
    markingDict = {}
    cycleCheckingGraph, transferList = buildReconciliation(HostTree,
                                                           ParasiteTree, reconciliation)
    Hroot = findRoot(HostTree)
    markingDict[Hroot] = ['check']
    cycleEdge = recurseChildren(cycleCheckingGraph, markingDict, Hroot)
    newCycleCheckingGraph, guiltyTransfer, transferList = deleteTransfer(
            cycleCheckingGraph, markingDict, transferList, cycleEdge)
    if guiltyTransfer:
        guiltyTransferList.append(guiltyTransfer)
    while cycleEdge is not None:
        markingDict = {}
        cycleEdge = recurseChildren(newCycleCheckingGraph,
                                    {Hroot: ['check']}, Hroot)
        if cycleEdge is None:
            for node in newCycleCheckingGraph:
                if not checked(markingDict, node):
                    check(markingDict, node)
                    cycleEdge = recurseChildren(newCycleCheckingGraph,
                                                markingDict, node)
                    if cycleEdge is not None:
                        break
        newCycleCheckingGraph, guiltyTransfer, transferList = deleteTransfer(
                newCycleCheckingGraph, markingDict, transferList, cycleEdge)
        if guiltyTransfer:
            guiltyTransferList.append(guiltyTransfer)

    return newCycleCheckingGraph, guiltyTransferList


def checked(markingDict, node):
    """This function takes as input a markingDict and a node, and checks the
    node in markingDict, marking it as visited."""

    return node in markingDict


def ticked(markingDict, node):
    """This function takes as input a markingDict and a checked node, and
    returns True if the node is already ticked, False if it is not."""

    return node in markingDict and len(markingDict[node]) == 2


def tick(markingDict, node):
    """This function takes as input a markingDict and node which is checked
    but not ticked, and ticks the node in markingDict."""

    markingDict[node] = markingDict[node] + ['tick']


def untick(markingDict, node):
    """This function takes as input a markingDict and a ticked node, and
    unticks the node in markingDict."""

    markingDict[node] = markingDict[node][:1]


def check(markingDict, node):
    """This function takes as input markingDict and a node, and checks the
    node in markingDict."""

    markingDict[node] = ['check']


def recurseChildren(cycleCheckingGraph, markingDict, node):
    """This function takes as input the cycle checking graph
    cycleCheckingGraph, markingDict, a dictionary that keeps track of all the
    childNodes that are marked or ticked, and node, the node that we will
    recurse on. Nodes are marked permanently as soon as a recursive call is
    made on them. Ticks are added to a node when a recursive call is made on
    them, and then removed as soon as the recursive call is finished. A cycle
    is detected whenever a child of the current node is already ticked. The
    function updates markingDict and returns a tuple where the first element
    is the parent of the ticked child and the second element is that child.
    If the function finds no cycles, it returns None."""

    tick(markingDict, node)
    for child in cycleCheckingGraph[node]:
        if not checked(markingDict, child) and child is not None:
            check(markingDict, child)
            cycleEdge = recurseChildren(cycleCheckingGraph, markingDict,
                                        child)

            if cycleEdge is not None:
                return cycleEdge

        elif child is not None:
            if ticked(markingDict, child):
                return node, child

    untick(markingDict, node)

    return None


def deleteTransfer(cycleCheckingGraph, markingDict, transferList, cycleEdge):
    """This function takes as input the cycle checking graph
    cycleCheckingGraph, a dictionary markingDict, a list transferList of all
    transfers in the reconciliation, and cycleEdge, which is either a tuple
    with two elements, or None. The function returns a new cycle checking
    graph newCycleCheckingGraph, from which the guilty transfer has been
    removed. It also returns the guilty transfer and transferList with the
    guilty transfer removed."""

    newCycleCheckingGraph = copy.deepcopy(cycleCheckingGraph)
    guiltyTransfer = []
    if cycleEdge is None:
        return newCycleCheckingGraph, guiltyTransfer, transferList
    node, cycleNode = cycleEdge

    for transfer in transferList:

        if cycleNode in transfer:

            guiltyTransfer = transfer
            transferList.remove(transfer)
            # remove all the edges that were added due to the guilty transfer
            removeChild(newCycleCheckingGraph, transfer[1], transfer[0])
            removeChild(newCycleCheckingGraph, transfer[0], transfer[2])
            removeChild(newCycleCheckingGraph, transfer[0], transfer[4])
            if transfer[1] != transfer[3]:
                removeChild(newCycleCheckingGraph, transfer[3], transfer[0])
            break

    return newCycleCheckingGraph, guiltyTransfer, transferList


def removeChild(cycleCheckingGraph, parent, child):
    """This function takes as input a graph cycleCheckingGraph, a parent, and
    its child. It removes the edge between the parent and child in
    cycleCheckingGraph."""

    childList = cycleCheckingGraph[parent]
    childList.remove(child)
    cycleCheckingGraph[parent] = childList


def updateReconciliation(guiltyTransferList, HostTree, ParasiteTree,
                         reconciliation):
    """This function takes as input a list guiltyTransferList of transfers
    that are responsible for cycles, a host tree and parasite tree, and the
    original reconciliation of host and parasite. It returns a new
    reconciliation in which all of the guilty transfers have been marked as a
    new event 'GT'."""

    parents = createParentsDict(HostTree, ParasiteTree)
    newReconciliation = copy.deepcopy(reconciliation)
    for transfer in guiltyTransferList:
        for key in newReconciliation:
            if reconciliation[key][0] == 'T':
                # check if transfer and key are the same event
                if transfer[0] == key[0] and (transfer[2] == key[1] or
                                                      transfer[4] == key[1]):
                    # create a new event 'GT' instead of 'T'
                    newValue = ['GT'] + newReconciliation[key][1:]
                    # replace the old event with newValue in the reconciliation
                    newReconciliation[key] = newValue
    return newReconciliation


def detectCyclesWrapper(HostTree, ParasiteTree, reconciliation):
    """This function takes in a host tree, parasite tree, and reconciliation.
    It returns an updated cycle checking graph where the edges that came from
    the guilty transfers have been removed, and it returns a new
    reconciliation where the guilty transfers have been marked."""

    newCycleCheckingGraph, guiltyTransferList = detectCycles(HostTree,
                                                             ParasiteTree, reconciliation)
    newReconciliation = updateReconciliation(guiltyTransferList, HostTree,
                                             ParasiteTree, reconciliation)
    return newReconciliation, newCycleCheckingGraph
