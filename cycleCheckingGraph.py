# cycleCheckingGraph.py

# Srinidhi Srinivasan, Juliet Forman
# July 2015

# This file contains function for building the cycle checking graph, 
# cycleCheckingGraph, which is in the form of a dictionary. This dictionary 
# has keys that are nodes and values that are a list of all the children. The 
# cycleCheckingGraph represents edges between nodes that show the temporal 
# relationship between the host tree and the parasite Tree. The main function
# in this file is buildReconstruction and the rest of the functions are 
# helper function that are used by buildReconstruction

# This cycle checking graph is known in the paper as a temporal feasability 
# graph

from Greedy import findRoot


def InitDicts(tree):
    """This function takes as input a tree dictionary and returns a dictionary
    with all of the bottom nodes of the edges as keys and empty lists as
    values."""

    treeDict = {}
    for key in tree:
        if key == 'pTop' or key == 'hTop':
            treeDict[tree[key][1]] = []
        else:
            treeDict[key[1]] = []
    return treeDict


def treeFormat(tree):
    """Takes a tree in the format that it comes out of newickFormatReader and
    converts it into a dictionary with keys which are the bottom nodes of the
    edges and values which are the children. The values are in the form of a
    list of names of nodes, each of which is in a string."""

    treeDict = InitDicts(tree)
    treeRoot = findRoot(tree)
    for key in tree:
        # deal with case where the key is not in tuple form
        if key == 'hTop' or key == 'pTop':  #check if key is a root
            #check if the first child of key is None
            if tree[key][-2] == None:
                #add None to the list of children of key in treeDict
                treeDict[treeRoot] = treeDict[treeRoot] + [tree[key][-2]]
            else:
                #add the host node of the first child to the children of key
                treeDict[treeRoot] = treeDict[treeRoot] + [tree[key][-2][1]]
            #check if the second child of key is None
            if tree[key][-1] == None:
                #add None to the list of children of key in treeDict
                treeDict[treeRoot] = treeDict[treeRoot] + [tree[key][-1]]
            else:
                #add the host node of the second child to the children of key
                treeDict[treeRoot] = treeDict[treeRoot] + [tree[key][-1][1]]

        else:  #where key is in tuple form
            #check if the first child of key is None
            if tree[key][-2] == None:
                #add None to the list of children of key in treeDict
                treeDict[key[1]] = treeDict[key[1]] + [tree[key][-2]]
            else:
                #add the host node of the first child to the children of key
                treeDict[key[1]] = treeDict[key[1]] + [tree[key][-2][1]]
            #check if the second child of key is None
            if tree[key][-1] == None:
                #add None to the list of children of key in treeDict
                treeDict[key[1]] = treeDict[key[1]] + [tree[key][-1]]
            else:
                #add the host node of the second child to the children of key
                treeDict[key[1]] = treeDict[key[1]] + [tree[key][-1][1]]
    return treeDict


def createParentsDict(H, P):
    """Takes a host and a parasite tree with edges as keys and returns a
    dictionary with keys which are the bottom nodes of those edges and values
    which are the top nodes of those edges."""

    parentsDict = {}
    for key in H:
        if key == 'hTop':
            parentsDict[H[key][1]] = H[key][0]
        else:
            parentsDict[key[1]] = H[key][0]
    for key in P:
        if key == 'pTop':
            parentsDict[P[key][1]] = P[key][0]
        else:
            parentsDict[key[1]] = P[key][0]
    return parentsDict


def uniquify(list):
    """Takes as input a list and returns a list containing only the unique
    elements of the input list."""

    holdDict = {}
    for thing in list:
        holdDict[thing] = 1
    return holdDict.keys()


def buildReconciliation(HostTree, ParasiteTree, reconciliation):
    """Takes as input a host tree, a parasite tree, and a reconciliation, and
    returns a graph where the keys are host or parasite nodes, and the values
    are a list of the children of a particular node. The graph represents
    temporal relationships between events."""

    # create a dictionary with a list of parents of each host and parasite node
    parents = createParentsDict(HostTree, ParasiteTree)
    H = treeFormat(HostTree)
    P = treeFormat(ParasiteTree)
    cycleCheckingGraph = H
    cycleCheckingGraph.update(P)
    for key in reconciliation:
        #deal with transfer case:
        if reconciliation[key][0] == 'T':
            #add the children of the parasite node to the list of children
            #of the host node in cycleCheckingGraph
            cycleCheckingGraph[key[0]] = P[key[0]] + \
                                         [reconciliation[key][1][1], reconciliation[key][2][1]]
            #find the parents of the take-off and landing host nodes
            parent1 = parents[reconciliation[key][1][1]]
            parent2 = parents[reconciliation[key][2][1]]
            #add the parasite node as a child of parent1 and parent2
            cycleCheckingGraph[parent1] = cycleCheckingGraph[parent1] + \
                                          [key[0]]
            cycleCheckingGraph[parent2] = cycleCheckingGraph[parent2] + \
                                          [key[0]]

        #deal with speciation case:
        elif reconciliation[key][0] == 'S':
            parent = parents[key[0]]
            if parent != 'Top':
                cycleCheckingGraph[parent] = cycleCheckingGraph[parent] + \
                                             [key[1]]
            cycleCheckingGraph[key[1]] = cycleCheckingGraph[key[1]] + \
                                         cycleCheckingGraph[key[0]]
        #deal with duplication case:
        elif reconciliation[key][0] == 'D':
            parent = parents[key[1]]
            if parent != 'Top':
                cycleCheckingGraph[parent] = cycleCheckingGraph[parent] + \
                                             [key[0]]
            cycleCheckingGraph[key[0]] = cycleCheckingGraph[key[0]] + [key[1]]
        #deal with contemporary case:
        elif reconciliation[key][0] == 'C':
            cycleCheckingGraph[key[1]] = [None]
            cycleCheckingGraph[key[0]] = [None]

    for key in cycleCheckingGraph:
        cycleCheckingGraph[key] = uniquify(cycleCheckingGraph[key])

    return cycleCheckingGraph

