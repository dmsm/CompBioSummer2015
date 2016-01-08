# dp.py
# Ran Libeskind-Hadas, June 2015
# The basic DP algorithm for reconciling pairs of trees

# Altered and expanded by Carter Slocum and Annalise Schweickart


# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".
# Edited by Annalise Schweickart and Carter Slocum, July 2015 to return
# the DTL reconciliation graph that uses frequency scoring, as well as the
# number of reconciliations of the host and parasite trees

import newickFormatReader
import greedy
import copy

Infinity = float('inf')


def preorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a
    list of the edges in that tree in preorder (high edges to low edges)"""

    value = tree[rootEdgeName]
    _, _, leftChildEdgeName, rightChildEdgeName = value

    # base case
    if leftChildEdgeName is None:  # then rightChildEdgeName is None also
        return [rootEdgeName]
    else:
        return [rootEdgeName] + \
               preorder(tree, leftChildEdgeName) + \
               preorder(tree, rightChildEdgeName)


def postorder(tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns a
    list of the edges in that tree in postorder (low edges to high edges)"""

    value = tree[rootEdgeName]
    _, _, leftChildEdgeName, rightChildEdgeName = value
    # base case
    if leftChildEdgeName is None:  # then rightChildEdgeName is None also
        return [rootEdgeName]
    else:
        return postorder(tree, leftChildEdgeName) + \
               postorder(tree, rightChildEdgeName) + \
               [rootEdgeName]


def DP(hostTree, parasiteTree, phi, D, T, L):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, and
        duplication cost (D), transfer cost (T), and loss cost (L) and
        returns the DTL graph in the form of a dictionary, as well as a
        the number of maximum parsimony reconciliations. The notation and
        dynamic programming algorithm are explained in the tech report.
        Cospeciation is assumed to cost 0. """

    A = {}  # A, C, O, and bestSwitch are all defined in tech report
    C = {}
    O = {}
    eventsDict = {}  # Dictionary to keep track of events, children, and scores
    bestSwitch = {}
    Minimums = {}  # Dictionary to keep track of minimum reconciliation costs
    oBest = {}  # Dictionary to keep track of the lowest costing events in O
    bestSwitchLocations = {}  # Dictionary to keep track of switch locations
    Score = {}  # Dictionary to calculate the frequency scoring of each event

    # Following logic taken from tech report
    for ep in postorder(parasiteTree, "pTop"):
        for eh in postorder(hostTree, "hTop"):
            _, vp, ep1, ep2 = parasiteTree[ep]
            _, vh, eh1, eh2 = hostTree[eh]
            eventsDict[(vp, vh)] = []
            oBest[(vp, vh)] = []
            # is vp a tip?
            if ep1 is None:  # then ep2 is None too and vp is a tip!
                vpIsATip = True
                pChild1 = None
                pChild2 = None
            else:
                vpIsATip = False
                pChild1 = parasiteTree[ep][2][1]
                pChild2 = parasiteTree[ep][3][1]

            # is vh a tip?
            if eh1 is None:  # then eh2 is None too and vh is a tip!
                vhIsATip = True
                hChild1 = None
                hChild2 = None
            else:
                vhIsATip = False
                hChild1 = hostTree[eh][2][1]
                hChild2 = hostTree[eh][3][1]

            # Compute A(ep, eh)

            if vhIsATip:
                if vpIsATip and phi[vp] == vh:
                    A[(ep, eh)] = 0
                    # Contemporary event to be added to eventsDict
                    Amin = [["C", (None, None), (None, None), 1.0]]
                    Score[(vp, vh)] = 1.0
                else:
                    Score[(vp, vh)] = Infinity
                    A[(ep, eh)] = Infinity
                    Amin = [Infinity]
            else:  # vh is not a tip
                # Compute S and create event list to add to eventsDict
                if not vpIsATip:
                    COepeh = min(C[(ep1, eh1)] + C[(ep2, eh2)],
                                 C[(ep1, eh2)] + C[(ep2, eh1)])
                    coMin = []  # List to keep track lowest cost speciation
                    if COepeh == C[(ep2, eh1)] + C[(ep1, eh2)]:
                        coMin.append(["S", (pChild2, hChild1),
                                      (pChild1, hChild2), (Score[(pChild2, hChild1)] *
                                                           Score[(pChild1, hChild2)])])
                    if COepeh == C[(ep1, eh1)] + C[(ep2, eh2)]:
                        coMin.append(["S", (pChild1, hChild1),
                                      (pChild2, hChild2), (Score[(pChild1, hChild1)]
                                                           * Score[(pChild2, hChild2)])])

                else:
                    COepeh = Infinity
                    coMin = [Infinity]
                    Score[(vp, vh)] = Infinity
                # Compute L and create event list to add to eventsDict
                LOSSepeh = L + min(C[(ep, eh1)], C[(ep, eh2)])
                lossMin = []  # List to keep track of lowest cost loss
                if LOSSepeh == L + C[(ep, eh1)]: lossMin.append(
                        ["L", (vp, hChild1), (None, None),
                         Score[(vp, hChild1)]])
                if LOSSepeh == L + C[(ep, eh2)]: lossMin.append(
                        ["L", (vp, hChild2), (None, None), Score[(vp, hChild2)]])

                # Determine which event occurs for A[(ep, eh)]
                A[(ep, eh)] = min(COepeh, LOSSepeh)
                # Record event occuring for A[(ep, eh)] as Amin
                if COepeh < LOSSepeh:
                    Amin = coMin
                elif LOSSepeh < COepeh:
                    Amin = lossMin
                else:
                    Amin = lossMin + coMin

            # Compute C(ep, eh)
            # First, compute D
            if not vpIsATip:
                DUPepeh = D + C[(ep1, eh)] + C[(ep2, eh)]
                # List to keep track of lowest cost duplication event
                dupList = ["D", (pChild1, vh), (pChild2, vh),
                           (Score[(pChild1, vh)] * Score[(pChild2, vh)])]
            else:
                DUPepeh = Infinity
                dupList = [Infinity]
            # Next, Compute T and create event list to add
            # to eventsDict using bestSwitchLocations
            if not vpIsATip:
                switchList = []  # List to keep track of lowest cost switch
                SWITCHepeh = T + min(C[(ep1, eh)] + bestSwitch[(ep2, eh)],
                                     C[(ep2, eh)] + bestSwitch[(ep1, eh)])
                # if ep2 switching has the lowest cost
                if (C[(ep1, eh)] + bestSwitch[(ep2, eh)]) < (C[(ep2, eh)] +
                                                                 bestSwitch[(ep1, eh)]):
                    for location in bestSwitchLocations[(pChild2, vh)]:
                        currentLoc = location[1]  # Switch landing site
                        if currentLoc is None:  # Switches to a leaf
                            Score[(pChild1, currentLoc)] = Infinity
                            Score[(pChild2, currentLoc)] = Infinity
                        switchList.append(["T", (pChild1, vh), (pChild2,
                                                                currentLoc), (Score[(pChild1, vh)] *
                                                                              Score[(pChild2, currentLoc)])])
                # if ep1 switching has the lowest cost
                elif (C[(ep2, eh)] + bestSwitch[(ep1, eh)]) < (C[(ep1, eh)] +
                                                                   bestSwitch[(ep2, eh)]):
                    for location in bestSwitchLocations[(pChild1, vh)]:
                        currentLoc = location[1]
                        if currentLoc is None:
                            Score[(pChild1, currentLoc)] = Infinity
                            Score[(pChild2, currentLoc)] = Infinity
                        switchList.append(["T", (pChild2, vh),
                                           (pChild1, currentLoc), (Score[(pChild2, vh)] *
                                                                   Score[(pChild1, currentLoc)])])
                # if ep1 switching has the same cost as ep2 switching
                else:
                    for location in bestSwitchLocations[(pChild2, vh)]:
                        currentLoc = location[1]
                        if currentLoc is not None:
                            switchList.append(["T", (pChild1, vh),
                                               (pChild2, currentLoc), (Score[(pChild1, vh)] *
                                                                       Score[(pChild2, currentLoc)])])
                        else:
                            switchList.append(["T", (pChild1, vh),
                                               (pChild2, currentLoc), Infinity])
                    for location in bestSwitchLocations[(pChild1, vh)]:
                        currentLoc = location[1]
                        if currentLoc is not None:
                            switchList.append(["T", (pChild2, vh),
                                               (pChild1, currentLoc), (Score[(pChild2, vh)] *
                                                                       Score[(pChild1, currentLoc)])])
                        else:
                            switchList.append(["T", (pChild1, vh),
                                               (pChild2, currentLoc), Infinity])

            else:
                SWITCHepeh = Infinity
                switchList = [Infinity]
            # Compute C[(ep, eh)] and add the event or events with that cost
            # to the dictionary eventsDict
            C[(ep, eh)] = min(A[(ep, eh)], DUPepeh, SWITCHepeh)
            Minimums[(vp, vh)] = C[(ep, eh)]
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == DUPepeh:
                eventsDict[(vp, vh)].append(dupList)
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == SWITCHepeh:
                eventsDict[(vp, vh)].extend(switchList)
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == A[(ep, eh)]:
                eventsDict[(vp, vh)].extend(Amin)
            for key in eventsDict:
                mapScore = 0  # initialize frequency scoring for each event
                for event in eventsDict[key]:
                    if type(event) is list:
                        mapScore += event[-1]
                Score[key] = mapScore
            if Minimums[(vp, vh)] == Infinity:
                del Minimums[(vp, vh)]
                del eventsDict[(vp, vh)]
            # Compute O(ep, eh)
            # Compute oBest[(vp, vh)], the source of O(ep, eh)
            if vhIsATip:
                O[(ep, eh)] = C[(ep, eh)]
                oBest[(vp, vh)] = [(vp, vh)]
            else:
                # finds the minimum switch locations for O
                oMin = [C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)]].index \
                    (min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)]))
                if oMin == 0:
                    oBest[(vp, vh)].append((vp, vh))
                if oMin == 1:
                    oBest[(vp, vh)].extend(oBest[(vp, hChild1)])
                if oMin == 2:
                    oBest[(vp, vh)].extend(oBest[(vp, hChild2)])

                    # finds Minimum Cost for O
                O[(ep, eh)] = min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)])

        # Compute bestSwitch values
        bestSwitch[(ep, "hTop")] = Infinity
        bestSwitchLocations[(vp, hostTree["hTop"][1])] = [(None, None)]
        for eh in preorder(hostTree, "hTop"):
            _, vp, ep1, ep2 = parasiteTree[ep]
            _, vh, eh1, eh2 = hostTree[eh]

            # is vp a tip?
            if ep1 == None:
                vpIsATip = True
                pChild1 = None
                pChild2 = None
            else:
                vpIsATip = False
                pChild1 = parasiteTree[ep][2][1]
                pChild2 = parasiteTree[ep][3][1]

            # is vh a tip?
            if eh1 is None:  # then eh2 is None too and vh is a tip!
                vhIsATip = True
                hChild1 = None
                hChild2 = None
            else:
                vhIsATip = False
                hChild1 = hostTree[eh][2][1]
                hChild2 = hostTree[eh][3][1]
            # find best place for a switch to occur (bestSwitch)
            # and the location to which the edge switches (bestSwitchLocations)
            if eh1 is not None and eh2 is not None:  # not a tip
                bestSwitchLocations[(vp, hChild1)] = []
                bestSwitchLocations[(vp, hChild2)] = []
                bestSwitch[(ep, eh1)] = min(bestSwitch[(ep, eh)], O[(ep, eh2)])
                bestSwitch[(ep, eh2)] = min(bestSwitch[(ep, eh)], O[(ep, eh1)])
                if bestSwitch[(ep, eh1)] == bestSwitch[(ep, eh)] and \
                                bestSwitchLocations[(vp, vh)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild1)].extend \
                        (bestSwitchLocations[(vp, vh)])
                if bestSwitch[(ep, eh1)] == O[(ep, eh2)] and \
                                oBest[(vp, hChild2)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild1)].extend \
                        (oBest[(vp, hChild2)])
                if bestSwitch[(ep, eh2)] == bestSwitch[(ep, eh)] and \
                                bestSwitchLocations[(vp, vh)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild2)].extend \
                        (bestSwitchLocations[(vp, vh)])
                if bestSwitch[(ep, eh2)] == O[(ep, eh1)] and \
                                oBest[(vp, hChild1)] != [(None, None)]:
                    bestSwitchLocations[(vp, hChild2)].extend \
                        (oBest[(vp, hChild1)])
    for key in bestSwitchLocations:
        if bestSwitchLocations[key][0] == (None, None):
            bestSwitchLocations[key] = bestSwitchLocations[key][1:]
    # Add the costs of each event to the corresponding eventsDict entry
    for key in eventsDict:
        eventsDict[key].append(Minimums[key])

    # Use findPath and findBestRoots to construct the DTL graph dictionary
    treeMin = findBestRoots(parasiteTree, Minimums)
    DTL = findPath(treeMin, eventsDict, {})
    for key in Score.keys():
        if key not in DTL:
            del Score[key]
    DTL, numRecon = addScores(treeMin, DTL, Score)
    return DTL, numRecon


def preorderDTLsort(DTL, ParasiteRoot):
    """This takes in a DTL dictionary and parasite root and returns a sorted list, orderedKeysL, that is ordered
    by level from largest to smallest, where level 0 is the root and the highest level has tips."""

    keysL = greedy.orderDTL(DTL, ParasiteRoot)
    orderedKeysL = []
    levelCounter = 0
    while len(orderedKeysL) < len(keysL):
        for mapping in keysL:
            if mapping[-1] == levelCounter:
                orderedKeysL = orderedKeysL + [mapping]
        levelCounter += 1

    return orderedKeysL


def preorderCheck(preOrderList):
    """This takes a list from preorderDTLsort and removes the duplicate tuples"""
    newList = []
    preDict = {}
    for root in preOrderList:
        if not root in newList:
            newList.append(root)
    for x in range(len(newList)):
        currentRoot = newList[x][0]
        currentLevel = newList[x][1]
        if currentRoot in preDict:
            if preDict[currentRoot][0] > currentLevel:
                newList[x] = (None, None)
            else:
                location = preDict[currentRoot][1]
                newList[location] = (None, None)
        else:
            preDict[currentRoot] = (currentLevel, x)
    return newList


def addScores(treeMin, DTLDict, ScoreDict):
    """Takes the list of reconciliation roots, the DTL reconciliation graph,
    a dictionary of parent nodes, and a dictionary of score values, and
    returns the DTL with the normalized frequency scores calculated."""
    newDTL = copy.deepcopy(DTLDict)
    parentsDict = {}
    preOrder1 = preorderDTLsort(DTLDict, treeMin[0][0])
    preOrder2 = preorderCheck(preOrder1)
    for root in preOrder2:
        if root != (None, None):
            vertices = root[0]
            if root[1] == 0:
                parentsDict[vertices] = ScoreDict[vertices]
            for n in range(len(DTLDict[vertices]) - 1):
                _, child1, child2, oldScore = DTLDict[vertices][n]
                newDTL[vertices][n][3] = parentsDict[vertices] * \
                                         (oldScore / ScoreDict[vertices])
                if child1 != (None, None):
                    if child1 in parentsDict:
                        parentsDict[DTLDict[vertices][n][1]] += \
                            newDTL[vertices][n][3]
                    else:
                        parentsDict[child1] = newDTL[vertices][n][3]
                if child2 != (None, None):
                    if child2 in parentsDict:
                        parentsDict[child2] += newDTL[vertices][n][3]
                    else:
                        parentsDict[child2] = newDTL[vertices][n][3]
    normalize = newDTL[preOrder2[-1][0]][0][-1]
    for key in newDTL:
        for event in newDTL[key][:-1]:
            event[-1] = event[-1] / normalize
    return newDTL, normalize


def findBestRoots(Parasite, MinimumDict):
    """Takes Parasite Tree and a dictionary of minimum reconciliation costs
    and returns a list of the minimum cost reconciliation tree roots"""
    treeTops = [key for key in MinimumDict if key[0] == Parasite['pTop'][1]]
    minimum = min([MinimumDict[root] for root in treeTops])
    treeMin = [pair for pair in treeTops if MinimumDict[pair] == minimum]
    return treeMin


def findPath(tupleList, eventDict, uniqueDict):
    """Takes as input tupleList, a list of minimum reconciliation cost roots,
     eventDict, the dictionary of events and children for each node, and
     uniqueDict, the dictionary of unique vertex mappings. This returns the
     completed DTL graph as a Dictionary"""
    for vertexPair in tupleList:
        if not vertexPair in uniqueDict:
            uniqueDict[vertexPair] = eventDict[vertexPair]
        for event in eventDict[vertexPair][:-1]:
            for location in event:
                if type(location) is tuple and location != (None, None):
                    findPath([location], eventDict, uniqueDict)
    return uniqueDict


def reconcile(fileName, D, T, L):
    """Takes as input a newick file, FileName, a dupliction cost, a transfer
    cost, and a loss cost. This uses newickFormatReader to extract the host
    tree, parasite tree and tip mapping from the file and then calls DP to
    return the DTL reconciliation graph of the provided newick file"""
    host, paras, phi = newickFormatReader.getInput(fileName)
    return DP(host, paras, phi, D, T, L)

