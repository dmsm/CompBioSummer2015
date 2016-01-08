# DP.py
# Ran Libeskind-Hadas, June 2015
# The basic DP algorithm for reconciling pairs of trees

# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".

import newickFormatReader
# import DrawDTL
import copy

Infinity = float('inf')


def preorder(Tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns
        a list of the edges in that tree in preorder (high edges to low edges)"""

    value = Tree[rootEdgeName]
    leftChildEdgeName = value[2]
    rightChildEdgeName = value[3]
    # base case
    if leftChildEdgeName == None:  # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:
        return [rootEdgeName] + \
               preorder(Tree, leftChildEdgeName) + \
               preorder(Tree, rightChildEdgeName)


def postorder(Tree, rootEdgeName):
    """ Takes a tree as input (see format description above) and returns
        a list of the edges in that tree in postorder (low edges to high edges)"""

    value = Tree[rootEdgeName]
    leftChildEdgeName = value[2]
    rightChildEdgeName = value[3]
    # base case
    if leftChildEdgeName == None:  # then rightChildEdgeName == None also
        return [rootEdgeName]
    else:
        return postorder(Tree, leftChildEdgeName) + \
               postorder(Tree, rightChildEdgeName) + \
               [rootEdgeName]


def DP(hostTree, parasiteTree, phi, D, T, L):
    """ Takes a hostTree, parasiteTree, tip mapping function phi, and
        duplication cost (D), transfer cost (T), and loss cost (L) and
        returns the DTL graph in the form of a dictionary, as well as a
        drawing of the DTL graph. Cospeciation is assumed to cost 0. """

    A = {}
    C = {}
    O = {}
    Dictionary = {}
    BestSwitch = {}
    Minimums = {}
    Obest = {}
    BestSwitchLocations = {}
    Score = {}
    Parents = {}

    for ep in postorder(parasiteTree, "pTop"):
        for eh in postorder(hostTree, "hTop"):
            vp = parasiteTree[ep][1]
            vh = hostTree[eh][1]
            ep1 = parasiteTree[ep][2]
            ep2 = parasiteTree[ep][3]
            eh1 = hostTree[eh][2]
            eh2 = hostTree[eh][3]
            Dictionary[(vp, vh)] = []
            Obest[(vp, vh)] = []
            # is vp a tip?
            if ep1 == None:  # then ep2 == None too and vp is a tip!
                vpIsATip = True
                pChild1 = None
                pChild2 = None
            else:
                vpIsATip = False
                pChild1 = parasiteTree[ep][2][1]
                pChild2 = parasiteTree[ep][3][1]

            # is vh a tip?
            if eh1 == None:  # then eh2 == None too and vh is a tip!
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
                    Amin = [["C", (None, None), (None, None), 1]]  # Contemporary event to be added to Dictionary
                    Score[(vp, vh)] = 1
                else:
                    Score[(vp, vh)] = Infinity
                    A[(ep, eh)] = Infinity
                    Amin = [Infinity]
            else:  #vh is not a tip
                # Compute S and create event list to add to Dictionary
                if not vpIsATip:
                    COepeh = min(C[(ep1, eh1)] + C[(ep2, eh2)], \
                                 C[(ep1, eh2)] + C[(ep2, eh1)])
                    coMin = []
                    if COepeh == C[(ep2, eh1)] + C[(ep1, eh2)]:
                        coMin.append(["S", (pChild2, hChild1), (pChild1, hChild2),
                                      (Score[(pChild2, hChild1)] * Score[(pChild1, hChild2)])])
                    if COepeh == C[(ep1, eh1)] + C[(ep2, eh2)]:
                        coMin.append(["S", (pChild1, hChild1), (pChild2, hChild2),
                                      (Score[(pChild1, hChild1)] * Score[(pChild2, hChild2)])])

                else:
                    COepeh = Infinity
                    coMin = [Infinity]
                    Score[(vp, vh)] = Infinity
                # Compute L and create event list to add to Dictionary
                LOSSepeh = L + min(C[(ep, eh1)], C[(ep, eh2)])
                lossMin = []
                if LOSSepeh == L + C[(ep, eh1)]: lossMin.append(
                    ["L", (vp, hChild1), (None, None), Score[(vp, hChild1)]])
                if LOSSepeh == L + C[(ep, eh2)]: lossMin.append(
                    ["L", (vp, hChild2), (None, None), Score[(vp, hChild2)]])

                # Determine which event occurs for A[(ep, eh)]
                A[(ep, eh)] = min(COepeh, LOSSepeh)
                if COepeh < LOSSepeh:
                    Amin = coMin
                elif LOSSepeh < COepeh:
                    Amin = lossMin
                else:
                    Amin = lossMin + coMin

            # Compute C(ep, eh)
            #   First, compute D
            if not vpIsATip:
                DUPepeh = D + C[(ep1, eh)] + C[(ep2, eh)]
                dupList = ["D", (pChild1, vh), (pChild2, vh), (Score[(pChild1, vh)] * Score[(pChild2, vh)])]
            else:
                DUPepeh = Infinity
                dupList = [Infinity]
            #   Next, Compute T and create event list to add to Dictionary using BestSwitchLocations
            if not vpIsATip:
                switchList = []
                SWITCHepeh = T + min(C[(ep1, eh)] + BestSwitch[(ep2, eh)], \
                                     C[(ep2, eh)] + BestSwitch[(ep1, eh)])
                if (C[(ep1, eh)] + BestSwitch[(ep2, eh)]) < (C[(ep2, eh)] + BestSwitch[(ep1, eh)]):
                    for item in BestSwitchLocations[(pChild2, vh)]:
                        if item[1] == None:
                            Score[(pChild1, item[1])] = Infinity
                            Score[(pChild2, item[1])] = Infinity
                        switchList.append(["T", (pChild1, vh), (pChild2, item[1]),
                                           (Score[(pChild1, vh)] * Score[(pChild2, item[1])])])
                elif (C[(ep2, eh)] + BestSwitch[(ep1, eh)]) < (C[(ep1, eh)] + BestSwitch[(ep2, eh)]):
                    for item in BestSwitchLocations[(pChild1, vh)]:
                        if item[1] == None:
                            Score[(pChild1, item[1])] = Infinity
                            Score[(pChild2, item[1])] = Infinity
                        switchList.append(["T", (pChild2, vh), (pChild1, item[1]),
                                           (Score[(pChild2, vh)] * Score[(pChild1, item[1])])])
                else:
                    for item in BestSwitchLocations[(pChild2, vh)]:
                        if item[1] != None:
                            switchList.append(["T", (pChild1, vh), (pChild2, item[1]),
                                               (Score[(pChild1, vh)] * Score[(pChild2, item[1])])])
                        else:
                            switchList.append(["T", (pChild1, vh), (pChild2, item[1]), Infinity])
                    for item in BestSwitchLocations[(pChild1, vh)]:
                        if item[1] != None:
                            switchList.append(["T", (pChild2, vh), (pChild1, item[1]),
                                               (Score[(pChild2, vh)] * Score[(pChild1, item[1])])])
                        else:
                            switchList.append(["T", (pChild1, vh), (pChild2, item[1]), Infinity])

            else:
                SWITCHepeh = Infinity
                switchList = [Infinity]
            # Compute C[(ep, eh)] and add it's source to the Dictionary
            C[(ep, eh)] = min(A[(ep, eh)], DUPepeh, SWITCHepeh)
            Minimums[(vp, vh)] = C[(ep, eh)]
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == DUPepeh:
                Dictionary[(vp, vh)].append(dupList)
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == SWITCHepeh:
                Dictionary[(vp, vh)].extend(switchList)
            if min(A[(ep, eh)], DUPepeh, SWITCHepeh) == A[(ep, eh)]: Dictionary[(vp, vh)].extend(Amin)
            for key in Dictionary.keys():
                mapScore = 0
                for item in Dictionary[key]:
                    if type(item) == list:
                        mapScore += item[-1]
                Score[key] = mapScore
            if Minimums[(vp, vh)] == Infinity:
                del Minimums[(vp, vh)]
                del Dictionary[(vp, vh)]
            # Compute O(ep, eh)
            # Compute Obest[(vp, vh)], the source of O(ep, eh)
            if vhIsATip:
                O[(ep, eh)] = C[(ep, eh)]
                Obest[(vp, vh)] = [(vp, vh)]
            else:
                omin = [C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)]].index(min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)]))
                if omin == 0:
                    Obest[(vp, vh)].append((vp, vh))
                if omin == 1:
                    Obest[(vp, vh)].extend(Obest[(vp, hChild1)])
                if omin == 2:
                    Obest[(vp, vh)].extend(Obest[(vp, hChild2)])
                O[(ep, eh)] = min(C[(ep, eh)], O[(ep, eh1)], O[(ep, eh2)])

        # Compute BestSwitch values
        BestSwitch[(ep, "hTop")] = Infinity
        BestSwitchLocations[(vp, hostTree["hTop"][1])] = [(None, None)]
        for eh in preorder(hostTree, "hTop"):
            vp = parasiteTree[ep][1]
            vh = hostTree[eh][1]
            ep1 = parasiteTree[ep][2]
            ep2 = parasiteTree[ep][3]
            eh1 = hostTree[eh][2]
            eh2 = hostTree[eh][3]
            #is vp a tip?
            if ep1 == None:
                vpIsATip = True
                pChild1 = None
                pChild2 = None
            else:
                vpIsATip = False
                pChild1 = parasiteTree[ep][2][1]
                pChild2 = parasiteTree[ep][3][1]

            # is vh a tip?
            if eh1 == None:  # then eh2 == None too and vh is a tip!
                vhIsATip = True
                hChild1 = None
                hChild2 = None
            else:
                vhIsATip = False
                hChild1 = hostTree[eh][2][1]
                hChild2 = hostTree[eh][3][1]
            # find best place for a switch to occur, and the location to which the edge switches    
            if eh1 != None and eh2 != None:
                BestSwitchLocations[(vp, hChild1)] = []
                BestSwitchLocations[(vp, hChild2)] = []
                BestSwitch[(ep, eh1)] = min(BestSwitch[(ep, eh)], O[(ep, eh2)])
                BestSwitch[(ep, eh2)] = min(BestSwitch[(ep, eh)], O[(ep, eh1)])
                if BestSwitch[(ep, eh1)] == BestSwitch[(ep, eh)] and BestSwitchLocations[(vp, vh)] != [(None, None)]:
                    BestSwitchLocations[(vp, hChild1)].extend(BestSwitchLocations[(vp, vh)])
                if BestSwitch[(ep, eh1)] == O[(ep, eh2)] and Obest[(vp, hChild2)] != [(None, None)]:
                    BestSwitchLocations[(vp, hChild1)].extend(Obest[(vp, hChild2)])
                if BestSwitch[(ep, eh2)] == BestSwitch[(ep, eh)] and BestSwitchLocations[(vp, vh)] != [(None, None)]:
                    BestSwitchLocations[(vp, hChild2)].extend(BestSwitchLocations[(vp, vh)])
                if BestSwitch[(ep, eh2)] == O[(ep, eh1)] and Obest[(vp, hChild1)] != [(None, None)]:
                    BestSwitchLocations[(vp, hChild2)].extend(Obest[(vp, hChild1)])

    for key in BestSwitchLocations.keys():
        if BestSwitchLocations[key][0] == (None, None):
            BestSwitchLocations[key] = BestSwitchLocations[key][1:]
    # Add the costs of each event to the corresponding Dictionary entry
    for key in Dictionary.keys():
        Dictionary[key].append(Minimums[key])

    # Use findPath and findBest to construct the DTL graph dictionary
    treeMin = findBest(parasiteTree, Minimums)
    DTL = {}
    DTL = findPath(treeMin, Dictionary, DTL)
    for key in Score.keys():
        Score[key] = Score[key] * 1.0
        if not key in DTL:
            del Score[key]

    newDTL = copy.deepcopy(DTL)
    DTL, numRecon = addScores(treeMin, DTL, Parents, Score, newDTL)
    # Draw the DTL reconciliation of this DTL Graph

    #DrawDTL.drawNodes(treeMin, DTL, 450, {})

    return DTL


def orderDTL(DTL, ParasiteRoot):
    """This function takes in a DTL graph and the ParasiteRoot. It outputs a list, keysL, that contains tuples. Each tuple 
    has two elements. The first is a mapping node of the form (p, h), where p is a parasite node and h is a host node. The 
    second element is a level representing the depth of that mapping node within the tree."""

    keysL = []
    topNodes = []
    for key in DTL:
        if key[0] == ParasiteRoot:
            topNodes.append(key)
    for vertex in topNodes:
        keysL.extend(orderDTLRoots(DTL, vertex, 0))
    return keysL


def orderDTLRoots(DTL, vertex, level):
    """this function takes a DTL graph, one node, vertex, of the DTL graph, a level, and a dictionary markingDict, and 
    returns a list, keysL, that contains tuples. Each tuple has two elements. The first is a mapping node of the form 
    (p, h), where p is a parasite node and h is a host node. The second element is a level representing the depth of 
    that mapping node within the tree. This function adds the input vertex to keysL and recurses on its children."""

    keysL = []
    for i in range(len(DTL[vertex]) - 1):  #loop through each event associated with key in DTL
        event = DTL[vertex][i]
        child1 = event[1]
        child2 = event[2]
        keysL = keysL + [(vertex, level)]
        if child1[0] != None:
            keysL.extend(orderDTLRoots(DTL, child1, level + 1))
        if child2[0] != None:
            keysL.extend(orderDTLRoots(DTL, child2, level + 1))
    return keysL


def preorderDTLsort(DTL, ParasiteRoot):
    """This takes in a DTL dictionary and parasite root and returns a sorted list, orderedKeysL, that is ordered
    by level from largest to smallest, where level 0 is the root and the highest level has tips."""

    keysL = orderDTL(DTL, ParasiteRoot)
    orderedKeysL = []
    levelCounter = 0
    while len(orderedKeysL) < len(keysL):
        for mapping in keysL:
            if mapping[-1] == levelCounter:
                orderedKeysL = orderedKeysL + [mapping]
        levelCounter += 1

    lastLevel = orderedKeysL[-1][1]
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


def addScores(treeMin, DTLDict, ParentsDict, ScoreDict, newDTL):
    """Takes the list of reconciliation roots, the DTL , a dictionary of parent nodes, and
    a dictionary of score values, and returns the DTL with the normalized frequency scores calculated."""
    preOrder = preorderDTLsort(DTLDict, treeMin[0][0])
    preOrderCheck = preorderCheck(preOrder)
    for root in preOrderCheck:
        if root != (None, None):
            vertices = root[0]
            if root[1] == 0:
                ParentsDict[vertices] = ScoreDict[vertices]
            for n in range(len(DTLDict[vertices])):
                if type(DTLDict[vertices][n]) == list:
                    child1 = DTLDict[vertices][n][1]
                    child2 = DTLDict[vertices][n][2]
                    oldScore = DTLDict[vertices][n][3]
                    newDTL[vertices][n][3] = ParentsDict[vertices] * (1.0 * oldScore / ScoreDict[vertices])
                    if child1 != (None, None):
                        if child1 in ParentsDict:
                            ParentsDict[DTLDict[vertices][n][1]] += newDTL[vertices][n][3]
                        else:
                            ParentsDict[child1] = newDTL[vertices][n][3]
                    if child2 != (None, None):
                        if child2 in ParentsDict:
                            ParentsDict[child2] += newDTL[vertices][n][3]
                        else:

                            ParentsDict[child2] = newDTL[vertices][n][3]
    normalize = newDTL[preOrderCheck[-1][0]][0][-1]
    for key in newDTL.keys():
        for n in range(len(newDTL[key]) - 1):
            newDTL[key][n][-1] = newDTL[key][n][-1] / normalize

    return newDTL, normalize


def findBest(Parasite, MinimumDict):
    """Takes Parasite Tree and a dictionary of minimum resolution costs and 
    returns a list of the minimum cost reconciliation tree roots"""
    treeMin = {}
    for key in MinimumDict.keys():
        if key[0] == Parasite['pTop'][1]:
            treeMin[key] = MinimumDict[key]

    minimum = treeMin[treeMin.keys()[0]]
    for key in treeMin.keys():
        if treeMin[key] < minimum:
            minimum = treeMin[key]
    for key in treeMin.keys():
        if treeMin[key] > minimum:
            del treeMin[key]
    return treeMin.keys()


def findPath(TupleList, eventDict, uniqueDict):
    """Takes as input TupleList, a list of minimum reconciliation cost roots, eventDict, the dictionary of events and children
     for each node, and uniqueDict, the dictionary of unique vertex mappings. This returns the completed DTL graph as a 
     Dictionary"""
    for Tuple in TupleList:
        if not Tuple in uniqueDict:
            uniqueDict[(Tuple)] = eventDict[(Tuple)]
        for thing1 in eventDict[Tuple]:
            if type(thing1) == list:
                for thing2 in thing1:
                    if type(thing2) == tuple and thing2 != (None, None):
                        findPath([thing2], eventDict, uniqueDict)
    return uniqueDict


def reconcile(fileName, D, T, L):
    """Takes as input a newick file, FileName, a dupliction cost, a transfer cost, and a loss cost. This uses
    newickFormatReader to extract the host tree, parasite tree and tip mapping from the file and then calls 
    DP to return the DTL reconciliation graph of the provided newick file"""
    host, paras, phi = newickFormatReader.getInput(fileName)
    return DP(host, paras, phi, D, T, L)