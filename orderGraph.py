# orderGraph.py
# Carter Slocum
# July 2015

# This file contains a function for topologically ordering a tree graph and detecting cyclic graphs


def date(cycleGraph):
    """takes a CycleCheckingGraph (frankenstien graph compination
        of host and parasite tree)and returns a dictionary representation
        of the ordering of the tree. If there is a cycle, the function returns timeTravel"""
    # ordering of the Nodes
    order = {}
    # InnerNodes
    innerNodes = {}
    # dict of Leaves
    Leaves = {}
    # List of Nodes with In-degree zero
    LonerList = []
    for key in cycleGraph.keys():
        if key is not None:
            innerNodes[key] = 0
        for child in cycleGraph[key]:
            if child is not None:
                innerNodes[child] = 0
            else:
                Leaves[key] = True
    for key in Leaves.keys():
        if key in innerNodes.keys():
            del innerNodes[key]
    for key in innerNodes.keys():
        for child in cycleGraph[key]:
            if child is not None and child in innerNodes:
                innerNodes[child] += 1
    place = 0
    for key in innerNodes.keys():
        if innerNodes[key] == 0:
            del innerNodes[key]
            LonerList.append(key)
    while LonerList:
        nodeZero = LonerList[0]
        del LonerList[0]
        order[nodeZero] = place
        place += 1
        for child in cycleGraph[nodeZero]:
            if child in innerNodes.keys():
                innerNodes[child] -= 1
                if innerNodes[child] == 0:
                    del innerNodes[child]
                    LonerList.append(child)
    if len(innerNodes.keys()) > 0:
        return "timeTravel"
    else:
        for leaf in Leaves:
            order[leaf] = len(Leaves)
        return order