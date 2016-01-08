# newickToVis.py
#July 2015
#Carter Slocum

#File contains function that creates separate newick files for the 
#parasite tree and the ultra-metric host tree.

from compbio.rasmus import treelib1, util
import copy
import sys
import cycleCheckingGraph
import newickFormatReader


def convert(fileName, HostOrder, n, writeParasite):
    """takes name of original .newick file and the dictionary of host tree branch lengths
    and creates files for the host + parasite trees. Parasite tree can
    be ommited if desired"""
    f = open(fileName, 'r')
    contents = f.read()
    host, paras, phi = newickFormatReader.getInput(fileName)
    hostRoot = cycleCheckingGraph.findRoot(host)
    f.close()
    H, P, phi = contents.split(";")
    P = P.strip()
    H = H.strip()
    H = H + ';'
    host = treelib1.parse_newick(H, HostOrder)
    for key in HostOrder:
        H = H.replace(str(key), str(key) + ':' + str(HostOrder[key]))
    f = open(fileName[:-7] + str(n) + ".stree", 'w')
    treelib1.write_newick(host, f, root_data=True)
    f.close()
    if writeParasite:
        f = open(fileName[:-7] + '.tree', 'w')
        f.write(P + ";")
        f.close()