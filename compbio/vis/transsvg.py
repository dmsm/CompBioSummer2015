# transsvg.py
# Edited by Annalise Schweickart
# July 2015

# This is the drawing tool of vistrans, drawing the host and parasite trees and the events that occur
# This code was edited to include loss events, transfer events causing cycles, and frequencies of events
from collections import defaultdict

from rasmus import treelib1, svg, util, stats
from compbio import phylo
import copy

def draw_stree(canvas, stree, slayout,
               yscale=100,
               stree_width=.8, 
               stree_color=(.4, .4, 1),
               snode_color=(.2, .2, .4)):
    '''Takes as input a canvas, a host tree, stree, created by treelib1, a slayout, also created by
    treelib1, a scaling, yscale, the width of the tree, sree_width, and the colors of the tree and its nodes
    and draws the host tree on which the species tree is built'''
    # draw stree branches
    for node in stree:
        x, y = slayout[node]
        px, py = slayout[node.parent]
        
        # draw branch
        w = yscale * stree_width / 2.0
        canvas.polygon([px, py-w,
                        x, y-w,
                        x, y+w,
                        px, py+w], strokeColor=(0, 0, 0, 0),
                       fillColor=stree_color)

    # draw stree nodes
    for node in stree:
        if not node.is_leaf():
            x, y = slayout[node]
            canvas.line(x, y-w, x, y+w, color=snode_color,
                        style="stroke-dasharray: 1, 1")
    
def draw_tree(tree, brecon, stree,
              xscale=100, yscale=100,
              leaf_padding=10, 
              label_size=None,
              label_offset=None,
              font_size=12,
              stree_font_size=20,
              canvas=None, autoclose=True,
              rmargin=10, lmargin=10, tmargin=0, bmargin=0,
              tree_color=(0, 0, 0),
              tree_trans_color=(0, 0, 0),
              stree_color=(.3, .7, .3),
              snode_color=(.2, .2, .7),
              loss_color = (1,1,1),
              loss_color_border=(.5,.5,.5),
              dup_color=(0, 0, 1),
              dup_color_border=(0, 0, 1),
              trans_color=(1, 1, 0),
              trans_color_border=(.5, .5, 0),
              gtrans_color=(1,0,0),
              gtrans_color_border=(.5,0,0),
              event_size=10,
              snames=None,
              rootlen=None,
              stree_width=.8,
              filename="tree.svg"
              ):
    '''Takes as input a parasite tree, tree, a reconciliation file, brecon, a host tree, stree, as well as
    sizes and colors of the trees components and returns a drawing of the reconciliation of the parasite 
    tree on the host tree with event nodes of specified colors'''
    # set defaults
    font_ratio = 8. / 11.    
    if label_size is None:
        label_size = .7 * font_size


    if sum(x.dist for x in tree.nodes.values()) == 0:
        legend_scale = False
        minlen = xscale

    if snames is None:
        snames = dict((x, x) for x in stree.leaf_names())

    # layout stree
    slayout = treelib1.layout_tree(stree, xscale, yscale)
    if rootlen is None:
        rootlen = .1 * max(l[0] for l in slayout.values())

    # setup slayout
    x, y = slayout[stree.root]
    slayout[None] =  (x - rootlen, y)
    for node, (x, y) in slayout.items():
        slayout[node] = (x + rootlen, y  - .5 * yscale)

    # layout tree
    ylists = defaultdict(lambda: [])
    yorders = {}
    # layout speciations and genes (y)
    #losses = copy.copy(loss)
    for node in tree.preorder():
        yorders[node] = []
        for ev in brecon[node]:
            snode, event, frequency = ev
            if event == "spec" or event == "gene" or event == "loss":
                yorders[node].append(len(ylists[snode]))
                ylists[snode].append(node)


    # layout dups and transfers (y)
    for node in tree.postorder():

        for ev in brecon[node]:
            snode, event, frequency = ev
            if event != "spec" and event != "gene" and event != "loss":
                v = [yorders[child]
                    for child in node.children
                    if brecon[child][-1][0] == snode]
                if len(v) == 0:
                    yorders[node].append(0)
                else:
                    yorders[node].append(stats.mean(flatten(v)))

    # layout node (x)
    xorders = {}
    xmax = defaultdict(lambda: 0)
    #losses = copy.copy(loss)
    for node in tree.postorder():
        xorders[node] = []
        for ev in brecon[node]:
            snode, event, frequency = ev
            if event == "spec" or event == "gene" or event == "loss":
                    xorders[node].append(0)


            else:
                v = [xorders[child] for child in node.children
                    if brecon[child][-1][0] == snode]
                if len(v) == 0:
                    xorders[node].append(1)
                else:
                    xorders[node].append(maxList(v) + 1)
            xmax[snode] = max(xmax[snode], maxList([xorders[node]]))

    # setup layout
    layout = {None: [slayout[brecon[tree.root][-1][0].parent]]}
    for node in tree:
        for n in range(len(brecon[node])):
            snode, event, frequency = brecon[node][n]
            nx, ny = slayout[snode]
            px, py = slayout[snode.parent]


        # calc x
            frac = (xorders[node][n]) / float(xmax[snode] + 1)
            deltax = nx - px
            x = nx - frac * deltax

        # calc y
            deltay = ny - py
            slope = deltay / float(deltax)
            deltax2 = x - px
            deltay2 = slope * deltax2
            offset = py + deltay2
            frac = (yorders[node][n] + 1) / float(max(len(ylists[snode]), 1) + 1)
            y = offset + (frac - .5) * stree_width * yscale

            if node in layout: layout[node].append((x, y))
            else:
                layout[node] = [(x, y)]
        
        # order brecon nodes temporally
        brecon[node] = orderLoss(node, brecon, layout)
        # order layout nodes temporally
        layout[node] = orderLayout(node, layout)

        if y > max(l[1] for l in slayout.values()) + 50:
            print nx, ny
            print px, py
            print offset, frac
            print ylists[snode], yorders[node]
            print brecon[node]
            print node, snode, layout[node]


    # layout label sizes
    max_label_size = max(len(x.name)
        for x in tree.leaves()) * font_ratio * font_size
    max_slabel_size = max(len(x.name)
        for x in stree.leaves()) * font_ratio * stree_font_size


    '''
    if colormap == None:
        for node in tree:
            node.color = (0, 0, 0)
    else:
        colormap(tree)
    
    if stree and gene2species:
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.label_events(tree, recon)
        losses = phylo.find_loss(tree, stree, recon)
    else:
        events = None
        losses = None
    
    # layout tree
    if layout is None:
        coords = treelib.layout_tree(tree, xscale, yscale, minlen, maxlen)
    else:
        coords = layout
    '''
    
    xcoords, ycoords = zip(* slayout.values())
    maxwidth = max(xcoords) + max_label_size + max_slabel_size
    maxheight = max(ycoords) + .5 * yscale
    
    
    # initialize canvas
    if canvas is None:
        canvas = svg.Svg(util.open_stream(filename, "w"))
        width = int(rmargin + maxwidth + lmargin)
        height = int(tmargin + maxheight + bmargin)
        
        canvas.beginSvg(width, height)
        canvas.beginStyle("font-family: \"Sans\";")
        
        if autoclose == None:
            autoclose = True
    else:
        if autoclose == None:
            autoclose = False

    canvas.beginTransform(("translate", lmargin, tmargin))
    
    draw_stree(canvas, stree, slayout,
               yscale=yscale,
               stree_width=stree_width, 
               stree_color=stree_color,
               snode_color=snode_color)

    # draw stree leaves
    for node in stree:
        x, y = slayout[node]
        if node.is_leaf():
            canvas.text(snames[node.name], 
                        x + leaf_padding + max_label_size,
                        y+stree_font_size/2., stree_font_size,
                        fillColor=snode_color)


    # draw tree

    for node in tree:

        containsL= containsLoss(node, brecon)
        for n in range(len(brecon[node])):
            x, y = layout[node][n]
            
            if containsL == False: # no loss event
                px, py = layout[node.parent][-1]       
            else: # loss event present
                if n == 0: # event is loss
                    px, py = layout[node.parent][-1]
                else: # event stems from loss
                    px, py = layout[node][n-1]
            

            trans = False

            if node.parent:
                snode, event, frequency =  brecon[node][n]
                if n == 0:
                    psnode, pevent, pfrequency = brecon[node.parent][-1]
                
                # Event stemming from a loss event
                else: psnode, pevent, pfrequency = brecon[node][n-1]
                if pevent == "trans" or pevent == "gtrans":
                    if psnode != snode:
                        trans = True
                else: trans = False

                if not trans:
                    canvas.line(x, y, px, py, color=tree_color)
                
                # draw the transfer dashed line        
                else:
                    arch = 20
                    x2 = (x*.5 + px*.5) - arch
                    y2 = (y*.5 + py*.5)
                    x3 = (x*.5 + px*.5) - arch
                    y3 = (y*.5 + py*.5)
                    
                    if pevent == "trans":
                        canvas.write("<path d='M%f %f C%f %f %f %f %f %f' %s />\n " %
                            (x, y, x2, y2,
                             x3, y3, px, py,
                            " style='stroke-dasharray: 4, 2' " +
                            svg.colorFields(tree_trans_color, (0,0,0,0))))
                    # draw guilty transfer line
                    else: canvas.write("<path d='M%f %f C%f %f %f %f %f %f' %s />\n " %
                            (x, y, x2, y2,
                             x3, y3, px, py,
                            " style='stroke-dasharray: 4, 2' " +
                            svg.colorFields(gtrans_color, (0,0,0,0))))


    # draw events
    for node in tree:
        for n in range(len(brecon[node])):
            snode, event, frequency =  brecon[node][n]
            x, y = layout[node][n]
            o = event_size / 2.0
            if event == "loss":
                canvas.rect(x - o, y - o, event_size, event_size,
                        fillColor=loss_color,
                        strokeColor=loss_color_border)
                canvas.text(frequency, x-o, y-o, font_size+2, fillColor = (1,1,1))

	
            if event == "spec":
                canvas.rect(x - o, y - o, event_size, event_size,
                        fillColor=(0,0,0),
                        strokeColor=(0,0,0))
                canvas.text(frequency, x-o, y-o, font_size+2, fillColor = (0,0,0))


            if event == "dup":
                canvas.rect(x - o, y - o, event_size, event_size,
                        fillColor=dup_color,
                        strokeColor=dup_color_border)
                canvas.text(frequency, x-o, y-o, font_size+2, fillColor=dup_color)

            elif event == "trans":
                canvas.rect(x - o, y - o, event_size, event_size,
                        fillColor=trans_color,
                        strokeColor=trans_color_border)
                canvas.text(frequency, x-o, y-o, font_size+2, fillColor=trans_color)
            
            elif event == "gtrans":
                canvas.rect(x-o, y-o, event_size, event_size,
                        fillColor=gtrans_color,
                        strokeColor=gtrans_color_border)
                canvas.text(frequency, x-o, y-o, font_size+2, fillColor=gtrans_color)

    # draw tree leaves
    for node in tree:
        for n in range(len(brecon[node])):
            x, y = layout[node][n]
            if node.is_leaf() and brecon[node][n][1] == "gene":
                canvas.text(node.name, 
                        x + leaf_padding, y+font_size/2., font_size+2,
                        fillColor=(0, 0, 0))

        
    canvas.endTransform()
    
    if autoclose:
        canvas.endStyle()
        canvas.endSvg()
    
    return canvas

def containsLoss(node, brecon):
    '''Takes as input a node and a reconciliation, brecon, and returns a boolean
    of whether that node leads to a loss event or not'''
    containsL = False
    for n in range(len(brecon[node])):
        if brecon[node][n][1] == "loss":
            containsL = True
    return containsL


def orderLoss(node, brecon, layout):
    '''Takes as input a node, a reconciliation, and a tree layout and changes the
    nodes reconciliation to be ordered based on location in the layout. The nodes
    occuring earlier in the graph occur earlier in the nodes reconciliation entry'''
    eventDict = {}
    for n in range(len(brecon[node])):
        eventDict[str(brecon[node][n])] = layout[node][n]
    return sorted(brecon[node], key=lambda xCoord: eventDict[str(xCoord)][0])

def orderLayout(node, layout):
    '''Takes as input a node and a tree layout and returns the layout of that
    nodes children in timed order, much like orderLoss'''
    return sorted(layout[node], key=lambda xCoord: xCoord[0])

def maxList(vList):
    '''Takes as input vList, a list of lists, and finds the maximum value in all
    of the lists'''
    maxVal = 0
    for item in vList:
        for value in item:
            if value >= maxVal:
                maxVal = value
    return maxVal

def flatten(vList):
    '''Takes as input a list of lists, vList, and returns a newList containing every item in the
    internal lists of vList'''
    newList = []
    for item in vList:
        newList.extend(item)
    return newList
