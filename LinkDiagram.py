r"""

LinkDiagram
===========

AUTHORS:

- Yokota, Hiroshi (2016-03-19): Ver. 0.1.0

ABSTRACT
========

This LinkDiagram is a class of the regular projection of the link
or the knot in S^3. This class needs one list of extended Gauss Codes
for a link or a knot, where the Extended Gauss Code is a list of
crossing points of the link regular porjection.

This class represents the semigroup structure of links for the addition
(i.e. the connected sum), so this class has the additon `+'.

This class needs NetworkX and matplotlib for drawing graphs, 
sqlite3 for storing results to DB(if needed).

ALGORITHMS
==========


Extended Gauss Code::

    -0. The empty list [] is used for the trivial knot component.
    -1. To each components of the link, give orientation.
    -2. Put start points to each components.
    -3. Do the following for each components:
        -a. Move from the start point along to the orientation.
        -b. When you reach a crossing point, if you reach this crossing
            at first time, then you specify your path, upper path(+) or
            lower path(-) as a sign of a crossing point number.
            if you passed this crossing, then you specify the sign of 
           the crossing.

VARIABLES
=========

This class has 6 variables:

    -1. ExGaussCodes:    A list of the Extended Gauss Code
    -2. GaussCodes:      A list of the GaussCodes, this is only path data.
    -3. SeifertCircles:  A list of Seifert Circles for a Seifert Circle
                         of the link.
    -4. Crossings:       A list of the Crossing points data, the sign is
                         of the crossings.
    -5. Components:      The number of the link components.
    -6. SQLite3_DB:      The path to the DB for SQLite3, that is used form
                         a caluculation of Kauffman bracket polynomial.

ETHODS(ABSTRACT)
================
    * __init__():             Initializes the instance.
    * __repr__():             Determins the formal representation
                              of the instance.
    * _add_():                Determins the addition for the connected sum
                              between links.
    * link_component():       Returns the components of the link as a knot.
    * del_link_component():   Returns the links without specified component.
    * mirror_image():         Returns the mirror image of the link(all or 
                                  a specified component).
    * gauss_code():           Returns a GaussCode form an Extend Gauss Code
                                  of the instance. 
    * seifert_circles():      Returns a list of Seifert System.
    * draw_seifert_circles(): Draws the graph of a Seifert System of the link.
                                  This method needs NetworkX and plot from matplotlib.
EXAMPLES
========

DEFINITIONS
***********
    sage: Trivial_Knot = LinkDiagram()
    sage: Trivial_Knot
     GaussCodes:[[]]
     Crossings:[]
     SeifertCircles:[[], []]
     Components: 1

    sage: K3_1 = LinkDiagram([[-1, 3, -2, 1, 3, 2]])
    sage: K3_1
     GaussCodes:[[-1, 3, -2, 1, -3, 2]]
     Crossings:[1, 3, 2]
     SeifertCircles:[[3, 2, 1], [1, 3, 2]]
     Components: 1


ADDITION
*********
    sage: K3_1 = LinkDiagram([[-1,3,-2,1,3,2]])
    sage: K4_1 = LinkDiagram([[-1,4,-2,-1,-3,-2,-4,-3]])
    sage: CS = K3_1 + K4_1        
    sage: CS
     GaussCodes:[[-1, 4, -2, 1, -3, 2, -4, 3]]
     Crossings:[-1, -4, -2, -3]
     SeifertCircles:[[-3, -1], [-2, -1, -4, -3], [-4, -2]]
     Components: 1

HISTORY::
=========

    - 19/03/2016      Version 0.1

"""

#*****************************************************************************
#       Copyright (C) 2016 YOUR NAME  Yokota, Hiroshi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import networkx as nx
import matplotlib.pyplot as plt
from sage.structure.element import FieldElement

class LinkDiagram(FieldElement):
    """
    CLASS: LinkDiagram
    *******************

    VARIABLES::
        - ExGaussCodes:   A list of the extended Gauss Codes for the link.
        - GaussCodes:     A list of the Gauss Codes for the link. 
        - SeifertCircles: A list of the circles of Seifert System for the link.
        - Crossings:      A list of the crossing points with their signs for the link.
        - Components:     An integer, the number of the components for the link.
        - SQLite3_DB:     A string, specifies the path of the database of SQLite3.
    """
    ExGaussCodes = [[]]
    GaussCodes = [[]]
    SeifertCircles = [[]]
    Crossings = []
    Components = 1
    SQLite3_DB = ''
    
    def __init__(self, ExGaussCodes=None):
        """
        ABSTRACT
        ========
        
        Initialization of an instance of LinkDiagram.
        The argument of ExGaussCodes is optionally.
        If no argument, then the instance is for a trivial knot.
        
        For the given link, encode the each components to the Extend Gauss
        Codes. In this initialization, all variables except SQLite3_DB are
        initialized. But SQLite3_DB is specified in calculation of
        a Kauffman bracket polynomial.
        """
        if ExGaussCodes is not None:
            self.ExGaussCodes = ExGaussCodes
            self.GaussCodes = []
            self.Crossings = []
            for ExGaussCode in ExGaussCodes:
                [gssc, self.Crossings] = self.gauss_code(ExGaussCode)
                self.GaussCodes.append(gssc)
            self.Components = len(self.GaussCodes)
            self.SeifertCircles = self.seifert_circles()
            
    def _repr_(self):
        """
        ABSTRACT
        ========
        
        Determin the formal expression of the instance of Link_Diagrams.
        In this expressions, there are all entities except ExGaussCodes
        and SQLite3_DB.
        """
        rstr = ' GaussCodes:%s\n Crossings:%s\n SeifertCircles:%s\n Components: %s\n'
        sgc = str(self.GaussCodes)
        scr = str(self.Crossings)
        ssc = str(self.SeifertCircles)
        scm = str(self.Components)
        return rstr %(sgc, scr, ssc, scm)

    def _add_(self, other):
        """
        ABSTRACT
        ========
        
        Define the addition `+' for LinkDiagram class.
    　    This addition is the connected sum for links.
        For links, the addition is applied on the first components:
        
            x = [ L1,  L2,  ..., Ln]      (x is n components link.)
            y = [ L'1, L'2, ..., L'm]     (y is m components link.)

            x + y = [ L1 + L'1, L2, ..., Ln, L'2, ..., L'm]

        Where L1 + L'1 means the connected sum of a knot L1 and a knot L'1,
        not the addition of LinkDiagram class.
        
        L1 + L'1 is represented with list's `+' between L1 and L2,
        but in this method, I add the procedure that moves the last entity of L2
        to the top of L2, and changes the signs of the crossing corresponding to
        the entity, then adding by `+' with L1. This procedure is more suitable to 
        understand the connected sum of L1 and L2 in the plane.
        
        This addition changes the crossing numbers of links and take
        the connected sum:
        If L1's crossing numbers are from n1 to n2 and L2's are from m1 to m2,
        then The addition changes the crossings n1 to 1, n2 to n2 - n1 + 1,
        m1 to n2 - n1 + 2 and m2 to m2 + n2 -n1 +1. But the operands of `+' are
        not changed.
        
        SYNTAX
        ======
        
        <x> + <y>
        
        INPUT:
        ======

        2 instances of LinkDiagram class.
        
        OUTPUT:
        =======
        
        an instance of LinkDiagram class.
        
        EXAMPLES
        ========
        sage: K3_1 = LinkDiagram([[-1,3,-2,1,3,2]])
        sage: K4_1 = LinkDiagram([[-1,4,-2,-1,-3,-2,-4,-3]])
        sage: CS = K3_1 + K4_1        
        sage: CS
         GaussCodes:[[-1, 4, -2, 1, -3, 2, -4, 3]]
         Crossings:[-1, -4, -2, -3]
         SeifertCircles:[[-3, -1], [-2, -1, -4, -3], [-4, -2]]
         Components: 1  
          
        """
        ExGaussCodes = []
        "The case that the first operand is a trivial knot" 
        if self.ExGaussCodes[0]==[]:
            return(LinkDiagram(other.number_shift(1)))
        else:
            "The case that the second operand is a trivial knot" 
            if other.ExGaussCodes[0]==[]:
                return(LinkDiagram(self.number_shift(1)))
            else:
                fa = map(abs, flatten(self.ExGaussCodes))
                bs = max(fa) - min(fa) + 2
                "1st link starts from 1."
                a = self.number_shift(1)
                "2nd link starts bs."
                b = other.number_shift(bs)
                ExGaussCodes.append(a[0] + b[0])
                for i in a[1:]:
                    ExGaussCodes.append(i)
                for i in b[1:]:
                    ExGaussCodes.append(i)
        return(LinkDiagram(ExGaussCodes))
    
    def number_shift(self, n=None):
        """
        ABSTRACT
        ========
        
        Changes the crossing numbers of the link.
        if n is None, the crossing number starts from 0.
        If n is specified, the crossing number starts from n.
        
        This methods returns the Extended Gauss Codes.
        
        SYNTAX
        ======
        
        obj.number_shift(<n>)
        
        obj.number_shift()
        
        
        INPUT:
        ======
        
        n: an integer>=0 or None.
        
        OUTPUT:
        =======
        
        an Extended Gauss Codes.
        
        """
        ExGaussCodes = []
        a = copy(self.ExGaussCodes)
        m = min(map(abs, flatten(a)))
        if n is not None and n>0:
            bs = m - n
        else:
            bs = m
        for xgausscode in self.ExGaussCodes:
            ExGaussCode = []
            for i in xgausscode:
                ExGaussCode.append(i - sign(i) * bs)
            ExGaussCodes.append(ExGaussCode)
        return(ExGaussCodes)

    def link_component(self, n=None):
        """
        ABSTRACT
        ========
        
        This methods returns the component of the link as a knot.
        If n is None, All components returns as a copy of the instance.
        
        Syntax
        ======
        
        obj.link_component(<n>)
        
        obj.link_component()
        
        INPUT:
        ======
        
        n: an integer>=0 or None. This n specifies the component of the link.
           
        OUTPUT:
        =======
        
        The componet of the link as a LinkDiagram object.
        """
        ExGaussCode = []
        if n is not None:
            k = 0
            if n<self.Components:
                xgssc = self.ExGaussCodes[n]
                axgssc = map(abs, xgssc)
                crssngs = list(set(axgssc))
                chck = map(lambda (x):axgssc.count(x)==2, axgssc)
                for i in chck:
                    if i:
                        ExGaussCode.append(xgssc[k])
                    k = k + 1
        return(Diagram([ExGaussCode]))

    def del_Link_component(self, n=None):
        """
        ABSTRACT
        ========
        
        Returns the link specified component is removed.
        
        SYNTAX
        ======
        
        obj.del_Link_component(<n>)

        obj.del_Link_component()

        INPUT:
        ======
        
        n:  an integer>=0 or None.
        
        OUTPUT:
        =======
        
        A LinkDiagram object, which is removed specified component from the link.
        If no argument, then returns a copy of the object.
        """
        dcrssng = []
        xgsscs = self.ExGaussCodes
        if n is not None:
            if n<self.Components:
                exgssc = xgsscs[n] 
                xgsscs.pop(n)
                if len(exgssc)>0:
                    aexgssc = map(abs,exgssc)
                    chck = map(lambda (x):aexgssc.count(x)==1, aexgssc)
                    k = 0
                    for i in chck:
                        if i:
                            dcrssng.append(aexgssc[k])
                        k = k + 1
                    for i in dcrssng:
                        for x in xgsscs:
                            if i in x:
                                x.remove(i)
                            if -i in x:
                                x.remove(-i)
        return(LinkDiagram(xgsscs))

    def mirror_image(self, n=None):
        """
        ABSTRACT
        ========
        
        Make the mirror image of the specified component.
        If n is None, this method makes the mirror image of all.
        If n>=0, it makes the mirror image of the n-th component only.
        So the crossings between n-th and (n-1)-th, n-th and (n+1)-th
        components are not changed. 
        
        SYNTAX
        ======
        
        obj.mirror_image(<n>)
        
        obj.mirror_image()
        
        INPUT:
        ======
        
        n: an integer>=0 or None.
        
        RETRUN
        ======
        
        None. This method changes the instance itself.
        
        EXAMPLES
        ========
        sage: d2 = LinkDiagram([[-1,3,-2,1,3,-4,5,2],[-4,-5,-6,-6]])
        sage: d2
         GaussCodes:[[-1, 3, -2, 1, -3, -4, 5, 2], [4, -5, -6, 6]]
         Crossings:[1, 3, 2, -4, -5, -6]
         SeifertCircles:[[-6, -4, -5], [-6], [-4, -5, 2, 1, 3], [3, 2, 1]]
         Components: 2
        
        sage: md20 = d2.mirror_iomage(0)
        sage: print md20
         GaussCodes:[[1, -3, 2, -1, 3, -4, 5, -2], [4, -5, -6, 6]]
         Crossings:[-1, -3, -2, -4, -5, -6]
         SeifertCircles:[[-6, -4, -5], [-6], [-4, -5, -2, -1, -3], [-1, -3, -2]]
         Components: 2
         
        sage: md21 = d2.mirror_image(1)
        sage: print md21
         GaussCodes:[[-1, 3, -2, 1, -3, -4, 5, 2], [4, -5, 6, -6]]
         Crossings:[1, 3, 2, -4, -5, 6]
         SeifertCircles:[[6], [-4, -5, 6], [-4, -5, 2, 1, 3], [3, 2, 1]]
         Components: 2
         
        """
        ExGaussCodes = []
        if n is None or n<0 or n>self.Components:
            for x in self.ExGaussCodes:
                ExGaussCodes.append(map(lambda(i):-i, x))
        else:
            xgausscode = map(lambda(z):-z,self.ExGaussCodes[n])
            xcrssngs = map(lambda(x):abs(x), xgausscode)
            crssngs = list(set(xcrssngs))
            cnt = map(lambda(z): xcrssngs.count(z)==1, crssngs)
            w = range(0,len(cnt))
            if sum(cnt)>0:
                crps = list(set(map(lambda(z):cnt[z]*crssngs[z],w)))
                for x in self.ExGaussCodes[0:n]:
                    ExGaussCodes.append(x)
                    for i in crps:
                        if i in xgausscode:
                            n = xgausscode.index(i) 
                            xgausscode[n] = -i
                        else:
                            if -i in xgausscode:
                                n = xgausscode.index(-i)
                                xgausscode[n] = i
                ExGaussCodes.append(xgausscode)
                for x in self.ExGaussCodes[n+1:]:
                    ExGaussCodes.append(x)
                    for i in crps:
                        if i in xgausscode:
                            n = xgausscode.index(i) 
                            xgausscode[n] = -i
                        else:
                            if -i in xgausscode:
                                n = xgausscode.index(-i)
                                xgausscode[n] = i
            else:
                for x in self.ExGaussCodes[0:n-1]:
                    ExGaussCodes.append(x)
                ExGaussCodes.append(map(lambda(z):-z, xgausscode))
                for x in self.ExGaussCodes[n+1:]:
                    ExGaussCodes.append(x)
        return(LinkDiagram(ExGaussCodes))
            
    def gauss_code(self, ExGaussCode):
        """
        ABSTRACT
        ========
        
        Generates the Gauss Code from an Extended Gauss Code.
        This method returns a list of the GaussCode.
        
        SYNTAX
        ======
        
        obj.gauss_code(<ExGaussCode>)

        INPUT:
        ======
        
        ExGaussCode: a list(the extended Gauss Code).
        
        OUTPUT:
        =======
        
        A list(ExGaussCode).
        """
        GaussCode = []
        crssngs = copy(self.Crossings)
        listCrossings = list(set(map(abs, ExGaussCode)))
        for x in ExGaussCode:
            i = abs(x)
            s = sign(x)
            if i in listCrossings:
                listCrossings.remove(i)
                if -i in self.Crossings:
                    crssngs[crssngs.index(-i)] = s*i
                    GaussCode.append(i)
                else:
                    if i in self.Crossings:
                        crssngs[crssngs.index(i)] = s*i
                        GaussCode.append(-i)
                    else:
                        crssngs.append(x)
                        GaussCode.append(x)
            else:
                if x in GaussCode:
                    GaussCode.append(-x)
                else:
                    GaussCode.append(x)
                    if -i in crssngs:
                        crssngs[crssngs.index(-i)] = s*i
                    else:
                        if i in crssngs:
                            crssngs[crssngs.index(i)] = s*i
        return([GaussCode, crssngs])
    
    def seifert_circles(self):
        """
        ABSTRACT
        ========
        
        Generates the Seifert System for the link. 
        Seifert Systrem is a collection of circles,
        which are boundary of a Seifert surface of the link.
        
        This method uses the GaussCodes and the Crossings of the object. 
        This process is a half of the computations for the Kauffman bracket 
        polynomial. Indeed, the skein relation is only used for the path's
        orientation reserving changes. This methods returns a list of
        Seifert cricles, which is contained of the numbers of crossings and
        their signs. 
        
        SYNTAX
        ======
        
        object.seifert_circles()
        
        INPUT:
        ======
        
        None
        
        OUTPUT:
        =======
        
        a list of Seifert Circles( = the Seifert System) for the link.
        
        EXAMPLES
        ========
        sage: d2 = LinkDiagram([[-1,3,-2,1,3,-4,5,2],[-4,-5,-6,-6]])
        sage: d2
         GaussCodes:[[-1, 3, -2, 1, -3, -4, 5, 2], [4, -5, -6, 6]]
         Crossings:[1, 3, 2, -4, -5, -6]
         SeifertCircles:[[-6, -4, -5], [-6], [-4, -5, 2, 1, 3], [3, 2, 1]]
         Components: 2
        """
        GaussCodes = self.GaussCodes
        Crossings = self.Crossings
        newCrossings = copy(Crossings)
        for i in Crossings:
            newGaussCodes = []
            newCrossings.remove(i)
            gcds = []
            gpr = []
            sgn = sign(i)
            ai = abs(i)
            mi = -ai
            for x in GaussCodes:
                if not(i in x or -i in x):
                    gcds.append(x)
                else:
                    gpr.append(x)
            if len(gpr)==1:
                GaussCode = gpr[0]
                p0 = GaussCode.index(mi)
                p1 = GaussCode.index(ai)
                if p0<p1:
                    newGaussCodes.append(GaussCode[0:p0] + [i] + GaussCode[p1+1:])
                    newGaussCodes.append(GaussCode[p0+1:p1] + [i])
                else:
                    newGaussCodes.append([i] + GaussCode[p1+1:p0])
                    newGaussCodes.append(GaussCode[p0+1:] + GaussCode[0:p1] + [i])
            else:
                if len(gpr)==2:
                    if ai in gpr[0]:
                        ga = gpr[1]
                        gb = gpr[0]
                    else:
                        ga = gpr[0]
                        gb = gpr[1]
                    p0 = ga.index(mi)
                    p1 = gb.index(ai)
                    newGaussCodes.append([i] + gb[p1+1:] + gb[0:p1] + [i] + ga[p0+1:] + ga[0:p0])
            for j in gcds:
                newGaussCodes.append(j)
            GaussCodes = newGaussCodes
        return(GaussCodes)

    def draw_seifert_circles(self, file=None, layout=None):
        """
        ABSTRACT
        ========
        
        Generates the directional graph of the Seifert System for the link and
        draws the graph with plot from matplotlib.
        This method needs NetworkX for the directional graph and matplotlib
        for plotting the graph.
        
        SYNTAX
        =======
        
        obj.draw_seifert_circles(file=<file>, layout=<layout>)
        
        obj.draw_seifert_circles(file=<file>)

        obj.draw_seifert_circles(layout=<layout>)

        obj.draw_seifert_circles()

        INPUT:
        ======
        
        Two arguments if needed.
        
        file: Default is None. The filename of the drawing for the Seifert System.
        
        layout: Default is shell. The layout of nodes and edges. Supports circlar, shell
                or spring layout. If None, shell_layout is used.
                
        OUTPUT:
        ========
        
        A directional graph of NetworkX. 
        
        """
        G = nx.MultiDiGraph()
        G.clear()
        scs = self.SeifertCircles
        nodes = []
        nodelist = []
        fsc = []
        edges = []
        a = ''
        b = ''
        asc = 97
        for x in scs:
            b = ''
            tmp = []
            for i in x:
                a = b
                ai = abs(i)
                fsc.append(ai)
                b = chr(asc) + str(ai)
                nodes.append(b)
                tmp.append(b)
                if len(a)>0:
                    edges.append((a,b))
                    G.add_edge(a,b,weight=0.5, sign=0)
            a = chr(asc) + str(abs(x[0]))
            nodelist.append(tmp)
            edges.append((b,a))
            G.add_edge(b,a,weight=0.5,sign=0)
            asc = asc + 1
        G.add_nodes_from(nodes)
        for j in self.Crossings:
            sj = sign(j)
            aj = abs(j)
            p0 = fsc.index(aj)
            p1 = fsc[p0+1:].index(aj) + p0 + 1
            G.add_edge(nodes[p0],nodes[p1],weight= -1, sign=sj)
        "The classifications of edges into 3 types with sign."
        sc=[(u,v) for (u,v,d) in G.edges(data=True) if d['sign']==0]
        pb=[(u,v) for (u,v,d) in G.edges(data=True) if d['sign']==1]
        mb=[(u,v) for (u,v,d) in G.edges(data=True) if d['sign']==-1]
        "The initialization of plt. No plt.clf(), Old drawings are reserved."
        plt.clf()
        if layout is None:
            pos = nx.shell_layout(G)
        else:
            if layout=="circular":
                pos = nx.circular_layout(G)
            else:
                if layout=="spring":
                    pos = nx.spring_layout(G)
                else:
                    pos = nx.shell_layout(G)
        for i in nodelist:
            nx.draw_networkx_nodes(G,pos,nodelist=i,node_size=200,alpha=0.8,color='r')
            nx.draw_networkx_edges(G,pos,edgelist=sc,width=1.5)
        nx.draw_networkx_edges(G,pos,edgelist=pb,style='dashed',edge_color='b',width=1)
        nx.draw_networkx_edges(G,pos,edgelist=mb,style='dotted',edge_color='g',width=1)
        nx.draw_networkx_labels(G,pos,font_size=8,font_family='sans-serif')
        plt.axis('off')
        plt.show()
        if file is not None:
            plt.savefig(file)
        return(G)


def crossing_change_a(connectedDiagram):
    """
    ABSTRACT::
    
       For the crossing change TYPE-A. This crossing change occures for the connected 
　     component only. By this crossing change, A link component are divided two components
　　　 or one component. This splitting is occure by the combination of the Lo or Loo 
       and +1 or -1.
       
    INPUT::
    
        - Diagram:  A list of Gauss code and Corssings.
        
    OUTPUT::
    
        - Diagrams: A list of Diagram.
    
    """
    [GaussCodes, Crossings] = connectedDiagram
    GaussCode = GaussCodes[0]
    n = Crossings[0]
    s = sign(n)
    i = abs(n)
    sa = Crossings[1:]
    sb = copy(sa)
    x = []
    p0 = GaussCode.index(-i)
    p1 = GaussCode.index(i)
    "splitting"
    if p1>p0:
        a1 = GaussCode[p0+1:p1]  
        a2 = GaussCode[p1+1:] + GaussCode[:p0]
    else:
        a1 = GaussCode[p0+1:] + GaussCode[:p1]
        a2 = GaussCode[p1+1:p0]
    "fusion"
    x = a1[-1::-1]
    sb = toggle_crossings(sb, x)
    a0 = [a2 + x]
    pa = [[a1, a2], sa]
    ma = [a0, sb]
    if s==1:
        z = [pa, ma]
    else:
        z = [ma, pa]
    return(z)

def crossing_change_b(disjointDiagram):
    """
    ABSTRACT
    ========
    
        For a crossing between another components.
    
    INPUT
    =====
    
        - Diagram: A list of Gauss code and Crossings for a diagram.
                   This Diagram is not an instance of LinkDiagram.
              
    OUTOUT::
    
        - A list of Diagram.
    
    """
    [GaussCodes, Crossings] = disjointDiagram
    z = []
    n = Crossings[0]
    s = sign(n)
    i = abs(n)
    sa = Crossings[1:]
    sb = copy(sa)
    if -i in GaussCodes[0] and i in GaussCodes[1]:
        CodeA = GaussCodes[0]
        CodeB = GaussCodes[1]
    else:
        CodeA = GaussCodes[1]
        CodeB = GaussCodes[0]
    p0 = CodeA.index(-i)
    p1 = CodeB.index(i)
    px = CodeB[p1+1:] + CodeB[:p1]  
    pa =[[CodeA[0:p0] + px + CodeA[p0+1::]], sa]
    x = px[-1::-1]
    sb = toggle_crossings(sb, x)
    pb =[[CodeA[0:p0] + x + CodeA[p0+1::]], sb]
    if s==1:
        z = [pa, pb]
    else:
        z = [pb, pa]
    return(z)

def crossing_change(Diagram):
    """
    INPUT
    =====
        Diagram:  A list of the Gauss Codes and the Crossings.
                  The Gauss Codes is a list of the link which is derived from 
                  the crossing changes of the link in LinkDiagram class.
                  This `Diagram' is different from LinkDiagram!
                  
    OUTPUT
    ======
        Diagrams: A list of `Diagram's
        
    """
    [GaussCodes, Crossings] = Diagram
    i = abs(Crossings[0])
    y = []
    cpl = []
    newGaussCodes = []
    Diagrams = []
    for x in GaussCodes:
        if not(-i in x or i in x):
            newGaussCodes.append(x)
        else:
            if -i in x and i in x:
                y = crossing_change_a([[x], Crossings])
            else:
                cpl.append(x)
    if len(cpl)>0:
        y = crossing_change_b([cpl, Crossings])
    Diagrams = [[y[0][0] + newGaussCodes, y[0][1]],[y[1][0] + newGaussCodes, y[1][1]]]
    return(Diagrams)

def toggle_crossings(Crossings, GaussCode):
    """
    ABSTRACT::
    
        This function returns a list ot new crossing points with signs after the crossing change.
        
    INPUT::
    
        - Crossings:  A list of the crossing point numbers with signs.
        - GaussCode:  A Gauss code, a list of crossing point numbers.
        
    OUTPUT::
    
        - newCrossings: A list of the crossing numbers with signs after the crossing change.
    """
    newCrossings = []
    list_crossings = list(map(abs, GaussCode))
    ordr = map(abs, Crossings)
    for n in Crossings:
        i = abs(n)
        m = n
        if i in list_crossings:
            if list_crossings.count(i)==1:
                m = -n
        newCrossings.append(m)
    return(newCrossings)



def number2position(n, stage):
    """
    ABSTRACT::
    
        This function returns the position of the link, which is produced with the crossing change.

    INPUT::
    
        - n:      A position of the link in the Diagrams as a list.
        - stage:  The counting of the crossing changes.
    
    OUTPUT::
    
        - Position: A string. This indicates the position of the link in whole crossing changes.
    
    """
    Position = ''
    bg = Integer(n).digits(2)
    m = len(bg)
    for i in bg:
        Position = str(i)+Position
    for j in range(0, stage - m):
        Position = '0' + Position
    return(Position)

def kauffman_bracket(linkDiagram, KnotName=None, DB=None):
    var('A')
    Diagram = [linkDiagram.GaussCodes, linkDiagram.Crossings]
    As = []
    Fs = []
    KBTK = []
    Stage = 0
    if DB is not None and KnotName is not None:
        insert_diagrams2table(KnotName, "diagrams", Stage, int(0), [Diagram])
    Stage = 1
    Crossing = abs(Diagram[1][0])
    Diagrams = crossing_change(Diagram)
    if DB is not None and KnotName is not None:
        insert_diagrams2table(KnotName, "diagrams", Stage, Crossing, Diagrams)
    cps = Diagrams[0][1]
    for i in cps:
        Stage = Stage + 1
        tmp = map(crossing_change, Diagrams)
        Diagrams = []
        for j in tmp:
            Diagrams = Diagrams + j
        Crossing = abs(i)
        if DB is not None and KnotName is not None:
            insert_diagrams2table(KnotName, "diagrams", Stage, Crossing, Diagrams)
    for i in Diagrams:
        j = len(i[0]) - 1
        KBTK.append((-A**2-A**(-2))**j)
    n = len(KBTK)
    m = len(Integer(n-1).digits(2))
    Polynomial = 0
    for i in range(0,n):
        As.append(sum(Integer(i).digits(2)))
    for i in As:
        Fs.append(A**(m-2*i))
    for i in range(0,n):
        Polynomial = Polynomial + KBTK[i]*Fs[i]
    Polynomial = expand(Polynomial)
    if DB is not None and KnotName is not None:
        insert_kauffman_bracket2table(KnotName, "kauffman_bracket_polynomial", Crossing, Diagram, Polynomial)
    return(Polynomial)


def kauffman_bracket_polynomial(LinkDiagram, KnotName=None, DB=None):
    var('A')
    Polynomial = kauffman_bracket(LinkDiagram, KnotName, DB)
    w = Integer(sum(map(sign, LinkDiagram.Crossings)))
    return(expand((-A**3)**(-w)*Polynomial))


def insert_diagrams2table(DBName, TableName, Rolfsen, Stage, CrossingPoint, Diagrams):
    cursor = sqlite3.connect(DBName)
    sql = "insert into " + TableName + "values (?,?,?,?,?,?)"
    m = 0
    for i in Diagrams:
        Position = number2position(m, Stage)
        m = m + 1 
        GaussCodes = str(i[0])
        Crossings = str(i[1])
        cursor.execute(sql,(Rolfsen, int(Stage), int(CrossingPoint), Position, GaussCodes, Crossings))
    cursor.commit()
    cursor.close()
    
def insert_kauffman_bracket2table(DBName, TableName, Rolfsen, CrossingPoint, Diagram, Polynomial):
    cursor = sqlite3.connect(DBName)
    sql = "insert into " + TableName + " values (?,?,?,?,?)",
    GaussCode = str(Diagram[0])
    Crossings = str(Diagram[1])
    cursor.execute(sql, (Rolfsen, int(CrossingPoint), GaussCode, Crossings, str(Polynomial))) 
    cursor.commit()
    cursor.close()






