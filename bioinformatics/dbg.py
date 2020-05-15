class DeBruijnGraph:
    ''' De Bruijn directed multigraph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes
        are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. '''

    @staticmethod
    def chop(st, k):
        ''' Chop string into k-mers of given length '''
        for i in range( len( st ) - (k - 1) ):
            yield (st[i:i + k], st[i:i + k - 1], st[i + 1:i + k])

    class Node:
        ''' Node representing a k-1 mer.  Keep track of # of
            incoming/outgoing edges so it's easy to check for
            balanced, semi-balanced. '''

        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0

        def isSemiBalanced(self):
            return abs( self.nin - self.nout ) == 1

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash( self.km1mer )

        def __str__(self):
            return self.km1mer

    def __init__(self, strIter, k, circularize=False):
        self.G = {}  # multimap from nodes to neighbors
        self.nodes = {}
        for st in strIter:
            if circularize:
                st += st[:k - 1]
            for kmer, km1L, km1R in self.chop( st, k ):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node( km1L )
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node( km1R )
                nodeL.nout += 1
                nodeR.nin += 1
                self.G.setdefault( nodeL, [] ).append( nodeR )
        # Iterate over nodes; tally # balanced, semi-balanced, neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian walk (not cycle)
        self.head, self.tail = None, None
        for node in iter( self.nodes.values() ):
            if node.isBalanced():
                self.nbal += 1
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def nnodes(self):
        ''' Return # nodes '''
        return len( self.nodes )

    def nedges(self):
        ''' Return # edges '''
        return len( self.G )

    def hasEulerianWalk(self):
        ''' Return true iff graph has Eulerian walk. '''
        return self.nneither == 0 and self.nsemi == 2

    def hasEulerianCycle(self):
        ''' Return true iff graph has Eulerian cycle. '''
        return self.nneither == 0 and self.nsemi == 0

    def isEulerian(self):
        ''' Return true iff graph has Eulerian walk or cycle '''
        # technically, if it has an Eulerian walk
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    def eulerianWalkOrCycle(self):
        ''' Find and return sequence of nodes (represented by
            their k-1-mer labels) corresponding to Eulerian walk
            or cycle '''
        #assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            g.setdefault( self.tail, [] ).append( self.head )
        # graph g has an Eulerian cycle
        tour = []
        src = next( iter( g.keys() ) )  # pick arbitrary starting node

        def __visit(n):
            while len( g[n] ) > 0:
                dst = g[n].pop()
                __visit( dst )
            tour.append( n )

        __visit( src )
        tour = tour[::-1][:-1]  # reverse and then take all but last node

        if self.hasEulerianWalk():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index( self.head )
            tour = tour[sti:] + tour[:sti]

        # Return node list
        return list( map( str, tour ) )

class DeBruijnPlot(DeBruijnGraph):
    def to_dot(self, weights=False):
        """ Write dot representation to given filehandle.  If 'weights'
            is true, label edges corresponding to distinct k-1-mers
            with weights, instead of writing a separate edge for each
            copy of a k-1-mer. """
        dot_str = []
        dot_str.append("digraph \"DeBruijn graph\" {\n")
        dot_str.append("  bgcolor=\"transparent\";\n")
        for node in self.G.keys():
            lab = node.km1mer
            dot_str.append("  %s [label=\"%s\"] ;\n" % (lab, lab))
        for src, dsts in self.G.items():
            srclab = src.km1mer
            if weights:
                weightmap = {}
                if weights:
                    for dst in dsts:
                        weightmap[dst] = weightmap.get(dst, 0) + 1
                for dst, v in weightmap.iteritems():
                    dstlab = dst.km1mer
                    dot_str.append("  %s -> %s [label=\"%d\"] ;\n" % (srclab, dstlab, v))
            else:
                for dst in dsts:
                    srclab = src.km1mer
                    dstlab = dst.km1mer
                    dot_str.append("  %s -> %s [label=\"\"] ;\n" % (srclab, dstlab))
        dot_str.append("}\n")
        return ''.join(dot_str)
