#!/usr/bin/env sage -python

import itertools
from sage.all import *

class MyGraph:
     
    def sort_leaves(self):
        """This is a helper that takes a graph G and returns tuple (L, I) where L
        is a list of the leaves and I is a list of the non-leaf vertices."""
        leaves = []
        nonleaves = []

        for v in self.graph.vertices():
            if len(self.graph.neighbors(v)) == 1:
                leaves.append(v)
            else:
                nonleaves.append(v)
    
        return (sorted(leaves, key=str), sorted(nonleaves, key=str))
    
    def make_ys(self, y = "y"):
        """This method makes the y-variables that correspond go edges in the graph."""
        (L, I) = (self.leaves, self.nonleaves)
    
        ys = {}
        for e in self.graph.edges(labels=False):
            (u, v) = e
            
            if u in I and v in I:
                u = str(u).split("_")[1]
                v = str(v).split("_")[1]
                
                if u < v:
                    edgename = y + "_" + u + v
                else:
                    edgename = y + "_" + v + u

            elif u in L:
                u = str(u).split("_")[1]
                edgename = y + "_" + str(u)
            else:
                v = str(v).split("_")[1]
                edgename = y + "_" + str(v)
           
            #ys[var(edgename)] = e
            ys[e] = var(edgename)

        self.ys_map  = ys
        self.ys_list = sorted(self.ys_map.values(), key=str)
        #self.ys_list = sorted(self.ys_map.keys(), key=str)

    def make_qs(self, q = "q", s = "s"):
        """This helper method creates the list of qs and the mapping from ys to qs."""
        #P = sorted(powerset(self.leaves), key=len)
        P = powerset(self.leaves)
        n = self.last_taxon
    
        self.qs_list = []
        self.ss_list = []

        for S in P:
            if n not in S:
                sorted_S = sorted(S, key=str)
                name = "".join(map(str, sorted_S))
                if name == "": 
                    name = "0"
                new_q = q + "_" + name
                new_s = s + "_" + name
                self.qs_list.append(var(new_q))
                self.ss_list.append(var(new_s))

        self.ss_map = dict(zip(self.ss_list, self.qs_list))
    
        self.qs_map = {}
        for (e, y) in self.ys_map.items():
            # take out the edge temporarily
            self.graph.delete_edge(e)
            [C, C2] = self.graph.connected_components()

            if n in C: C = C2
            
            sorted_C = sorted(filter(lambda x: x in self.leaves, C), key=str)   
            matching_q = q + "_" + "".join(map(str, sorted_C)) 
            
            self.qs_map[y] = var(matching_q) 
            
            # put the edge back into the graph
            self.graph.add_edge(e)
        
        # print sorted(self.qs_map.items(), key=str)

    def __init__(self, name, E):
        """An instance of a MyGraph object is initalized by passing a name and
        an edge set. The graph will then initialize the following data members:
            self.graph          self.ys_map
            self.name           self.ys_list
            self.last_taxon     self.qs_map
            self.leaves         self.qs_list
            self.nonleaves
        Then it will also register itself with the global dictionary my_graphs."""

        #self.graph = Graph(map(lambda (x,y): (var(x), var(y)), E),)
        
        # making an effort to not overwrite e this time
        edges = map(lambda (x,y): (var("v_" + str(x)), var("v_" + str(y))), E)
        
        self.graph = Graph(edges)
        self.name  = name
        self.taxa  = sorted(self.graph.vertices(), key=str)
        self.graph.name(name)


        (self.leaves, self.nonleaves) = self.sort_leaves()
        self.last_taxon = sorted(self.leaves, key=str, reverse=True)[0]
       
        # make other edge/split representation things
        self.make_ys()
        # self.make_qs()
    
    def __str__(self):
        return self.graph.__str__()

    def characters(self):
        P = powerset(self.leaves)
        n = self.last_taxon

        chars = []
        for S in P:
            if n not in S:
                vector = [int(x in S) for x in self.leaves]
                chars.append(zip(self.leaves, vector))

        return chars

    def compatible_chars(self):
        """This is a method that returns the compatible characters of a graph."""
        L, I = self.leaves, self.nonleaves 
        
        # never trust python dictionaries
        chars, labels = [], []
        
        for (e, y) in sorted(self.ys_map.items(), key=lambda x: str(x[1])):
            char, labeling = [], []

            self.graph.delete_edge(e)
            components = self.graph.connected_components()
            [c1, c2] = components
            
            for v in c1:
                if v in L: 
                    char += [(v,0)]
                else:
                    labeling += [(v,0)]

            for v in c2:
                if v in L: 
                    char += [(v,1)]
                else:
                    labeling += [(v,1)]

            self.graph.add_edge(e)
            chars.append(char)
            labels.append(labeling)

            #print char, labeling

        s = lambda x: str(x[0])

        chars  = [zip(L, [0]*len(L))] + map(lambda char: sorted(char, key=s), chars)
        labels = [zip(I, [0]*len(I))] + map(lambda label: sorted(label, key=s), labels)
        
        return chars, labels

def labelings(n):
    """
    This is a terrible example of me trying to code iterated Cartesian product
    at 2 am and getting it to a working state. It generates the internal labelings
    for n nodes. 

    For example, labelings(2) returns [(0,0), (0,1), (1,0), (1,1)]."""
    b = [0, 1]

    # this is the most ridiculous thing ever.
    result = itertools.product(b, b)
    
    if n == 1:
        return map(lambda x: [x], b)
    else:
        for i in range(0, n-2):
            result = itertools.product(list(result), b)

    L = map(flatten , list(result))
    return L

def labeled_likelihood(G, char, label):

    labeled_char = dict(label + char)

    p = 1
    
    for (e, y) in sorted(G.ys_map.items(), key=str):
        (u, v) = e

        if labeled_char[u] == labeled_char[v]:
            p = p * (1 + y)
        else:
            p = p * (1 - y)

    return p


def average_likelihood(G, data, normalizing=False):
    """Given a graph G and character char, this function returns the likelihood
    function of char on G.""" 

    (L, I) = (G.leaves, G.nonleaves)
    labels = map( lambda x: zip(I, x), labelings(len(I)))
    
    print labels
    
    chars = G.characters()

    if len(data) != len(chars):
        print data, "does not have all characters; trying compatible characters."
        
        c, _ = G.compatible_chars()
        if len(data)==len(c):
            print "Using compatible chars."
            chars = c
        else:
            print "That didn't work."
            return
    
    avg_likelihood = 1 
    
    for i in range(len(data)):
        exp, char = data[i], chars[i]
        
        if exp == 0: continue
        
        char_avg_likelihood = 0

        for label in labels:
            char_avg_likelihood += labeled_likelihood(G, char, label)
            
        
        if normalizing: char_avg_likelihood /= 2

        #print char_avg_likelihood.expand()
        #print
        
        avg_likelihood *= char_avg_likelihood**exp

    #return avg_likelihood.canonicalize_radical() 
    return avg_likelihood

    # wtf sage?
    # f = f.canonicalize_radical()



