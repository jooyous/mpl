from graph import MyGraph
from graph import likelihood    


def generate_partials(G, chars=None, normalizing=False):
    """This function takes a graph G and an optional argument chars and
    generates the output (T, ps, exps), where T is the product of arbitrary
    powers of the likelihood functions for all the characters in chars, and ps
    is the list of partials of T with respect to all of its edges and exps is a
    list of the arbitrary exponent variables used in T.  If chars is
    unspecified, the function uses the compatible characters of G."""
    i = 0
    T = 1

    exps  = []
    if chars == None: 
        print "No chars given. Using compatible characters of " + g.name
        chars = G.compatible_chars()
    
    print "Generating likelihood function T."
    for x in chars:
        exp = "k" + str(i)
        var(exp)
        exps.append(exp)
    
        f = likelihood(G, x, normalizing)
    
        T = T * f^var(exp)
        i = i + 1
    
    print "Computing partials."
    ps = []
    for y in G.ys_list:
        print "  partial with respect to " + str(y)
        ps.append(T.diff(y) == 0)
    
    return (T, ps, exps)

def analyze_partials(ps, ks):
    """This function constructs and row-reduces the matrix of partials
    generated by the generate_partials function with respect to the exps."""

    M = []

    print "Building matrix."
    for p in ps:
        row = []
        for k in ks:
            row.append(p.lhs().coefficient(var(k)))
        M.append(row)

    M = matrix(M)
    print type(M)
    print "Dimensions: " + str(M.dimensions())

    #print "Simplifying."
    #M = M.factor()

    print "Row-reducing."
    MM = M.echelon_form()
    
    print "Simplifying."
    MM = MM.factor()

    return M

#def second_point(G, M):
#       """This function generates a new set of edgenames for G, generates the equation matrix
#       for the two 'distinct' critical points and returns the row-reduced matrix."""
#       edges = my_edgenames[G.name()]
#       other_edges = make_edgenames(G, "z")
#       
#       print "Building new matrix."
#   
#       D = {}
#   
#       for key in edges:
#           D[var(edges[key])] = var(other_edges[key])
#       
#       MM = map(list, M.rows())
#       Mz = M.substitute(D)
#       
#       for row in Mz.rows():
#           MM.append(list(row))
#       
#       MM = matrix(MM)
#   
#       print MM.dimensions()
#       print "Row-reducing."
#       return MM.echelon_form()

def build_exponents(G, point, M = None, normalizing=False, chars=None):
    """M needs to be the output of corrected_matrix or it will be wrong."""
    if chars == None: chars = G.compatible_chars()
    if M == None: 
        M = corrected_matrix(G, chars)

    n = normalizing

    likes = [likelihood(G, char, normalizing=n) for char in chars]
    ks    = [var("k" + str(i)) for i in range(0, len(likes))]

    #p = dict(zip(G.ys_list, list(point))) 
    #p = Point(g, point) 
    
    print "Substituting point: " + str(point.ys)
    M = M.subs(point.ys)
    print M
    M = M.rref()
    print M

    eqs = []
    for row in M:
        row = list(row)
        eq = 0
        for i in range(0, len(row)): 
            eq = eq + row[i] * var(ks[i])
        
        eqs.append(eq.factor())
   
    if M.nrows() < len(likes): 
        overlap = M.nrows() - len(likes)
    else:
        overlap = M.nrows()
    
    solns = solve(eqs, ks[:overlap])
    choices = {x:1 for x in ks[overlap:]}

    expvalues  = []

    for soln in solns:
        soln = list(soln)
       
        for i in range(0, len(soln)):
            exp = soln[i].rhs().simplify_exp() 
            expvalues.append(exp)
            
            for (kk, v) in choices.items():
                d = exp.coefficient(var(kk)).denominator()
                choices[kk] = lcm(d, v)
       
        print "Solution to exponent equations: " + str(expvalues)
        print "Reasonable choices for free variables " + str(choices)
        
        ks = [k.subs(choices) for k in expvalues] + [choices[k] for k in ks[overlap:]]
        print "ks vector: " + str(ks)

        print "Computing likelihood function."
        
        T = 1
        for i in range(0, len(ks)):
            T = T * likes[i]^ks[i]
        
        #for i in range(0, len(choices)):
        #     T = T * likes[i + M.nrows()]^var(ks[i + M.nrows()])

        T = T.function(*G.ys_list)

        print "Computing the Hessian for T at given point."
        hess = T.diff(2)(*point.y_vector)
       
        # check to see if p is a rational point to choose the ring to switch to
        #if reduce(operator.and_, map(lambda x: x in QQ, list(p))):
        #    print "  Changing ring to QQ."
        #    hess = hess.change_ring(QQ)
        #else:
        #    print "  Changing ring to RR."
        #    hess = hess.change_ring(RR)

        print "Computing eigenvalues."
        eigs = hess.eigenvalues()
        
        eigstring = ""
        for eig in eigs:
            if    eig < 0: eigstring += "-"
            elif  eig > 0: eigstring += "+"
            else: eigstring += "0"
        
        print eigstring
        return T

def fix_matrix(M):
    rows = M.nrows()
    L = [x.simplify_exp() for x in M.list()]
    return matrix(rows, L)


def correction_matrix(G, chars=None, normalizing=False):

    if chars == None:
        chars = G.compatible_chars()

    likes = [likelihood(G, char, normalizing) for char in chars] 
    return diagonal_matrix(map(lambda x: 1/x, likes))

def corrected_matrix(G, chars=None, normalizing=False):
    
    M = chargrad_matrix(G, chars, normalizing)

    CM = M * correction_matrix(g, chars, normalizing) 

    print "Simplifying corrected matrix."
    return CM.factor()
    
def chargrad_matrix(G, chars=None, normalizing=False):
    """Because it's a handmade matrix! Hahahaha."""
    (L, I) = (G.leaves, G.nonleaves)

    if chars == None:
        print "No chars given. Using compatible characters of " + g.name
        chars = G.compatible_chars()

    likes = {}
    for char in chars:
        likes[str(char)] = likelihood(G, char, normalizing)

    M = []
    for edge in G.ys_list:
        row = []
        for char in chars:
            F  = likes[str(char)]
            dF = F.diff(var(edge))
            row.append(dF)
        M += [row]

    M = matrix(M)
    #print "Row reducing."
    #MM = M.echelon_form()
    # print "Simplifying."
    # M = M.factor()
    return M 
    

# this is a global dictionary of graphs
my_graphs = {}

# These are node names; might have to add extra ones.
# IMPORTANT: DON'T OVERWRITE e
# var("a b c d u v w x y z")

# Some trees we've looked at before.
my_graphs["tree3"] = MyGraph("tree3", [('a','u'), ('b','u'), ('c','u')])
my_graphs["tree4"] = MyGraph("tree4", [('a','u'), ('b','u'), ('u','v'), ('v','c'), ('v','d')])


class Point:
    def __init__(self, graph, point):
        """This method initializes a point by passing the graph it corresponds
        to and a list or tuple of values. The class figures out if point is in
        y-space or split space."""
       
        self.graph = graph

        num_splits = 2^(len(graph.leaves)-1)
        H = hadamard_matrix(num_splits)
        point = list(point)

        if len(point) == len(graph.ys_list):
            print "Point in edge space given."
            
            self.y_vector = vector(point)
            self.ys = dict(zip(graph.ys_list, list(point))) 
            
            self.qs = {q:0 for q in self.graph.qs_list}
           
            # this is taken from Rob's cheat sheet on sharelatex
            for (yy, qq) in graph.qs_map.items():
                self.qs[qq] = -1/2 * ln(self.ys[yy])
                self.qs[qq] = self.qs[qq].simplify_full()
           
            # fix up the empty-set entry to make the Hadamard stuff work out.
            if var("q_0") in self.qs.keys():
                self.qs[var("q_0")] = 1 - sum(self.qs.values())
            else:
                print "Missing key q_0. Did you use a different letter?"
                return
            
            self.q_vector = vector([self.qs[q] for q in graph.qs_list])

            self.s_vector = H.inverse() * vector([exp(x) for x in (H * self.q_vector)])
            self.s_vector = vector([x.simplify_full() for x in self.s_vector])

            self.ss = dict(zip(graph.ss_list, self.s_vector))

        elif len(list(point)) == num_splits:
            print "Point in split space given."
            
            self.s_vector = vector(point)
            self.ss = dict(zip(graph.ss_list, self.s_vector))
            
            self.q_vector = H.inverse() * vector([ln(x) for x in (H * self.s_vector)])
            self.q_vector = vector([x.simplify_full() for x in self.q_vector])
            
            self.qs = dict(zip(graph.qs_list, self.q_vector))

            self.ys = {}
            for (yy, qq) in graph.qs_map.items():
                self.ys[yy] = exp(-2 * self.qs[qq]).simplify_exp() 

            self.y_vector = vector([self.ys[y] for y in graph.ys_list])
       
    def __str__(self):
        s  = "\n".join(map(str, list(self.y_vector)))
            
        return s

    def subs(self, d):
        return Point(self.graph, self.y_vector.subs(d))

class Julia:
    def __init__(self):
        print "Initializing Julia."
        self.graph = my_graphs["tree4"]

        self.cs = list(var("c1 c2 c3 c4 c5 c6"))
        
        # my version of chor's solution C
        self.solC = Point(self.graph,
            [1/4*sqrt(-(c1 - c2 + c3 - c4)*(c2 + c3 + c5 + c6 - 16)/(c1 + c2 - c3 - c4)), \
             1/4*sqrt(-(c1 + c2 - c3 - c4)*(c2 + c3 + c5 + c6 - 16)/(c1 - c2 + c3 - c4)), \
             y_c, \
             1/16*(c1 + c2 + c3 + c4 - 16)/y_c, \
             4*y_c*sqrt(-((c1 - c4)^2 - (c2 - c3)^2)/(c2 + c3 + c5 + c6 - 16)) \
             /(16*y_c^2 + c1 + c2 + c3 + c4 - 16)] \
        ) 

        # my solution for 4 compatible characters    
        self.comp4 = Point(self.graph,[ \
             sqrt((c1 - c2 - c3 + c4)/(c1 + c2 + c3 + c4) * \
                  (c1 - c2 + c3 - c4)/(c1 + c2 - c3 - c4 - 16)), \
             sqrt((c1 - c2 - c3 + c4)/(c1 + c2 + c3 + c4) * \
                  (c1 + c2 - c3 - c4 - 16)/(c1 - c2 + c3 - c4)), \
             y_c, \
             (c1 + c2 + c3 + c4 - 16)/(16 * y_c), \
             sqrt((c1 + c2 + c3 + c4)/(c1 - c2 - c3 + c4)) * \
             sqrt((c1 + c2 - c3 - c4 + 16)*(c1 - c2 + c3 - c4)) * \
             y_c/(16 * y_c^2 + c1 + c2 + c3 + c4 - 16)] \
        )
        
        S = sqrt((c2 + c3 - c4 - c5)^2 - 64 * (c2 + c3 + c4 + c5) + 1024)   
        self.S = S
        #var('S')

        # five compatible characters, only separate out the leaves.
        self.comp5 =  Point(self.graph, [ \
            1/4 *sqrt(1/2 * (c2 + c3 - c4 - c5 + S) * \
                (4*c1 - 3*c2 + c3 + 3*c4 + 3*c5 - 32 + S)/(4*c1 + c2 - 3*c3 + 3*c4 + 3*c5 - 32 + S)), \
            1/4 *sqrt(1/2 * (c2 + c3 - c4 - c5 + S) * \
                (4*c1 + c2 - 3*c3 + 3*c4 + 3*c5 - 32 + S)/(4*c1 - 3*c2 + c3 + 3*c4 + 3*c5 - 32 + S)), \
            1/4 *sqrt((1/2 * 3*(c2 + c3 - c4 - c5) + S) * \
                (4*c1 - c2 - c3 + 3*c4  - c5 - 8 + S)/(4*c1 - c2 - 3*c3  - c4 + 3*c5 - 8 + S)), \
            1/4 *sqrt((1/2 * 3*(c2 + c3 - c4 - c5) + S) * \
                (4*c1 - c2 - 3*c3  - c4 + 3*c5 - 8 + S)/(4*c1 - c2 - c3 + 3*c4  - c5 - 8 + S)), \
            sqrt(2)/4 * (4*c1 - c2 - c3 + 3*c4  - c5 - 8 + S) * \
            sqrt((4*c1 - 3*c2 + c3 + 3*c4 + 3*c5 - 32 + S)*(4*c1 + c2 - 3*c3 + 3*c4 + 3*c5 - 32 + S)) / \
            (sqrt(c2 + c3 - c4 - c5 + S) * (4*c1 - c2 - c3 + 3*c4 + 3*c5 - 32 + S)) \
        ])
        
    # this must have been wrong    
    def comp6(self):

        return Point(self.graph,
            [ 1/4 * sqrt((c1 + c4 + c5 + c6 - 16) * (c1 - c2 + c3 - c4)/(c1 + c2 - c3 -c4)), \
              1/4 * sqrt((c1 + c4 + c5 + c6 - 16) * (c1 + c2 - c3 - c4)/(c1 - c2 + c3 -c4)), \
              y_c, \
              (c1 + c2 + c3 + c6 -16)/(16 * y_c), \
              sqrt((c1 - c2 + c3 - c6)*(c1 + c2 - c3 -c6)/(c1 + c4 + c5 + c6 - 16)) \
                  * (4 * y_c) / (16 * y_c^2 + c1 + c2 + c3 + c6 -16)]
        )

class Chor:
    
    def __init__(self):
        print "Initializing Chor." 
        self.graph = my_graphs["tree4"]

        self.dataA = [ 7, 0, 0, 1, 0, 1, 1, 0] 
        self.dataB = [14, 0, 0, 3, 0, 2, 1, 0]
        self.dataC = [10, 2, 2, 4, 0, 1, 1, 0]

        self.sol1 = Point(self.graph, [28/50, 4/50, 4/50, 4/50, (4+z)/50, 1/50, 1/50, (4-z)/50])
        self.sol2 = Point(self.graph, [28/50, (4+z)/50, (4-z)/50, 4/50, 4/50, 1/50, 1/50, 4/50])
        
        self.solB = Point(self.graph, [476/800, (51-z)/800, (51+z)/800, 102/800, \
                                       51/800*(1-11/z), 12/800, 6/800, 51/800*(1 + 11/z)])

        self.solC = Point(self.graph, [90/200, 27/200, 27/200, 36/200, z/200, 3/200, 3/200, (14-z)/200])

        self.chars = \
        [[(a, 0), (b, 0), (c, 0), (d, 0)], \
         [(a, 1), (b, 1), (c, 0), (d, 0)], \
         [(a, 1), (b, 0), (c, 1), (d, 0)], \
         [(a, 0), (b, 1), (c, 1), (d, 0)]]

        self.charsC = \
        [[(a, 0), (b, 0), (c, 0), (d, 0)], \
         [(a, 1), (b, 0), (c, 0), (d, 0)], \
         [(a, 0), (b, 1), (c, 0), (d, 0)], \
         [(a, 1), (b, 1), (c, 0), (d, 0)], \
         [(a, 1), (b, 0), (c, 1), (d, 0)], \
         [(a, 0), (b, 1), (c, 1), (d, 0)]]

        self.cs = list(var("c1 c2 c3 c4 c5 c6"))
                  
    def likelihood(self, data, sol):
    
        lnL = 0
        
        L = 1 

        for i in range(0, len(data)):
            lnL = lnL + data[i] * ln(self.graph.ss_list[i])        
            L = L * self.graph.ss_list[i] ^ data[i]

        H = hadamard_matrix(2^(len(self.graph.leaves) - 1))
        
        q = vector(self.graph.qs_list)
        ss_in_q_terms = H.inverse() * vector([exp(x) for x in (H * q)])
       
        qsub = dict(zip(self.graph.ss_list, ss_in_q_terms))

        print lnL.subs(qsub)
        print lnL.subs(qsub).simplify_exp()
        print
        print lnL.subs(qsub).subs(sol.qs)
        print
        print L.subs(qsub)
        print
        print L.subs(qsub).subs(sol.qs)

#chor = Chor()
#julia = Julia()
