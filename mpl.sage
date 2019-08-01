from graph import MyGraph, labelings, labeled_likelihood, average_likelihood
from sage.all import *

g = MyGraph("tree4", [('1','u'), ('2','u'), ('u','v'), ('v','3'), ('v','4')])
#g = MyGraph("tree5", [('1','x'), ('2','x'), ('x','y'), ('5','y'), ('y','z'), ('z','3'), ('z','4')])
#g = MyGraph("tree4", [('1','u'), ('2','u'),('u','3')])

chars, parsimony_labels = g.compatible_chars()
# this is hard-coded for printing reasons
# char_labels = ["0","1","2","3","4", "uv"]

m = len(chars)

for x in range(m):
    print chars[x], parsimony_labels[x]

L, I = g.sort_leaves()
labels = labelings(len(I))
num_labels = len(labels)

NUM_VALUES = 20

def label_per_char(copies):
    # the function is going to be encoded by vector of length m
    # where the ith entry tells us the index of which label we use
    ll = [[-1]] * m
    for i in range(m):
        if copies[i] != 0: 
            ll[i] = range(len(labels))
   
    #funcs = CartesianProduct(*ll) 
    funcs = cartesian_product(ll) 
    
    print len(funcs.list()), "functions"

    exps = set()
    cache = set()
    
    for f in funcs.list(): 
        ps, ns = get_exps(f, copies)
        ps, ns = tuple(ps), tuple(ns)

        print ps, ns,
        
        key = ",".join(map(str, ps + ns))
        
        if key in cache:
            print "skipped."
            continue
        else:
            print
            cache.add(key)

        exps.add((ps,ns))
    
    return exps
        
def sort_functions(exps):

    interior_max, border_max = [],[]

    cache = set()
    for (ps, ns) in exps: 
         
        pt = get_point(ps, ns) 
        key = ",".join(map(str, pt))
        
        if key in cache:
            continue
        else:
            cache.add(key)

        # throw out point in negative orthant
        if len(filter(lambda x: 0 > round(x,3), pt)) > 0: continue

        max_val = eval(ps, ns, pt)
 
        if 1 in pt or 0 in pt:
            border_max.append((max_val, pt, ps, ns)) 
        else:
            interior_max.append((max_val, pt, ps, ns)) 
    
    interior_max = sorted(interior_max, reverse=True)
    border_max = sorted(border_max, reverse=True)

    return interior_max, border_max

def find_lumps(interior_max, border_max): 

    im, bm = len(interior_max), len(border_max)
    point_list = interior_max + border_max

    print im, "interior maxima,", bm, "border maxima."
  
    print "border maxima"
    for x in border_max:
        print x
    print

    #for x in point_list: print x

    interior_nonlumps, border_nonlumps = {}, {}
    
    min_parsimony_score, min_parsimony_index = m * len(g.ys_list), 0

    for i in range(im):
        (max1, pt1, ps1, ns1) = point_list[i]

        for j in range(i+1, im):
            
            # j has already been covered
            if j in interior_nonlumps: continue

            (max2, pt2, ps2, ns2) = point_list[j]
            
            # function 1 covers max2
            if eval(ps1, ns1, pt2) > max2: interior_nonlumps[j]=i

    interior_lumps = set(range(im)) - set(interior_nonlumps.keys())

    print "\nCurrent interior lumps:"
    for x in interior_lumps:
        print point_list[x]
    print

    for i in range(im, im+bm):

        (max1, pt1, ps1, ns1) = point_list[i]
      
        for j in range(i+1, im+bm):

            # j has already been covered
            if j in border_nonlumps: continue

            (max2, pt2, ps2, ns2) = point_list[j]

            # function 1 does cover max2
            if eval(ps1, ns1, pt2) > max2: border_nonlumps[j] = i

    border_lumps = set(range(im, im+bm)) - set(border_nonlumps.keys())
    #print interior_lumps, border_lumps

    print "\nCurrent border lumps:"
    for x in border_lumps:
        print point_list[x]
    print

    print "Interior lumps: %s, border lumps: %s." % (len(interior_lumps), len(border_lumps))
    print "\nCross-checking interior against border maxima ..."

    for i in interior_lumps:
        (max1, pt1, ps1, ns1) = point_list[i]
        
        for j in border_lumps:
            (max2, pt2, ps2, ns2) = point_list[j]
      
            if i in interior_nonlumps: continue
            
            if j in border_nonlumps: continue

            if max1 > max2:
                if eval(ps1, ns1, pt2) > max2: 
                    border_nonlumps[j]=i
                    print "   Border max %s at %s covered by func with max %s" \
                          % (max2, pt2, max1)
            
            elif max2 > max1:
                if eval(ps2, ns2, pt1) > max1: 
                    interior_nonlumps[i]=j
                    print "   Interior max %s at %s covered by func with max %s" \
                          % (max1, pt1, max2)
    print "Done.\n"

    interior_lumps = interior_lumps - set(interior_nonlumps.keys())
    border_lumps = border_lumps - set(border_nonlumps.keys())
   
    print "Interior lumps: %s, border lumps: %s.\n" % (len(interior_lumps), len(border_lumps))
    #print interior_lumps, border_lumps

    print "FINAL INTERIOR LUMPS:", len(interior_lumps)
    interior_pts = sorted([point_list[i] for i in interior_lumps], reverse=True) 
    for x in interior_pts: print x

    print
    print "FINAL BORDER LUMPS:", len(border_lumps)
    border_pts = sorted([point_list[j] for j in border_lumps], reverse=True) 
    for x in border_pts: print x

    return interior_pts, border_pts 

def craziness():
  
    NUM = 7
    CHAR = 5

    print "Suppose we have %s copies of the %s character, \nbut each one \
can have different labels." % (NUM, chars[CHAR])
    
    ll = [[-1]] * NUM
    for i in range(NUM):
        ll[i] = range(len(labels))
   
    #funcs = CartesianProduct(*ll) 
    funcs = cartesian_product(*ll) 
    print len(funcs.list()), "functions"

    exps = set()
    cache = set()
    
    for f in funcs.list(): 

        e = len(g.ys_list)
        ps, ns  = [0] * e, [0] *e
        
        for x in f: 
            l, copies = [-1] * m, [0] * m
            l[CHAR] = x
            copies[CHAR] = 1 

            p_frag, n_frag = get_exps(l, copies) 
            
            for x in range(e):
                ps[x] += p_frag[x]
                ns[x] += n_frag[x]

        ps, ns = tuple(ps), tuple(ns)

        print ps, ns,
        
        key = ",".join(map(str, ps + ns))
        
        if key in cache:
            print "skipped."
            continue
        else:
            print
            cache.add(key)

        exps.add((ps,ns))
    
    ip, bp = sort_functions(exps)
    im, bm = find_lumps(ip, bp)


def summarize_lumps(data):
    exps = label_per_char(data)
    ip, bp = sort_functions(exps)
    il, bl = find_lumps(ip, bp)
    
    return il, bl
    #print_lumps(il, bl)

def print_lumps(interior_lumps, border_lumps):
    """Prints two latex tables for lumps."""
    
    ip, bp = len(interior_lumps), len(border_lumps)
    pts = interior_lumps + border_lumps

    print "\n\\begin{tabular}{"+ "c" * 5 +"}"
    print " interior log max v & $\\vec{y^*}$ & $\\vec{p}$ & $\\vec{n}$ & parsimony score \\\\"

    for i in range(len(pts)): 
        max_val, pt, ps, ns, = pts[i]
        print max_val, 
        print "& $\Big(%s\Big)$ " % ",".join(map(latex, pt))
        print "& (%s)" % ",".join(map(str, ps)),
        print "& (%s)" % ",".join(map(str, ns)),
        print "& %s" % sum(ns),
        print "\\\\"

        if i == ip-1:
            print "\\end{tabular}"
            print "\n\\begin{tabular}{"+ "c" * 5 +"}"
            print " border log max & $\\vec{y^*}$ & $\\vec{p}$ & $\\vec{n}$ & parsimony score \\\\"
    
    print "\\end{tabular}"


def eval(ps, ns, pt):
    """Evaluates the likelihood function for given p/n vectors at point pt."""
    c = ps[0] + ns[0]
    e = len(g.ys_list)

    #print ps, ns, pt

    #TODO: 
    # need to look up the sage documentation
    # which is better for precision?? don't know

    logval = 0
    for i in range(e):
        logval += ps[i] * log(1 + pt[i]) 
        
        if ns[i] != 0:
            logval += ns[i] * log(1 - pt[i])

    return round(logval, 4)

def get_point(ps, ns):
    """Returns the max point for the function with given exponent vectors. Note that this
    isn't the critical point; it might be on the boundary."""        
    pt = []

    for i in range(len(ps)):
        if ps[i] == 0:
            pt.append(0)
        elif ns[i] == 0:
            pt.append(1)
        else:
            x,y = ps[i], ns[i]
            pt.append((x - y)/(x + y))

    return pt

def walk(copies):

    L = map(unzip, parsimony_labels)
    L = map(lambda x: labels.index(x[1]), L)

    max_ps, max_ns = map(tuple, get_exps(L, copies))
    max_pt = get_point(max_ps, max_ns)

    orig_entry = []
    for j in range(m):
        X = [0]* num_labels
        X[L[j]] = copies[j]
        orig_entry += [X]   

    print orig_entry

    old_list = [orig_entry]
    new_list = []
    pn_cache = set()
    hits = 0
  
    from sys import stdout
    
    lumps = []
    covered_lumps = set()
    #lumps = [(max_ps, max_ns)]
    steps = 0
    
    while len(old_list) > 0:
        
        print "\nSteps: %s, Lumps: %s." % (steps, len(lumps) - len(covered_lumps))
        
        for z in range(len(old_list)):
            
            entry = old_list[z]

            ps, ns = get_exps_from_entry(copies, entry)
            pt  = get_point(ps, ns)
            val = eval(ps, ns, pt) 

            # ya can't be a lump if you're not in the right quadrant
            if not any(filter(lambda x: x < 0, pt)): 
                is_lump = True
                covered_someone = False
               
                for l in range(len(lumps)):
                    
                    if not l in covered_lumps:
                    
                        lump_ps, lump_ns = lumps[l]
                        
                        if is_lump and val < eval(lump_ps, lump_ns, pt):
                            is_lump = False
    
                        lump_pt = get_point(lump_ps, lump_ns)
                        
                        if eval(ps, ns, lump_pt) > eval(lump_ps, lump_ns, lump_pt):
                            covered_lumps.add(l)
                            covered_someone = True

                if is_lump:
                    if covered_someone: 
                        print '+', pt, val
                    else:
                        print '*', pt, val
                    lumps.append((ps, ns))
                
                elif covered_someone:
                    print '-', pt
                else:
                    print ' ', pt


            for j in range(1, m):

                if sum(entry[j]) == 0: continue

                for x in toggle(entry[j], opposite_corner(orig_entry[j])):

                    new_entry = deepcopy(entry)
                    new_entry[j] = x
                    
                    new_ps, new_ns = get_exps_from_entry(copies, new_entry)
                    
                    if (new_ps, new_ns) in pn_cache: 
                        hits +=1
                        continue
                    
                    #new_pt = get_point(new_ps, new_ns)
                    #if any(filter(lambda x: x < 0, new_pt)): continue

                    pn_cache.add((new_ps, new_ns))
                    new_list.append(new_entry)
        steps+=1

        old_list = new_list
        new_list = []
    
    print "Sorting."
    sorted_lumps = []
    for l in range(len(lumps)):
        
        if l in covered_lumps: continue
        
        ps, ns = lumps[l]
        pt = get_point(ps, ns)
        val = eval(ps,ns,pt)
        sorted_lumps.append((val, pt))
    
    sorted_lumps = sorted(sorted_lumps, reverse=True)
    
    print "\nLUMPS:"
    for (v, pt) in sorted_lumps: print v, pt

    print "\np,n cache is size %s, hits: %s" % (len(pn_cache), hits)

    return sorted_lumps

def neg(i):
    L = []
    for bit in labels[i]:
        L.append(1 - bit)
    
    return labels.index(L)

def neighbors(i):
    neighbors = []
    for x in range(len(labels[i])):
        LL = copy(labels[i])
        LL[x] = 1 - LL[x]
        neighbors.append(labels.index(LL))
    return neighbors

def distance(i,j):
    d = 0
    for x in range(len(labels[i])):
        d += labels[i][x] != labels[j][x]
    
    return d
  
def opposite_corner(L):
    
    assert(len(L) == num_labels)
    
    M = [0] * num_labels

    for x in range(num_labels): M[x] = L[neg(x)]
    
    return M

def toggle(changes, goal):
    """Helper for walk function that generates all 1-label flips of changes."""
    L = []

    nonzeros = []
    for i in range(num_labels):
        if goal[i] != 0: nonzeros.append(i)

    #print nonzeros

    for i in range(num_labels):
       
        if changes[i] == 0: continue

        for j in neighbors(i):

            backwards=False
            for n in nonzeros:
                if distance(n, i) < distance(n, j):
                    backwards=True
                    break

            if backwards: continue

            new_change = copy(changes)
            new_change[i] -= 1
            new_change[j] += 1
            
            L.append(new_change)
    
    return L     


def get_exps_from_entry(copies, entry):
    """This is a helper for the walk function."""
    
    y = len(g.ys_list)
    ps, ns = get_char_exps(make_char(0,0), copies[0])

    for j in range(1, m): #constant char gets no label changes
        
        changes = entry[j]

        assert(sum(changes) == copies[j])
        
        for label_index in range(num_labels):
            
            char = make_char(j, label_index)
            
            cps, cns = get_char_exps(char, copies=changes[label_index])
            
            for i in range(y):
                ps[i] += cps[i]
                ns[i] += cns[i]

    return tuple(ps), tuple(ns)

def make_char(char_index, label_index):

    return dict(chars[char_index] + zip(I, labels[label_index]))

def get_char_exps(char, copies=1):
    
    l = len(g.ys_list)
    ps = [0] * l
    ns = [0] * l
    
    for ((a,b), y) in g.ys_map.items():
        
        i = g.ys_list.index(y)

        if(char[a] == char[b]):
            ps[i] += copies
        else:
            ns[i] += copies
    
    return ps, ns


def avg_likelihood_info(data, pt):
    """This will generate the average likelihood function for your data and
    evaluate it at pt, and return the log value of the likelihood and the gradient
    at pt."""

    ML = average_likelihood(g, data, normalizing=True)
    
    pt_dict = dict(zip(g.ys_list, pt))
    
    ml = round(log(ML).subs(pt_dict),2)
    
    grad = [round(log(ML).diff(y).subs(pt_dict),2) for y in g.ys_list]

    return ml, grad


def get_exps(func, copies):
    """Converts the function label encoding and data set into two exponent vectors,
    one corresponding to the likelihood that there is a change on the edge and the
    other corresponding to no change."""

    l = len(g.ys_list)
    ps, ns = [0] * l, [0] * l
    
    for j in range(m):
       
        if copies[j] == 0: continue

        #char = dict(chars[j] + zip(I, labels[func[j]]))
        #print char, "x", copies[j]

        ps1, ns1 = get_char_exps(make_char(j, func[j]), copies=copies[j])
        for i in range(l):
            ps[i] += ps1[i]
            ns[i] += ns1[i]

    return ps, ns

def unzip(L):
    L1, L2 = [], []
    for (x, y) in L:
        L1.append(x)
        L2.append(y)
    return L1, L2

random_data = [] 
for i in range(len(chars)):
    random_data.append(randint(2,7))

data = [3,1,1,1,1,2]

print "The usual example we use is in the variable data =", data
print "Here's some random data: random_data =", random_data

print """
Usage:
    summarize_lumps(data)
    int_max, border_max = sort_functions(data); find_lumps(int_max, border_max)
    
The sort_functions procedure takes a long time, so you can get intermediate values.

    print_lumps(find_lumps(int_max, border_max))  # prints LaTeX tables. 
"""

print "You can make your own data, it just needs to be a list of length %s." % m


