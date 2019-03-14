# starter code imports
import sys
import networkx as nx

# my imports
import numpy as np
import matplotlib as plt
from collections import Counter, defaultdict
import math

from utils import *

# memoization
old_m = None

def main():
    inputfilename = sys.argv[1]
    outputfilename = sys.argv[2]
    compute(inputfilename, outputfilename)

#infile should be a csv
def compute(infile, outfile):
    x, x_labels = read_csv(infile)
    write_gph(k2_search(x), x_labels, outfile)

def k2_search(x):
    no_edge_dag = nx.DiGraph()
    for i in range(x.shape[1]):
        no_edge_dag.add_node(i)

    # find best neighbor; each state respresents a new base, with transitions
    # created to disease or not based on features in infile
    best_dag = no_edge_dag
    best_score = compute_score(no_edge_dag, x)
    for i in range(x.shape[1]):
        parents = range(0,i) + range(i+1,x.shape[1])
        for p in parents:
            # try a new edge on for size
            best_dag.add_edge(p,i)
            curr_score = compute_score(best_dag, x, x_changed=[i])
            # check if cycles were created
            try:
                nx.find_cycle(best_dag, source=p, orientation='original')
                creates_cycle = True
                # print("cycle found")
            except:
                creates_cycle = False
            # found new best graph, move on to next variable in x
            if curr_score > best_score and not creates_cycle:
                best_score = curr_score
            else:
                best_dag.remove_edge(p,i)
    print(best_score)
    return best_dag

def compute_score(dag, x, x_changed=None):

    # log of Bayesian probability of states given the dataset, ie. log(P(G|D))
    score = 0.0

    m = compute_m(dag, x, x_changed) # counts

    # iterate across n variables
    for i in range(len(m)):
        m_ij = defaultdict(int) # array of m_ij0, which is the sum of m_ijk from k=1 to r
        a_ij0 = 0.0
        # iterate through parental instantiations (q) and variable instantiations (r)
        for jk in m[i].keys():
            parents = jk[0]
            value = jk[1]
            score += math.lgamma(1.0 + m[i][jk]) - math.lgamma(1.0)
            m_ij[parents] += value
            a_ij0 += 1

        # iterate through parental instantiations (q) only
        for m_ij0 in m_ij.values():
            score += math.lgamma(a_ij0) - math.lgamma(a_ij0 + m_ij0)
    return -score


# creates a list of counters m that records the number of occurences of each possible scenario
# for the genomic region size
def compute_m(dag, x, x_changed=None):
    global old_m
    if not old_m:
        x_changed = None # make sure m is computed the first time
    # check whether to use memoization
    if x_changed:
        m = old_m
        for i in x_changed:
            # clear out counters that need to be updated
            m[i] = Counter()
    else:
        m = []
        # m is a list of counters, m[i] is counter for x_i
        # counters are indexed by tuples
        # (j, k) = (x index, state of x's parents as a string, value of x)
        # ('0110', 1) stores how many times the variable has value 1 with parental values 0110
        for i in range(x.shape[1]):
            m.append(Counter())
    # iterate through data points
    for data_point in x:
        # iterate through variables
        for i in range(len(data_point)):
            # memoization: filled in from old m, save time
            if x_changed and i not in x_changed:
                continue
            parent_state_string = ''
            # capture parent states for variable, if they exist
            if i in dag:
                parents = dag.predecessors(i)
                for p in parents:
                    parent_state = data_point[p]
                    parent_state_string += str(parent_state)
            # increment counter
            variable_state = data_point[i]
            m[i][(parent_state_string, variable_state)] += 1
    # print(m)
    old_m = m # memoization
    return m

if __name__ == '__main__':
    main()
