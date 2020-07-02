# -*- coding: utf-8 -*-

import itertools, re, string, sys
from collections import namedtuple
from . import fst_config as config
verbosity = 0

FST = namedtuple('FST', ['Q', 'T', 'q0', 'qf'])

class Transition():
    def __init__(self,
        src=None, ilabel=None, olabel=None, weight=None, dest=None):
        self.src = src
        self.ilabel = ilabel if ilabel is not None else olabel
        self.olabel = olabel if olabel is not None else ilabel
        self.weight = weight
        self.dest = dest
    
    def __str__(self):
        return '('+ str(self.src) +','+ str(self.ilabel) \
                +','+ str(self.olabel) +','+ str(self.dest) +')'
    
    def __repr__(self):
        return self.__str__()

#Transition = namedtuple('Transition', ['src', 'label', 'dest'])

def suffix(alpha, l):
    """
    Extract suffix of length (l-1)
    """
    if len(alpha) < l:
        return alpha
    else:
        return alpha[1:]


def prefix(alpha, l):
    """
    Extract prefix of length (l-1)
    """
    if len(alpha) < l:
        return alpha
    else:
        return alpha[:-1]


def left_context_acceptor(Sigma, length=1):
    """
    Construct acceptor (identity transducer) for segments 
    in immediately preceding contexts
    """
    begin_delim, end_delim = \
        config.begin_delim, config.end_delim

    # Initial state and outgoing transition
    q0, q1 = ('λ',), (begin_delim,)
    Q = {q0, q1}
    T = { Transition(src=q0, olabel=begin_delim, dest=q1) }

    # Interior transitions
    # xα -- y --> αy for each y
    Qnew = set(Q)
    for l in range(length+1):
        Qold = set(Qnew)
        Qnew = set()
        for q1 in Qold:
            if q1 == q0: continue
            for x in Sigma:
                q2 = suffix(q1,length) + (x,)
                T.add(Transition(src=q1, olabel=x, dest=q2))
                Qnew.add(q2)
        Q |= Qnew
    
    # Final states and incoming transitions
    qf = set()
    for q1 in Q:
        if q1 == q0: continue
        q2 = suffix(q1,length) + (end_delim,)
        T.add(Transition(src=q1, olabel=end_delim, dest=q2))
        qf.add(q2)
    Q |= qf

    A = FST(Q, T, q0, qf)
    A_trim = trim(A)
    return A_trim


def right_context_acceptor(Sigma, length=1):
    """
    Construct acceptor (identity transducer) for segments 
    in immediately following contexts
    """
    begin_delim, end_delim = \
        config.begin_delim, config.end_delim

    # Final state and incoming transition
    qf, qp = ('λ2',), (end_delim,)
    Q = {qf, qp}
    T = { Transition(src=qp, olabel=end_delim, dest=qf) }

    # Interior transitions
    # xα -- x --> αy for each y
    Qnew = set(Q)
    for l in range(length+1):
        Qold = set(Qnew)
        Qnew = set()
        for q2 in Qold:
            if q2 == qf: continue
            for x in Sigma:
                q1 = (x,) + prefix(q2, length)
                T.add(Transition(src=q1, olabel=x, dest=q2))
                Qnew.add(q1)
        Q |= Qnew

    # Initial state and outgoing transitions
    q0 = ('λ',)
    for q in Q:
        if q == qf: continue
        T.add(Transition(src=q0, olabel=begin_delim, dest=q))
    Q.add(q0)

    A = FST(Q, T, q0, {qf})
    A_trim = trim(A)
    return A_trim


def intersect(M1, M2):
    """
    Intersect two FSTs
    todo: speed up with label/state indexing; 
    generalize with matcher
    """
    Q1, T1, q0_1, qf_1 = \
        M1.Q, M1.T, M1.q0, M1.qf
    Q2, T2, q0_2, qf_2 = \
        M2.Q, M2.T, M2.q0, M2.qf

    Q = set()
    T = set()
    q0 = (q0_1, q0_2)
    qf = {(q1, q2) for q1 in qf_1 for q2 in qf_2}
    Q.add(q0)
    Q |= qf

    Qnew = {q0}
    while len(Qnew) != 0:
        Qold = set(Qnew)
        Qnew.clear()
        #print(Qold)
        for q in Qold:
            q1, q2 = q
            if verbosity>0: print(q, '-->', q1, 'and', q2)
            for t1 in filter(lambda t: t.src == q1, T1):
                x = t1.olabel
                for t2 in filter(lambda t: t.src == q2 and t.olabel == x, T2):
                    r = (t1.dest, t2.dest)
                    T.add(Transition(src=q, olabel=x, dest=r))
                    if r not in Q:
                        Qnew.add(r)
        Q.update(Qnew)
    #print(len(Q), len(T))
    M = trim(FST(Q, T, q0, qf))
    return M


def trim(M):
    """
    Remove dead states and transitions from FST M
    """
    # Forward pass
    Qforward = {M.q0}
    Qnew = {M.q0}
    while len(Qnew) != 0:
        Qold = set(Qnew)
        Qnew.clear()
        for q in Qold:
            q_T = [t for t in M.T if t.src == q]
            for t in q_T:
                r = t.dest
                if r not in Qforward:
                    Qforward.add(r)
                    Qnew.add(r)

    Q = Qforward.copy()
    T = {t for t in M.T \
            if (t.src in Q) and (t.dest in Q)}

    # Backward pass
    Qbackward = {q for q in M.qf}
    Qnew = {q for q in M.qf}
    while len(Qnew) != 0:
        Qold = set(Qnew)
        Qnew.clear()
        for r in Qold:
            r_T = [t for t in T if t.dest == r]
            for t in r_T:
                q = t.src
                if q not in Qbackward:
                    Qbackward.add(q)
                    Qnew.add(q)

    Q &= Qbackward
    T = {t for t in T \
             if t.src in Q and t.dest in Q}

    q0 = M.q0 if M.q0 in Q else None
    qf = {q for q in M.qf if q in Q}
    M_trim = FST(Q, T, q0, qf)
    return M_trim


def reverse(M):
    """
    Reverse FST, retaining delimiter semantics
    """
    begin_delim, end_delim, epsilon = \
        config.begin_delim, config.end_delim, \
        config.epsilon

    # Reverse
    q0 = ('λ2',)
    qf = {M.q0}
    Q = {q0, qf}

    T = { Transition(src=t.dest, olabel=t.olabel, dest=t.src) for t in M.T }
    for q1 in M.qf:
        T.add( Transition(src=q0, olabel=epsilon, dest=q1) )
    
    return FST(Q, T, q0, qf)


def flatten(M):
    """
    Flatten states of FST (e.g., after intersection)
    see http://stackoverflow.com/questions/3204245/how-do-i-convert-a-tuple-of-tuples-to-a-one-dimensional-list-using-list-comprehe
    xxx define flatten_state() !
    """
    T = { Transition(src=flatten_state(t.src), olabel=t.olabel, dest=flatten_state(t.dest)) for t in T }
    q0 = flatten_state(q0)

    return sum(q[0:-1], ()) + (q[-1],)


def linear_acceptor(x):
    """
    Linear acceptor for space-delimited string
    """
    Q = {0}
    T = set()
    x = x.split(' ')
    for i in range(len(x)):
        Q.add(i+1)
        T.add(Transition(src=i, olabel=x[i], dest=i+1))
    M = FST(Q, T, 0, {len(x)})
    return M


def trellis(max_len):
    """
    Trellis for strings of length 0 to max_len
    (begin/end delimiters not included in lengths)
    """
    begin_delim = config.begin_delim
    end_delim = config.end_delim
    Sigma = config.Sigma

    Q, T = set(), set()
    q0 = 0; Q.add(q0)
    q1 = 1; Q.add(q1)
    T.add(Transition(src=q0, olabel=begin_delim, dest=q1))

    qe = max_len+1; Q.add(qe)
    qf = max_len+2; Q.add(qf)
    T.add(Transition(src=qe, olabel=end_delim, dest=qf))

    for i in range(max_len):
        q = i + 1
        r = i + 2
        for x in Sigma:
            Q.add(r)
            T.add(Transition(src=q, olabel=x, dest=r))
        T.add(Transition(src=q, olabel=end_delim, dest=qf))
    return FST(Q, T, q0, {qf})


def map_states(M, f):
    """
    Apply function f to each state
    """
    Q = { f(q) for q in M.Q }
    T = { Transition(src=f(t.src), olabel=t.olabel, dest=f(t.dest)) \
            for t in M.T }
    q0 = f(M.q0)
    qf = {f(q) for q in M.qf}
    return FST(Q, T, q0, qf )


def accepted_strings(M, max_len):
    """
    Accepted strings up to maximum length 
    (not including begin/end delimiters)
    """
    prefixes = { (M.q0, '') }
    prefixes_new = prefixes.copy()
    for i in range(max_len+2):
        prefixes_old = set(prefixes_new)
        prefixes_new = set()
        for prefix in prefixes_old:
            for t in filter(lambda t: t.src == prefix[0], M.T):
                prefixes_new.add( (t.dest, prefix[1]+' '+t.olabel) )
        prefixes |= prefixes_new
        #print(i, prefixes_new); print()

    accepted = { prefix for (state, prefix) in prefixes \
                    if re.search(config.end_delim+'$', prefix) }
    accepted = { re.sub('^ ', '', prefix) for prefix in accepted }
    return accepted


def draw(M, fname):
    """
    Print FST in dot/graphviz format
    xxx return string
    """
    Q, T, q0, qf = M.Q, M.T, M.q0, M.qf
    stateid = {q:i for i,q in enumerate(Q)}
    f = open(fname, 'w')
    f.write('digraph G {\n')
    f.write('rankdir=LR;\n')
    f.write('node [shape=circle]\n')
    for q in stateid:
        q_str = ' ['
        if q == q0:
            q_str += 'style=bold '
        if q in qf:
            q_str += 'shape=doublecircle '
        q_str += 'label=\"' + ''.join(str(q)) +'\"'
        q_str += ']\n'
        f.write(str(stateid[q]) + q_str)
    for t in T:
        try:
            label = str(t.olabel) if t.ilabel == t.olabel \
                else str(t.ilabel) +':'+ str(t.olabel)
            f.write(str(stateid[t.src]) +' -> '+ str(stateid[t.dest]) +' [label=\"' + label + '\"]\n')
        except:
            pass
    f.write('}\n')


# xxx testing
if 0:
    Sigma = ['a', 'b']
    config.Sigma = Sigma
    A_left = left_context_acceptor(Sigma, length=1)
    draw(A_left, 'A_left.dot')
    A_right = right_context_acceptor(Sigma, length=1)
    draw(A_right, 'A_right.dot')


if 0:
    verbosity = 1
    Sigma = ['a', 'b', 'c', 'd']
    config.Sigma = Sigma
    #projections = [Sigma, ['a','b']]    # segmental projection must be first!
    projections = [Sigma, ['a','b']]
    [Q, T, q0, qf] = FST(Sigma, projections[1])