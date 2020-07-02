# -*- coding: utf-8 -*-

import itertools, string
from collections import namedtuple
import pywrapfst

FST = namedtuple('FST', ['Q', 'T', 'q0', 'qf'])
Transition = namedtuple('Transition', ['src', 'label', 'dest'])


def context_acceptor(Sigma, length=1, begin_delim='>', end_delim='<', epsilon='ε'):
    """
    Construct acceptor (identity transducer) for segments 
    in immediately preceding and following contexts
    xxx '⋊', '⋉'
    """
    Sigma_plus = Sigma + [begin_delim, end_delim, epsilon]
    # States are segment histories
    Q = { q for q in itertools.permutations(Sigma_plus, length) }
    # Transitions are defined by one-segment history updates
    T = { Transition(q1, x, q1[1:]+(x,)) for q1 in Q for x in Sigma }
    # Initial state and outgoing transitions
    q0 = (epsilon,)*length
    T |= {Transition(q0, begin_delim, q0[1:]+(begin_delim,))}
    # Final states and incoming transitions
    qf = set(filter(lambda q: q[-1] == end_delim, Q))
    T |= {Transition(q, end_delim, q[1:]+(end_delim,)) for q in Q}
    return (FST(Q, T, q0, qf))


def trim(M):
    """
    Remove dead states and transitions from FST M
    """
    # Forward pass
    Qforward = set([M.q0])
    Qnew = set([M.q0])
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
    T = set([t for t in M.T \
            if (t.src in Q) and (t.dest in Q)])

    # Backward pass
    Qbackward = set([q for q in M.qf])
    Qnew = set([q for q in M.qf])
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
    T = set([t for t in T \
             if t.src in Q and t.dest in Q])

    q0 = M.q0 if M.q0 in Q else None
    qf = {q for q in M.qf if q in Q}
    M_trim = FST(Q, T, q0, qf)
    return M_trim


# Test
Sigma = ['a']
#Sigma = list(string.ascii_lowercase)
A = context_acceptor(Sigma, 2)
print(A.Q)
print(A.q0)
print(A.qf)
for t in A.T:
    print(t)
print()
A_trim = trim(A)
print(A_trim.Q)
print(A_trim.q0)
print(A_trim.qf)
for t in A_trim.T:
    print(t)
#print(A_trim)
