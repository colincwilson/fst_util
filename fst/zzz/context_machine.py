#!/usr/bin/env python
# -*- coding: utf-8 -*-

import string

# construct acceptor (identity transducer) in which each state codes for the immediately preceding and following context (currently, only length-1 contexts)
def LR_context_acceptor(Sigma, word_begin, word_end, epsilon):
    global Q, T, stateid
    q_initial = ('-',word_begin)
    Q = {q_initial : 0} # state set
    T = []  # transition set
    stateid = 1
    Sigma.append(word_end)
    def build(Q_old):
        global Q, T, stateid
        Q_new = {}
        for (x,y) in Q_old:
            if y==word_end: continue
            q_src_id = Q[(x,y)]
            for z in Sigma:
                q_dest = (y,z)
                if q_dest in Q:
                    q_dest_id = Q[q_dest]
                else:
                    Q_new[q_dest] = stateid; stateid += 1
                    q_dest_id = Q_new[q_dest]
                T.append([q_src_id, y, q_dest_id])
        if not len(Q_new)==0:
            Q.update(Q_new)
            build(Q_new)
        return (Q, T)
    (Q,T) = build(Q)

    # add final state and transitions
    q_final = (word_end, '-')
    q_final_id = stateid
    Q[q_final] = stateid
    for q_src in [(x,y) for (x,y) in Q if y==word_end]:
        q_src_id = Q[q_src]
        T.append([q_src_id, word_end, q_final_id])

    # add explicit output epsilon transitions
    for (x,y) in Q:
        if y==word_begin or x==word_end: continue
        q = Q[(x,y)]
        T.append([q, epsilon, q])

    return (Q, T, q_initial, q_final)

# Test
Sigma = ['a', 'b', 'c', 'd']
#Sigma = list(string.ascii_lowercase)
(Q, T, q_initial, q_final) = LR_context_acceptor(Sigma, '<#', '#>', 'ε')

# print in dot format
print('digraph G {')
for q in Q:
    style = 'shape=circle,'
    if q == q_initial: style = 'shape=circle, style=\"bold\",'
    if q == q_final: style = ' shape=doublecircle,'
    print(Q[q], '['+style+'label=\"('+q[0]+','+q[1]+')''\"]')


for q in Q:
    q_id = Q[q]
    for t in [t for t in T if t[0]==q_id]:
        print(t[0], '->', t[2], '[label=\"'+t[1]+'\"]')
print('}')