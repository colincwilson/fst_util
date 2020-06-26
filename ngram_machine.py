#! /usr/bin/python

# alphabet and n-gram order
S = ['a', 'b', 'c']
k = 3

# histories up to length k
H1 = [(x) for x in S]
H2 = [(x,y) for x in S for y in S]
H3 = [(x,y,z) for x in S for y in S for z in S]

# map histories to unique state ids
H = H1 + H2
H = {H[i]:i+2 for i in xrange(len(H))}
H1 = {x:H[x] for x in H if len(x)==1}
H2 = {x:H[x] for x in H if len(x)==2}
H2 = {x:H[x] for x in H if len(x)==2}
print H1, H2

# make transition set
T = [(0, '<#', 1)]
for x in H1:
    T.append((1,x,H1[x]))
for x in H2:
    q_src = H[x[0]]
    label = x[1]
    q_dest = H2[x]
    T.append((q_src, label, q_dest))
print T
