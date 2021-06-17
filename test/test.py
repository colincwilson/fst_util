import sys
from pathlib import Path

sys.path.append(str(Path.home() / 'Code/Python/fst_util'))
from fst_util import fst_config
from fst_util.fst import *


def test():
    # State labels
    fst_config.init(sigma_syms=['a', 'b'], special_syms=['λ'])
    fst = Fst(fst_config.symtable)
    for sym in ['λ', '⋊', 'A', 'B']:
        fst.add_state(sym)
    fst.set_start('λ')
    for sym in ['A', 'B']:
        fst.set_final(sym)
    fst.add_arc(src='λ', ilabel='⋊', dest='⋊')
    fst.add_arc(src='⋊', ilabel='a', dest='A')
    fst.add_arc(src='⋊', ilabel='b', dest='B')
    fst.add_arc(src='A', ilabel='a', olabel='a', dest='A')
    fst.add_arc(src='A', ilabel='a', olabel='b', dest='B')
    #print(fst.print(acceptor=True, show_weight_one=True))
    fst.draw('tmp.dot')

    # Left- and right- context acceptors
    fst_config.init(sigma_syms=['a', 'b'], special_syms=['λ'])
    L = left_context_acceptor(context_length=2)
    L.draw('L.dot')
    R = right_context_acceptor(context_length=2)
    R.draw('R.dot')

    # Accepted strings
    print(accepted_strings(L, 'input', 4))

    # Connect with state labels preserved
    C = Fst(fst_config.symtable)
    qf = C.add_state('0')
    q = C.add_state('1')
    q0 = C.add_state('2')
    C.set_start(q0)
    C.set_final(qf)
    #C.add_arc(src=q0, ilabel='a', dest=q)
    C.add_arc(src=q0, ilabel='a', dest=qf)
    C.add_arc(src=q, ilabel='b', dest=qf)
    print(C._state2label)
    C_trim = C.connect()
    print(C_trim._state2label)
    C_trim.draw('C_trim.dot')

    # Composition
    fst_config.init(sigma_syms=['a', 'b'])
    M1 = Fst(fst_config.symtable)  # a*b*
    for q in [0, 1]:
        M1.add_state(q)
    M1.set_start(0)
    M1.set_final(1)
    M1.add_arc(src=0, ilabel='a', dest=0)
    M1.add_arc(src=0, ilabel='b', dest=1)
    M1.add_arc(src=1, ilabel='b', dest=1)

    M2 = Fst(fst_config.symtable)  # ab*
    for q in [0, 1]:
        M2.add_state(q)
    M2.set_start(0)
    M2.set_final(1)
    M2.add_arc(src=0, ilabel='a', dest=1)
    M2.add_arc(src=1, ilabel='b', dest=1)
    M = compose(M1, M2)
    M.draw('M.dot')

    # Arc deletion
    fst_config.init(sigma_syms=['a', 'b'])
    fst = Fst(fst_config.symtable)
    for q in [0, 1]:
        fst.add_state(q)
    fst.set_start(0)
    fst.set_final(1)
    fst.add_arc(src=0, ilabel='a', dest=1)
    fst.add_arc(src=0, ilabel='b', dest=1)
    print(fst.print())


if __name__ == '__main__':
    test()