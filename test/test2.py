import sys
from pathlib import Path

sys.path.append(str(Path.home() / 'Code/Python/fst_util'))
from fst_util.fst_config import *
from fst_util import proc2
from fst_util import ostia2
#from fst_util.fst import *
from fst_util.simple_fst import *


def test():
    D = [('a', '1'), ('b', '1'), ('a a', '0 1'), ('a b', '0 1'),
         ('a a a', '0 0 1'), ('a b a b', '0 1 0 1')]
    #D = [('a a', '1')]
    fst1, Sigma, Lambda = proc2.prefix_tree(D)
    print(fst1.print())
    fst1_ = fst1.pynini()
    print(fst1_.print(acceptor=False, show_weight_one=True, missing_sym='λ'))
    fst1_.draw('prefix_tree.dot')
    #sys.exit(0)

    fst2 = proc2.onward_tree(fst1, fst1.q0, λ)
    print(fst2.print())
    fst2_ = fst2.pynini()
    print(fst2_.print(acceptor=False, show_weight_one=True, missing_sym='λ'))
    fst2_.draw('onward_tree.dot')

    fst3 = ostia2.ostia(D)


if __name__ == '__main__':
    test()