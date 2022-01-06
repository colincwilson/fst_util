import sys
from pathlib import Path

sys.path.append(str(Path.home() / 'Code/Python/fst_util'))
from fst_util.fst_config import *
from fst_util import proc
from fst_util import ostia
from fst_util.fst import *


def test():
    D = [('a', '1'), ('b', '1'), ('a a', '0 1'), ('a b', '0 1'),
         ('a a a', '0 0 1'), ('a b a b', '0 1 0 1')]
    #D = [('a a', '1')]
    fst1, Sigma, Lambda = proc.prefix_tree(D)
    print(fst1.print(acceptor=False, show_weight_one=True, missing_sym='λ'))
    fst1.draw('prefix_tree.dot')

    fst2, _, _ = proc.onward_tree(fst1, fst1.start(), λ)
    print(fst2.print(acceptor=False, show_weight_one=True, missing_sym='λ'))
    fst2.draw('onward_tree.dot')

    fst3 = ostia.ostia(D)


if __name__ == '__main__':
    test()