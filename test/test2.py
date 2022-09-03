import sys
from pathlib import Path

sys.path.append(str(Path.home() / 'Code/Python/fst_util'))
from fst_util import config as fst_config
from fst_util.fst import *


def test():
    # N-gram machine with one-symbol history
    config = {'sigma': ['a', 'b', 'c']}
    fst_config.init(config)
    M = left_context_acceptor(context_length=1)
    print(M.print(acceptor=True, show_weight_one=True))

    M = M.map_weights(map_type='to_log')
    print(M.print(acceptor=True, show_weight_one=True))

    M = pynini.push(M, push_weights=True)
    print(M.print(acceptor=True, show_weight_one=True))


if __name__ == '__main__':
    test()
