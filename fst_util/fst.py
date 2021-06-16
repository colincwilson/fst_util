# -*- coding: utf-8 -*-

import sys
from pathlib import Path
import pynini

sys.path.append(str(Path('../../phon')))
from phon.str_util import prefix, suffix
import fst_config as config


class Fst(pynini.Fst):
    """
    Pynini Fst with labeled inputs/outputs/states
    - Input/output labels must be strings
    - State labels must be hashable (strings, tuples, etc.) 
    """

    def __init__(self, symtable):
        super(Fst, self).__init__()
        super(Fst, self).set_input_symbols(symtable)
        super(Fst, self).set_output_symbols(symtable)
        self.state2label = {}  # State id -> label
        self.label2state = {}  # Label -> state id

    def add_state(self, state_label=None):
        # Enforce unique labels
        if state_label is not None:
            if state_label in self.label2state:
                return self.label2state[state_label]
        # Create new state
        state = super(Fst, self).add_state()
        # Self-labeling by default
        if state_label is None:
            state_label = state
        # Register state with label
        self.state2label[state] = state_label
        self.label2state[state_label] = state
        return state

    def add_arc(self, src, ilabel, olabel=None, weight=None, dest=None):
        if not isinstance(src, int):
            src = self.label2state[src]
        if not isinstance(dest, int):
            dest = self.label2state[dest]
        if isinstance(ilabel, str):
            ilabel = self.input_symbols().find(ilabel)
        if olabel is None:
            olabel = ilabel
        elif isinstance(olabel, str):
            olabel = self.output_symbols().find(olabel)
        if weight is None:
            weight = 0
        #print(src, ilabel, olabel, weight, dest)
        arc = pynini.Arc(ilabel, olabel, weight, dest)
        return super(Fst, self).add_arc(src, arc)

    def set_start(self, state):
        if not isinstance(state, int):
            state = self.label2state[state]
        return super(Fst, self).set_start(state)

    def set_final(self, state):
        if not isinstance(state, int):
            state = self.label2state[state]
        return super(Fst, self).set_final(state)

    def mutable_arcs(self, state_label):
        pass

    def print(self, **kwargs):
        # Stringify state labels
        ssymbols = pynini.SymbolTable()
        for q, label in self.state2label.items():
            ssymbols.add_symbol(str(label), q)
        return super(Fst, self).print(
            isymbols=self.input_symbols(),
            osymbols=self.output_symbols(),
            ssymbols=ssymbols,
            **kwargs)

    def draw(self, source, acceptor=True, portrait=True, **kwargs):
        # Stringify state labels
        ssymbols = pynini.SymbolTable()
        for q, label in self.state2label.items():
            ssymbols.add_symbol(str(label), q)
        return super(Fst, self).draw(
            source,
            isymbols=self.input_symbols(),
            osymbols=self.output_symbols(),
            ssymbols=ssymbols,
            acceptor=acceptor,
            portrait=portrait,
            **kwargs)


def left_context_acceptor(context_length=1, sigma_tier=None):
    """
    Acceptor (identity transducer) for segments in immediately precedin contexts (histories) of specified length. If Sigma_tier is specified as  a subset of Sigma, only contexts over Sigma_tier are tracked (other member of Sigma are skipped, i.e., label self-loops on each interior state)
    """
    epsilon = config.epsilon
    bos = config.bos
    eos = config.eos
    if sigma_tier is None:
        sigma_tier = set(config.sigma)
        sigma_skip = set()
    else:
        sigma_skip = set(config.sigma) - sigma_tier
    fst = Fst(config.symtable)

    # Initial and peninitial states
    q0 = ('λ',)
    q1 = (epsilon,) * (context_length - 1) + (bos,)
    fst.add_state(q0)
    fst.add_state(q1)
    fst.add_arc(src=q0, ilabel=bos, dest=q1)

    # Interior arcs
    # xα -- y --> αy for each y
    Q = {q0, q1}
    Qnew = set(Q)
    for l in range(context_length + 1):
        Qold = set(Qnew)
        Qnew = set()
        for q1 in Qold:
            if q1 == q0:
                continue
            for x in sigma_tier:
                q2 = suffix(q1, context_length - 1) + (x,)
                fst.add_state(q2)
                fst.add_arc(src=q1, ilabel=x, dest=q2)
                Qnew.add(q2)
        Q |= Qnew

    # Final state and incoming arcs
    qf = (eos,)
    fst.add_state(qf)
    for q1 in Q:
        if q1 == q0:
            continue
        fst.add_arc(src=q1, ilabel=eos, dest=qf)
    Q.add(qf)

    # Self-transitions labeled by skipped symbols
    # on interior states
    # xxx move outside
    for q in Q:
        if (q == q0) or (q == qf):
            continue
        for x in sigma_skip:
            fst.add_arc(src=q, ilabel=x, dest=q)

    fst.set_start(q0)
    fst.set_final(qf)
    #fst = fst.connect() # xxx handle state relabeling!
    return fst


def right_context_acceptor(context_length=1, sigma_tier=None):
    """
    Acceptor (identity transducer) for segments in immediately following contexts (futures) of specified length. If Sigma_tier is specified as a subset of Sigma, only contexts over Sigma_tier are tracked (other members of Sigma are skipped, i.e., label self-loops on each interior state)
    """
    epsilon = config.epsilon
    bos = config.bos
    eos = config.eos
    if sigma_tier is None:
        sigma_tier = set(config.sigma)
        sigma_skip = set()
    else:
        sigma_skip = set(config.sigma) - sigma_tier
    fst = Fst(config.symtable)

    # Final and penultimate state
    qf = ('λ',)
    qp = (eos,) + (epsilon,) * (context_length - 1)
    fst.add_state(qf)
    fst.add_state(qp)
    fst.add_arc(src=qp, ilabel=eos, dest=qf)

    # Interior transitions
    # xα -- x --> αy for each y
    Q = {qf, qp}
    Qnew = set(Q)
    for l in range(context_length + 1):
        Qold = set(Qnew)
        Qnew = set()
        for q2 in Qold:
            if q2 == qf:
                continue
            for x in sigma_tier:
                q1 = (x,) + prefix(q2, context_length - 1)
                fst.add_state(q1)
                fst.add_arc(src=q1, ilabel=x, dest=q2)
                Qnew.add(q1)
        Q |= Qnew

    # Initial state and outgoing transitions
    q0 = (bos,)
    fst.add_state(q0)
    for q in Q:
        if q == qf:
            continue
        fst.add_arc(src=q0, ilabel=bos, dest=q)
    Q.add(q0)

    # Self-transitions labeled by skipped symbols
    # on interior states
    for q in Q:
        if (q == q0) or (q == qf):
            continue
        for x in sigma_skip:
            fst.add_arc(src=q, ilabel=x, dest=q)

    fst.set_start(q0)
    fst.set_final(qf)
    #fst = fst.connect() # xxx handle state relabeling!
    return fst


def test():
    config.init(sigma_syms=['a', 'b'], special_syms=['λ'])
    fst = Fst(config.symtable)
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

    config.init(sigma_syms=['a', 'b'], special_syms=['λ'])
    L = left_context_acceptor(context_length=2)
    L.draw('L.dot')
    R = right_context_acceptor(context_length=2)
    R.draw('R.dot')


if __name__ == '__main__':
    test()