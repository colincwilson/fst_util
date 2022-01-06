# -*- coding: utf-8 -*-

import sys
from copy import copy
from pynini import SymbolTable

from . import fst_config  # xxx
from .fst import Fst


class SimpleFst():
    """
    Bare-bones unweighted FST implementation.
    """

    def __init__(self, Q=None, q0=None, F=None, T=None):
        self.Q = set(Q) if Q is not None else set()  # States
        self.q0 = q0 if q0 is not None else -1  # Initial state
        self.F = set(F) if F is not None else set()  # Final states
        self.T = set(T) if T is not None else set()  # Transitions

    # Operate on attributes directly
    # Add state
    # fst.Q.add(q)

    # Set initial state
    # fst.q0 = q

    # Add transition
    #fst.T.add(t)

    def delete_states(self, dead_states):
        # [nondestructive]
        fst = SimpleFst(
            Q=self.Q.difference(dead_states),
            q0=self.q0 if self.q0 not in dead_states else -1,
            F=self.F.difference(dead_states),
            T={
                t for t in self.T
                if t.src not in dead_states and t.dest not in dead_states
            })
        return fst

    def copy(self):
        Q = {q for q in self.Q}
        q0 = self.q0
        F = {q for q in self.F}
        T = {copy(t) for t in self.T}
        return SimpleFst(Q, q0, F, T)

    def print(self):
        val = f'Q {self.Q}\n'
        val += f'q0 {self.q0}\n'
        val += f'F {self.F}\n'
        val += f'T {[str(t) for t in self.T]}\n'
        return val

    def pynini(self):
        """
        Convert to state-labeled pynini FST
        """
        fst = Fst(SymbolTable(), SymbolTable())
        for q in self.Q:
            fst.add_state(q)
        fst.set_start(self.q0)
        for q in self.F:
            fst.set_final(q)
        for t in self.T:
            fst.add_arc(t.src, t.ilabel, t.olabel, None, t.dest)
        return fst


class Transition():
    """
    Transition of SimpleFST
    """

    def __init__(self, src, ilabel, olabel, dest):
        self.src = src
        self.ilabel = ilabel
        self.olabel = olabel
        self.dest = dest

    def __copy__(self):
        q, a, w, r = \
            self.src, self.ilabel, self.olabel, self.dest
        return Transition(copy(q), copy(a), copy(w), copy(r))

    def __eq__(self, other):
        return isinstance(other, self.__class__) \
            and (self.src == other.src) \
                and (self.ilabel == other.ilabel) \
                    and (self.olabel == other.olabel) \
                        and (self.dest == other.dest)

    #def __lt__(self, other):
    #    if not isinstance(other, self.__class__):
    #       raise Error('Incorrect type for Transition lt')
    #   return (self.src, self.ilabel, self.olabel, self.dest) < \
    #                (other.src, other.ilabel, other.olabel, other.dest)

    def __hash__(self):
        return hash((self.src, self.ilabel, self.olabel, self.dest))

    def __str__(self):
        return f'({self.src}, {self.ilabel}, {self.olabel}, {self.dest})'
