# -*- coding: utf-8 -*-

import sys
import pynini
from . import fst_config  # xxx


class Fst(pynini.Fst):
    """
    Pynini Fst with labeled inputs/outputs/states
    - Input/output labels must be strings
    - State labels must be hashable (strings, tuples, etc.)
    todo: deepcopy; destructive operations
    """

    def __init__(self, symtable):
        super(Fst, self).__init__()
        super(Fst, self).set_input_symbols(symtable)
        super(Fst, self).set_output_symbols(symtable)
        self._state2label = {}  # State id -> label
        self._label2state = {}  # Label -> state id

    def add_state(self, state_label=None):
        """ Add new state, optionally specifying its label """
        # Enforce unique labels
        if state_label is not None:
            if state_label in self._label2state:
                return self._label2state[state_label]
        # Create new state
        state = super(Fst, self).add_state()
        # Self-labeling by default
        if state_label is None:
            state_label = state
        # State <-> label
        self._state2label[state] = state_label
        self._label2state[state_label] = state
        return state

    def add_arc(self, src, ilabel, olabel=None, weight=None, dest=None):
        """ Add arc (accepts int or string attributes) """
        if not isinstance(src, int):
            src = self._label2state[src]
        if not isinstance(dest, int):
            dest = self._label2state[dest]
        if not isinstance(ilabel, int):
            ilabel = self.input_symbols().find(ilabel)
        if olabel is None:
            olabel = ilabel
        elif not isinstance(olabel, int):
            olabel = self.output_symbols().find(olabel)
        if weight is None:
            weight = 0
        #print(src, ilabel, olabel, weight, dest)
        arc = pynini.Arc(ilabel, olabel, weight, dest)
        return super(Fst, self).add_arc(src, arc)

    def set_start(self, state):
        if not isinstance(state, int):
            state = self._label2state[state]
        return super(Fst, self).set_start(state)

    def set_final(self, state, weight=None):
        if not isinstance(state, int):
            state = self._label2state[state]
        return super(Fst, self).set_final(state, weight)

    def input_label(self, sym_id):
        return self.input_symbols().find(sym_id)

    def output_label(self, sym_id):
        return self.output_symbols().find(sym_id)

    def state_label(self, state):
        return self._state2label[state]

    def arcs(self, src):
        if not isinstance(src, int):
            src = self._label2state[src]
        return super(Fst, self).arcs(src)

    def mutable_arcs(self, src):
        if not isinstance(src, int):
            src = self._label2state[src]
        return super(Fst, self).mutable_arcs(src)

    def connect(self):
        """
        Remove states and arcs not on successful paths [nondestructive]
        """
        accessible = self.accessible(forward=True)
        coaccessible = self.accessible(forward=False)
        live_states = accessible & coaccessible
        dead_states = filter(lambda q: q not in live_states, self.states())
        fst = self._delete_states(dead_states)
        return fst

    def accessible(self, forward=True):
        """
        States accessible from initial state -or- coaccessible from final states
        """
        if forward:
            # Initial state and forward transitions
            Q = set([self.start()])
            T = {}
            for src in self.states():
                T[src] = set()
                for t in self.arcs(src):
                    dest = t.nextstate
                    T[src].add(dest)
        else:
            # Final states and backward transitions
            Zero = pynini.Weight.zero(self.weight_type())
            Q = set([q for q in self.states() if self.final(q) != Zero])
            T = {}
            for src in self.states():
                for t in self.arcs(src):
                    dest = t.nextstate
                    if dest not in T:
                        T[dest] = set()
                    T[dest].add(src)

        # Find (co)accessible states
        Q_old = set()
        Q_new = set(Q)
        while len(Q_new) != 0:
            Q_old, Q_new = Q_new, Q_old
            Q_new.clear()
            for src in filter(lambda src: src in T, Q_old):
                for dest in filter(lambda dest: dest not in Q, T[src]):
                    Q.add(dest)
                    Q_new.add(dest)
        return Q

    def delete_states(self, dead_states):
        """
        Remove states while preserving labels, trim result [nondestructive]
        """
        fst = self._delete_states(dead_states)
        return fst.connect()

    def _delete_states(self, dead_states):
        """
        Remove states while preserving labels [nondestructive]
        """
        fst = Fst(self.input_symbols())
        # Reindex live states, copying labels
        state_map = {}
        q0 = self.start()
        for q in filter(lambda q: q not in dead_states, self.states()):
            q_label = self._state2label[q]
            q_new = fst.add_state(q_label)
            state_map[q] = q_new
            if q == q0:
                fst.set_start(q_new)
            fst.set_final(q_new, self.final(q))
        # Copy transitions between live states
        for q in filter(lambda q: q not in dead_states, self.states()):
            src = state_map[q]
            for t in filter(lambda t: t.nextstate not in dead_states,
                            self.arcs(q)):
                dest = state_map[t.nextstate]
                fst.add_arc(src, t.ilabel, t.olabel, t.weight, dest)
        return fst

    def delete_arcs(self, dead_arcs):
        """
        Delete arcs [destructive]
        see: http://www.openfst.org/twiki/bin/view/Forum/FstForumArchive2014#Deleting%20a%20specific%20arc%20in%20an%20FS
        """
        # Group dead arcs by source state xxx do outside
        dead_arcs_ = {}
        for (src, t) in dead_arcs:
            if src in dead_arcs_:
                dead_arcs_[src].append(t)
            else:
                dead_arcs_[src] = [t]
        # Process states with dead arcs
        for q in filter(lambda q: q in dead_arcs_, self.states()):
            # Remove all arcs from state
            arcs = [t for t in self.arcs(q)]
            super(Fst, self).delete_arcs(q)
            # Add back live arcs
            for t1 in arcs:
                live = True
                for t2 in dead_arcs_[q]:
                    if arc_equal(t1, t2):  # xxx
                        live = False
                        break
                if live:
                    self.add_arc(q, t1.ilabel, t1.olabel, t1.weight,
                                 t1.nextstate)
        return self

    def num_arcs(self):
        """
        Total count of arcs from all states
        """
        val = 0
        for state in self.states():
            val += super(Fst, self).num_arcs(state)
        return val

    def print(self, **kwargs):
        # Stringify state labels
        ssymbols = pynini.SymbolTable()
        for q, label in self._state2label.items():
            ssymbols.add_symbol(str(label), q)
        return super(Fst, self).print(
            isymbols=self.input_symbols(),
            osymbols=self.output_symbols(),
            ssymbols=ssymbols,
            **kwargs)

    def draw(self, source, acceptor=True, portrait=True, **kwargs):
        # Stringify state labels
        ssymbols = pynini.SymbolTable()
        for q, label in self._state2label.items():
            ssymbols.add_symbol(str(label), q)
        return super(Fst, self).draw(
            source,
            isymbols=self.input_symbols(),
            osymbols=self.output_symbols(),
            ssymbols=ssymbols,
            acceptor=acceptor,
            portrait=portrait,
            **kwargs)


def arc_equal(arc1, arc2):
    """
    Arc equality (missing from pynini?)
    """
    val = (arc1.ilabel == arc2.ilabel) and \
        (arc1.olabel == arc2.olabel) and \
            (arc1.weight == arc2.weight) and \
                (arc1.nextstate == arc2.nextstate)
    return val


def compose(fst1, fst2):
    """
    FST composition, retaining contextual info from original machines by labeling each state q = (q1, q2) with (label(q1), label(q2))
    todo: matcher options
    """
    fst = Fst(fst_config.symtable)
    One = pynini.Weight.one(fst.weight_type())
    # xxx arcsort(), mutable_arcs(), final(), start(), states()

    q0_1 = fst1.start()
    q0_2 = fst2.start()
    q0 = (fst1.state_label(q0_1), fst2.state_label(q0_2))
    fst.add_state(q0)
    fst.set_start(q0)

    # Lazy state and transition construction
    Q = set([q0])
    Q_old, Q_new = set(), set([q0])
    while len(Q_new) != 0:
        Q_old, Q_new = Q_new, Q_old
        Q_new.clear()
        for src in Q_old:
            src1, src2 = src  # State labels in M1, M2
            for t1 in fst1.arcs(src1):
                for t2 in fst2.arcs(src2):  # xxx use arc sort
                    if t1.olabel != t2.ilabel:
                        continue
                    dest1 = t1.nextstate
                    dest2 = t2.nextstate
                    dest = (fst1.state_label(dest1), fst2.state_label(dest2))
                    fst.add_state(dest)  # No change if state already exists
                    fst.add_arc(
                        src=src, ilabel=t1.ilabel, olabel=t2.olabel, dest=dest)
                    if fst1.final(dest1) == fst2.final(dest2) == One:  # xxx
                        fst.set_final(dest)
                    if dest not in Q:
                        Q.add(dest)
                        Q_new.add(dest)

    return fst.connect()


def accepted_strings(fst, side='output', max_len=10):
    """
    Strings accepted by fst on designated side, up to max_len (not including bos/eos); cf. pynini for paths through acyclic fst
    todo: handle epsilons
    """
    q0 = fst.start()
    Zero = pynini.Weight.zero(fst.weight_type())

    accepted = set()
    prefixes_old = {(q0, None)}
    prefixes_new = set(prefixes_old)
    for i in range(max_len + 2):
        prefixes_old, prefixes_new = prefixes_new, prefixes_old
        prefixes_new.clear()
        for (src, prefix) in prefixes_old:
            for t in fst.arcs(src):
                dest = t.nextstate
                if side == 'input':
                    tlabel = fst.input_label(t.ilabel)
                else:
                    tlabel = fst.output_label(t.olabel)
                if prefix is None:
                    prefix_ = tlabel
                else:
                    prefix_ = prefix + ' ' + tlabel
                prefixes_new.add((dest, prefix_))
                if fst.final(dest) != Zero:
                    accepted.add(prefix_)
                    #print(prefix_)

    return accepted


def left_context_acceptor(context_length=1, sigma_tier=None):
    """
    Acceptor (identity transducer) for segments in immediately precedin contexts (histories) of specified length. If Sigma_tier is specified as  a subset of Sigma, only contexts over Sigma_tier are tracked (other member of Sigma are skipped, i.e., label self-loops on each interior state)
    """
    epsilon = fst_config.epsilon
    bos = fst_config.bos
    eos = fst_config.eos
    if sigma_tier is None:
        sigma_tier = set(fst_config.sigma)
        sigma_skip = set()
    else:
        sigma_skip = set(fst_config.sigma) - sigma_tier
    fst = Fst(fst_config.symtable)

    # Initial and peninitial states
    q0 = ('λ',)
    q1 = (epsilon,) * (context_length - 1) + (bos,)
    fst.add_state(q0)
    fst.set_start(q0)
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
                q2 = _suffix(q1, context_length - 1) + (x,)
                fst.add_state(q2)
                fst.add_arc(src=q1, ilabel=x, dest=q2)
                Qnew.add(q2)
        Q |= Qnew

    # Final state and incoming arcs
    qf = (eos,)
    fst.add_state(qf)
    fst.set_final(qf)
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

    #fst = fst.connect() # xxx handle state relabeling
    return fst


def right_context_acceptor(context_length=1, sigma_tier=None):
    """
    Acceptor (identity transducer) for segments in immediately following contexts (futures) of specified length. If Sigma_tier is specified as a subset of Sigma, only contexts over Sigma_tier are tracked (other members of Sigma are skipped, i.e., label self-loops on each interior state)
    """
    epsilon = fst_config.epsilon
    bos = fst_config.bos
    eos = fst_config.eos
    if sigma_tier is None:
        sigma_tier = set(fst_config.sigma)
        sigma_skip = set()
    else:
        sigma_skip = set(fst_config.sigma) - sigma_tier
    fst = Fst(fst_config.symtable)

    # Final and penultimate state
    qf = ('λ',)
    qp = (eos,) + (epsilon,) * (context_length - 1)
    fst.add_state(qf)
    fst.set_final(qf)
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
                q1 = (x,) + _prefix(q2, context_length - 1)
                fst.add_state(q1)
                fst.add_arc(src=q1, ilabel=x, dest=q2)
                Qnew.add(q1)
        Q |= Qnew

    # Initial state and outgoing transitions
    q0 = (bos,)
    fst.add_state(q0)
    fst.set_start(q0)
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

    #fst = fst.connect() # xxx handle state relabeling
    return fst


def _prefix(x, l):
    """ Length-l prefix of tuple x """
    if l < 1:
        return ()
    if len(x) < l:
        return x
    return x[:l]


def _suffix(x, l):
    """ Length-l suffix of tuple x """
    if l < 1:
        return ()
    if len(x) < l:
        return x
    return x[-l:]
