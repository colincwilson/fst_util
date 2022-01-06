# -*- coding: utf-8 -*-

import os, re, sys
from pathlib import Path
from progress.bar import Bar as ProgressBar

from . import fst_config as config
from .fst_config import *
from fst_util.simple_fst import *
from fst_util.proc2 import *

config.verbosity = 10


def ostia(D, Sigma=None, Lambda=None):
    """
    OSTIA (de la Higuera, Algorithm 18.7)
    """
    # Input and output alphabets
    print('[OSTIA] Input and output alphabets')
    if Sigma is None or Lambda is None:
        Sigma, Lambda = set([]), set([])
        for (x, y) in D:
            Sigma |= set(x.split(' '))
            Lambda |= set(y.split(' '))
    config.Sigma = Sigma
    config.Lambda = Lambda
    print(f'Sigma: {Sigma}')
    print(f'Lambda: {Lambda}')

    # Prefix tree
    print('[OSTIA] Build prefix tree')
    fst, _, _ = prefix_tree(D, Sigma, Lambda)
    fst.pynini().draw('prefix_tree.dot')

    # Onward prefix tree
    print('[OSTIA] Make prefix tree onward')
    fst = onward_tree(fst, fst.q0, λ)
    fst.pynini().draw('onward_tree.dot')

    # OSTIA proper
    print('[OSTIA] Merge states')
    print(f'|Q| = {len(fst.Q)}')
    progress = ProgressBar('Processing', max=len(fst.Q))
    red_states = [fst.q0]
    blue_states = set()
    for t in filter(lambda t: t.src == fst.q0, fst.T):
        if len(t.ilabel.split()) == 1:
            blue_states.add(t.dest)
    blue_states = list(blue_states)
    blue_states.sort(key=lambda q: (len(q.split()), q))

    while len(blue_states) != 0:
        report(
            f'\nred_states = {red_states}, blue_states = {blue_states}',
            level=5)

        # Select blue state (arbitrary?)
        q = blue_states.pop(0)
        # Attempt to merge blue state with each red state
        merge_flag = False
        for p in red_states:
            fst_old = fst.copy()
            fst = ostia_merge(fst, p, q)
            if fst:  # Merge succeeded, keep modified transducer
                report('merge accepted', level=5)
                merge_flag = True
                break
            else:  # Merge failed, undo changes to transducer
                fst = fst_old
                report('merge rejected', level=5)

        # All merge attempts failed
        if not merge_flag:
            if q not in red_states:
                red_states.append(q)

        # Merge attempt succeeded
        else:
            print(f'delete {q}')
            fst = fst.delete_states([q])
            #print(fst.print(missing_sym='λ'))
            #sys.exit(0)

        # Update blue states
        blue_states = set()
        for t in filter(lambda t: t.src in red_states, fst.T):
            q = t.dest
            if q in red_states:
                continue
            if re.search(eos, q):  # Keep all final states
                continue
            blue_states.add(q)
        blue_states = list(blue_states)
        blue_states.sort(key=lambda q: (len(q.split()), q))

        # Report progress
        progress.next()

    fst = fst.pynini()
    fst = fst.connect()
    print(fst.print(missing_sym='λ'))
    fst.draw('ostia.dot')
    print()
    print(f'|Q| {fst.num_states()}')
    print(f'|T| {fst.num_arcs()}')
    return fst


def ostia_merge(fst, q1, q2):
    """
    Replace q2 destination states with q1
    "... merge the states and, at the same time, ensure that the result is onward. The merging algorithm is a merge-and-fold variant: it first computes the longest common prefix of every two outputs it is going to have to merge, and then makes the necessary merges."
    (de la Higeura, Algorithm 18.5)
    """
    report(f'ostia_merge: q1 = {q1}, q2 = {q2}')
    fst = fst.copy()
    fst.T = list(fst.T)  # Modifiable transitions
    for t in fst.T:
        if t.dest == q2:
            t.dest = q1
    return ostia_fold(fst, q1, q2)


def ostia_fold(fst, q1, q2):
    """
    Subtree in q2 is "folded into" q1
    (de la Higuera, Algorithm 18.6)
    """
    fst = fst.copy()
    fst.T = list(fst.T)  # Modifiable transitions

    q1_out, q2_out = unk, unk
    for t in fst.T:
        if t.src == q1 and t.dest in fst.F:
            q1_out = t.olabel
        if t.src == q2 and t.dest in fst.F:
            q2_out = t.olabel

    w = ostia_outputs(q1_out, q2_out)
    report(f'w = {w}')
    if w is False:
        return False  # fail

    for t in fst.T:
        if t.src == q1 and t.dest in fst.F:
            t.olabel = w

    for a in config.Sigma:
        report(f'\ta = {a}')
        t_new = None

        for t2 in filter(lambda t: t.src == q2 and t.ilabel == a, fst.T):
            report(f'\tt2 = {t2}')
            t1_flag = False

            for t1 in filter(lambda t: t.src == q1 and t.ilabel == a, fst.T):
                report(f'\tt1 = {t1}')
                t1_flag = True

                #if t1_olabel != t2_olabel:  # de la Higuera, "due to loops"
                #    return False
                #if t1.dest in config.red_states:  # de la Higuera, errata (p. 6)
                #   return False
                if t1.olabel not in prefixes(
                        t2.olabel):  # Aksënova, pip release
                    return False
                #if t2.olabel not in prefixes(t1.olabel):
                #    return False  # Aksënova, github

                fst = ostia_pushback(fst, q1, q2, a)
                fst = ostia_fold(fst, t1.dest, t2.dest)
                if fst is False:  # Aksënova
                    return False

            if not t1_flag:
                t_new = Transition(q1, a, t2.olabel, t2.dest)

        if t_new is not None:  # Aksënova
            fst.T.add(t_new)

    return fst


def ostia_outputs(w1, w2):
    """
    Check if two state outputs are equal modulo ⊥ (unknown)
    (de la Higuera, Algorithm 18.3)
    """
    report(f'ostia_outputs: {w1}, {w2}', level=5)
    if w1 == unk:
        return w2
    if w2 == unk:
        return w1
    if w1 == w2:
        return w1
    return False  # fail


def ostia_pushback(fst, q1, q2, a):
    """
    "... takes as inputs two transitions both labelled by the same symbol a (one starting in state q1 and the other in state q2) and returns a transducer equivalent to the initial one in which the output has been unified. This is done by _pushing back_ whatever uncommon suffixes the algorithm finds. There are cases where, due to loops, the result is that we don't have tau1(q1,a) = tau1(q2,a)."
    "Typically unk just absorbs any pushed-back suffix."
    (de la Higeura, Algorithm 18.4)
    """
    report(f'ostia_pushback: q1 = {q1}, q2 = {q2}, a = {a}', level=5)
    fst = fst.copy()
    fst.T = list(fst.T)  # Modifiable transitions

    # Transitions out of q1 and q2 with input label a
    # (transitions assumed to be unique)
    t1, t2 = None, None
    for t in fst.T:
        if t.src == q1 and t.ilabel == a:
            t1 = t
        if t.src == q2 and t.ilabel == a:
            t2 = t
    #report(f't1 = {t2}, t2 = {t1}')

    # Remove lcp from output labels
    u = lcp([t1.olabel, t2.olabel])
    #report(f'lcp({t1.olabel}, {t2.olabel}) = {u}')
    u1 = delete_prefix(t1.olabel, u)
    u2 = delete_prefix(t2.olabel, u)
    t1.olabel = u
    t2.olabel = u

    # Push lcp onto labels of subsequent transitions
    # (also pushes lcp onto output strings)
    for t in filter(lambda t: t.src == t1.dest, fst.T):
        t.olabel = concat(u1, t.olabel)
    for t in filter(lambda t: t.src == t2.dest, fst.T):
        t.olabel = concat(u2, t.olabel)

    return fst