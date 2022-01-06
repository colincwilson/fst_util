# -*- coding: utf-8 -*-

import os, re, sys
from pathlib import Path
from progress.bar import Bar as ProgressBar
import pynini

from . import fst_config as config
from .fst_config import *
from fst_util.fst import *
from fst_util.proc import *

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
    fst.draw('prefix_tree.dot')

    # Onward prefix tree
    print('[OSTIA] Make prefix tree onward')
    fst, _, _ = onward_tree(fst, fst.start(), λ)
    fst.draw('onward_tree.dot')

    # OSTIA proper
    print('[OSTIA] Merge states')
    print(f'|Q| = {fst.num_states()}')
    progress = ProgressBar('Processing', max=fst.num_states())
    q0 = fst.start()
    config.red_states = red_states = [q0]
    config.blue_states = blue_states = []
    for t in fst.arcs(q0):
        t_ilabel = fst.input_label(t.ilabel)
        if len(t_ilabel.split()) == 1:
            blue_states.append(t.nextstate)

    while len(blue_states) != 0:
        report(
            f'\nred_states = {config.red_states}, blue_states = {config.blue_states}',
            level=5)

        # Select blue state (arbitrary?)
        q = config.blue_states.pop(0)
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
            print(f'delete {q} = {fst.state_label(q)}')
            #fst.delete_arcs(q)
            fst = fst.delete_states([q])
            #print(fst.print(missing_sym='λ'))
            #sys.exit(0)

        # Update blue states
        config.blue_states = blue_states = []
        for p in red_states:
            for t in fst.arcs(p):
                q = t.nextstate
                if q in red_states or q in blue_states:
                    continue
                if re.search(eos, fst.state_label(q)):
                    continue
                blue_states.append(q)
        print('blue_states:', [fst.state_label(q) for q in blue_states])

        # Report progress
        progress.next()

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
    report(
        f'ostia_merge: q1 = {fst.state_label(q1)}, q2 = {fst.state_label(q2)}')
    fst = fst.copy()
    for q in fst.states():
        q_arcs = fst.mutable_arcs(q)
        for t in q_arcs:
            if t.nextstate == q2:
                t.nextstate = q1
                q_arcs.set_value(t)
    return ostia_fold(fst, q1, q2)


def ostia_fold(fst, q1, q2):
    """
    Subtree in q2 is "folded into" q1
    (de la Higuera, Algorithm 18.6)
    """
    fst = fst.copy()

    q1_out, q2_out = unk, unk
    for t1 in fst.arcs(q1):
        if fst.is_final(t1.nextstate):
            q1_out = fst.output_symbols().find(t1.olabel)
    for t2 in fst.arcs(q2):
        if fst.is_final(t2.nextstate):
            q2_out = fst.output_symbols().find(t2.olabel)

    w = ostia_outputs(q1_out, q2_out)
    report(f'w = {w}')
    if w is False:
        return False  # fail

    q1_arcs = fst.mutable_arcs(q1)
    for t1 in q1_arcs:
        if fst.is_final(t1.nextstate):
            t1.olabel = fst.mutable_output_symbols().add_symbol(w)
            q1_arcs.set_value(t1)

    for a in config.Sigma:
        report(f'\ta = {a}')
        t_new = None

        for t2 in fst.arcs(q2):
            if fst.input_label(t2.ilabel) != a:
                continue

            report(f'\tt2 = {t2}')
            t1_flag = False

            for t1 in fst.arcs(q1):
                if fst.input_label(t1.ilabel) != a:
                    continue

                report(f'\tt1 = {t1}')
                t1_flag = True

                t1_olabel = fst.output_label(t1.olabel)
                t2_olabel = fst.output_label(t2.olabel)
                #if t1_olabel != t2_olabel:  # de la Higuera, "due to loops"
                #    return False
                #if t1.dest in config.red_states:  # de la Higuera, errata (p. 6)
                #   return False
                if t1_olabel not in prefixes(
                        t2_olabel):  # Aksënova, pip release
                    return False
                #if t2_olabel not in prefixes(t1_olabel):
                #    return False  # Aksënova, github

                fst = ostia_pushback(fst, q1, q2, a)
                fst = ostia_fold(fst, t1.nextstate, t2.nextstate)
                if fst is False:  # Aksënova
                    return False

            if not t1_flag:
                t_new = {
                    'src': q1,
                    'ilabel': a,
                    'olabel': t2.olabel,
                    'dest': t2.nextstate
                }

        if t_new is not None:  # Aksënova
            fst.add_arc(
                src=t_new['src'],
                ilabel=t_new['ilabel'],
                olabel=t_new['olabel'],
                dest=t_new['dest'])

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

    # Transitions out of q1 and q2 with input label a
    # (transitions assumed to be unique)
    t1, t2 = None, None
    q1_arcs = fst.mutable_arcs(q1)
    q2_arcs = fst.mutable_arcs(q2)
    for t in q1_arcs:
        if fst.input_label(t.ilabel) == a:
            t1 = t
            break
    for t in q2_arcs:
        if fst.input_label(t.ilabel) == a:
            t2 = t
            break
    #report(f't1 = {t2}, t2 = {t1}')

    # Remove lcp from output labels
    t1_olabel = fst.output_label(t1.olabel)
    t2_olabel = fst.output_label(t2.olabel)
    u = lcp([t1_olabel, t2_olabel])
    #report(f'lcp({t1.olabel}, {t2.olabel}) = {u}')
    u1 = delete_prefix(t1_olabel, u)
    u2 = delete_prefix(t2_olabel, u)
    u = fst.mutable_output_symbols().add_symbol(u)
    t1.olabel = u
    t2.olabel = u
    q1_arcs.set_value(t1)
    q2_arcs.set_value(t2)

    # Push lcp onto labels of subsequent transitions
    # (also pushes lcp onto output strings)
    t1_next_arcs = fst.mutable_arcs(t1.nextstate)
    for t in t1_next_arcs:
        t_olabel = concat(u1, fst.output_label(t.olabel))
        t.olabel = fst.mutable_output_symbols().add_symbol(t_olabel)
        t1_next_arcs.set_value(t)
    t2_next_arcs = fst.mutable_arcs(t2.nextstate)
    for t in t2_next_arcs:
        t_olabel = concat(u2, fst.output_label(t.olabel))
        t.olabel = fst.mutable_output_symbols().add_symbol(t_olabel)
        t2_next_arcs.set_value(t)

    return fst