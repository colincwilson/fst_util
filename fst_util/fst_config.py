# -*- coding: utf-8 -*-

from pynini import SymbolTable

epsilon = 'ϵ'  # <eps>
bos = '⋊'  # '>' | <s>
eos = '⋉'  # '<' | </s>
syms = None  # All symbols in symtable
sigma = None  # Ordinary symbols
symtable = None


def init(sigma_syms, special_syms=[]):
    global syms, sigma, symtable
    sigma = sigma_syms[:]
    symtable = SymbolTable()
    symtable.add_symbol(epsilon)
    symtable.add_symbol(bos)
    symtable.add_symbol(eos)
    for sym in special_syms:
        symtable.add_symbol(sym)
    for sym in sigma:
        symtable.add_symbol(sym)
    syms = [sym for (sym_id, sym) in symtable]
    #print(syms)