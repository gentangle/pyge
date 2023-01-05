"""
Functions to compute the fraction of entangled contacts

The function takes as argument the list of contacts together with info
about the respective thread and G' value. This list is the one
produced by gent.ge_loops
"""
import numpy as np

# funzione che prende la lista gia calcolata
# gradi di liberta': funzione di attivazine, G0, w

def hill_fun(inp, threshold, hill_coeff):
    return 1/( 1 + (threshold/inp)**hill_coeff )

def qentangled(contact_thr_ge_list, activation_fun, activation_params):
    pass
    