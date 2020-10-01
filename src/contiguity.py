# -*- coding: utf-8 -*-
"""
Functions related to the contiguity constraint

- calc_penalty_contig_ss
    returns the total number of violations of the contiguity constraint
    in a stacking sequence

- calc_penalty_contig_mp
    returns the total number of violations of the contiguity constraint
    by a multipanel structure

- penalty_contig
    returns the penalty value associated to the contiguity constraint

- is_contig
    returns True if a panel stacking sequence does not have a block
    with too many adjacent plies at the same fibre direction

- calc_n_viola_contig
    returns the number of violations of the contiguity constraint in each panel
    of a laminate design

- calc_contig_vector
    returns the vector of violations of the contiguity constraint for a
    stacking sequence

- calc_matrix_viol_contig
    returns the matrix of violations of the contiguity constraint by a
    multipanel structure

    viola_contig[ind_panel, ind_ply] == 0 when the plies of indices
    {ind_ply, ind_ply + 1, ..., ind_ply + n_contig} does not violates the
    contiguity constraint

    viola_contig[ind_panel, ind_ply] == 2 when the plies of indices
    {ind_ply, ind_ply + 1, ..., ind_ply + n_contig} violate the contiguity
    constraint and also the plies of indices
    {ind_ply - 1, ind_ply, ..., ind_ply + n_contig - 1}

    viola_contig[ind_panel, ind_ply] == 1 when the plies of indices
    {ind_ply, ind_ply + 1, ..., ind_ply + n_contig} violate the contiguity
    constraint but not the plies of indices
    {ind_ply - 1, ind_ply, ..., ind_ply + n_contig - 1}
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints

def penalty_contig(n, viola_contig=np.array(()), contig=False,
                   coeff_contig=1, norm_diso_contig=1):
    """
    returns the penalty value associated to the disorientation and contiguity
    constraints

    INPUTS

    - n: number of stacking sequences tested
    - contig if the contiguity constraint is active
    - viola_contig: matrix of violations of the contiguity constraint
    - coeff_diso: weight of the penalty for the cdisorientation constraint
    - norm_diso_contig: sum of the ply counts of each panel related to the
    current group of plies
    """
    if contig:
        return (coeff_contig * np.sum(
            viola_contig, axis=1)) / norm_diso_contig
    return np.zeros((n,), int)

def is_contig(ss, n_contig):
    """
    returns True if the stacking sequence ss does not have a block of more
    than n_contig adjacent plies at the same fibre direction

    INPUTS

    - ss: stacking sequence (array of ints)
    - n_contig: maximum number of adajcent plies with the same fibre
    orientation allowed (int)
    """
    n_plies = ss.size
    if n_contig < n_plies:
        diff = n_plies - n_contig
        if n_contig == 2:
            for ind in np.arange(diff):
                if ss[ind] == ss[ind + 1] \
                and ss[ind + 2] == ss[ind + 1]:
                    return False
        elif n_contig == 3:
            for ind in np.arange(diff):
                if ss[ind] == ss[ind + 1] \
                and ss[ind + 2] == ss[ind + 1] \
                and ss[ind + 3] == ss[ind + 1]:
                    return False
        elif n_contig == 4:
            for ind in np.arange(diff):
                if ss[ind] == ss[ind + 1] \
                and ss[ind + 2] == ss[ind + 1] \
                and ss[ind + 2] == ss[ind + 3] \
                and ss[ind] == ss[ind + 4]:
                    return False
        elif n_contig == 5:
            for ind in np.arange(diff):
                if ss[ind] == ss[ind + 1] \
                and ss[ind] == ss[ind + 2] \
                and ss[ind] == ss[ind + 3] \
                and ss[ind] == ss[ind + 4] \
                and ss[ind] == ss[ind + 5]:
                    return False
        elif n_contig == 6:
            for ind in np.arange(diff):
                if ss[ind] == ss[ind + 1] \
                and ss[ind] == ss[ind + 2] \
                and ss[ind] == ss[ind + 3] \
                and ss[ind] == ss[ind + 4] \
                and ss[ind] == ss[ind + 5] \
                and ss[ind] == ss[ind + 6]:
                    return False
        elif n_contig != 1:
            raise Exception('n_contig must be 2, 3, 4 or 5.')
    return True


def calc_n_viola_contig(
        mother_n_viola_contig, mother_ss_top, child_ss, n_panels, constraints,
        level, n_plies, pdl, middle_ply=False, reflexion=False,
        mother_ss_bot=None, last_level=False, odd_level=False):

    """
    returns the number of violations of the contiguity constraint in each panel
    of a laminate design

    INPUTS

    - mother_n_viola_contig: matrix of violations of the contiguity
    constraints by the mother node
    - mother_ss: mother stacking sequences
    - child_ss: possible fibre orientation for the new ply
    - level: level in the search tree
    - n_plies: number of plies in the group for the thickest panel
    - pdl: matrix of ply drops
    - n_panels: number of panels in the laminate structure
    - constraints: set of constraints
    - reflexion: true if the group is at the middle of a symmetric laminate
    - mother_ss_top: top stacking sequences of the incomplete laminate design
    - mother_ss_bot: bottom stacking sequence  of the incomplete laminate
    design
    - odd_level = True for a level of the beam search corresponding to a ply in
    the top of the laminate
    """
    n_viola_contig = np.copy(mother_n_viola_contig)

    ### last ply in asymmetric laminate
    if not constraints.sym and last_level:
        for ind_panel in range(n_panels):
            if pdl[ind_panel, level - 1] >= 0:
                new_stack = np.hstack((
                    mother_ss_top[ind_panel][-constraints.n_contig:],
                    np.array([child_ss[0]]),
                    mother_ss_bot[ind_panel][:constraints.n_contig]))
            else:
                new_stack = np.hstack((
                    mother_ss_top[ind_panel][-constraints.n_contig:],
                    mother_ss_bot[ind_panel][:constraints.n_contig]))
            vector = calc_contig_vector(new_stack, constraints.n_contig)
            n_viola_contig[ind_panel] += np.sum(vector)
        return n_viola_contig

    ### top
    if constraints.sym or not odd_level:

        if reflexion and last_level:
            for ind_panel in range(n_panels):
                if child_ss.size == 2 and pdl[ind_panel, level - 1] >= 0 \
                and pdl[ind_panel, level] >= 0:
                    new_stack = np.hstack((
                        mother_ss_top[ind_panel][-constraints.n_contig:],
                        np.array([child_ss[0]]),
                        np.array([child_ss[1]]),
                        np.array([child_ss[0]]),
                        np.flip(mother_ss_top[
                            ind_panel][-constraints.n_contig:], axis=0)))
                elif pdl[ind_panel, level - 1] >= 0:
                    new_stack = np.hstack((
                        mother_ss_top[ind_panel][-constraints.n_contig:],
                        np.array([child_ss[0]]),
                        np.array([child_ss[0]]),
                        np.flip(mother_ss_top[
                            ind_panel][-constraints.n_contig:], axis=0)))
                elif child_ss.size == 2 and pdl[ind_panel, level] >= 0 :
                    new_stack = np.hstack((
                        mother_ss_top[ind_panel][-constraints.n_contig:],
                        np.array([child_ss[1]]),
                        np.flip(mother_ss_top[
                            ind_panel][-constraints.n_contig:], axis=0)))
                else:
                    new_stack = np.hstack((
                        mother_ss_top[ind_panel][-constraints.n_contig:],
                        np.flip(mother_ss_top[
                            ind_panel][-constraints.n_contig:], axis=0)))
                vector = calc_contig_vector(new_stack, constraints.n_contig)
                n_viola_contig[ind_panel] += np.sum(vector)

        else:
            for ind_panel in range(n_panels):
                if pdl[ind_panel, level - 1] >= 0 \
                and  mother_ss_top[ind_panel].size > 0:
                    new_stack = np.hstack((
                        mother_ss_top[ind_panel][-constraints.n_contig:],
                        np.array([child_ss[0]])))
                    if not is_contig(new_stack, constraints.n_contig):
                        if constraints.sym:
                            n_viola_contig[ind_panel] += 2
                        else:
                            n_viola_contig[ind_panel] += 1

    else: ### bottom
        for ind_panel in range(n_panels):
            if pdl[ind_panel, level - 1] >= 0 \
            and mother_ss_bot[ind_panel].size > 0:
                new_stack = np.hstack((
                    np.array([child_ss[0]]),
                    mother_ss_bot[ind_panel][:constraints.n_contig]))
                vector = calc_contig_vector(new_stack, constraints.n_contig)
                if not is_contig(new_stack, constraints.n_contig):
                    n_viola_contig[ind_panel] += 1

    return n_viola_contig


def calc_contig_vector(ss, n_contig):
    """
    returns the vector of violations of the contiguity constraint for a
    stacking sequence
    """
    step = ss.size - n_contig
    if step > 0:
        result = np.zeros((step,), dtype='int16')
        for ind in range(step):
            if is_contig(ss[ind:ind + n_contig + 1], n_contig): # if violation
                pass
            else:
                result[ind] = 1
        return result
    return np.array(())

def calc_matrix_viol_contig(ss, n_contig):
    """
    returns the matrix of violations of the contiguity constraint by a multipanel
    structure

    - viola_contig[ind_panel, ind_ply] == 0 when the plies of indices
    {ind_ply, ind_ply + 1, ..., ind_ply + n_contig} does not violates the
    contiguity constraint

    - viola_contig[ind_panel, ind_ply] == 2 when the plies of indices
    {ind_ply, ind_ply + 1, ..., ind_ply + n_contig} violate the contiguity
    constraint and also the plies of indices
    {ind_ply - 1, ind_ply, ..., ind_ply + n_contig - 1}

    - viola_contig[ind_panel, ind_ply] == 1 when the plies of indices
    {ind_ply, ind_ply + 1, ..., ind_ply + n_contig} violate the contiguity
    constraint but not the plies of indices
    {ind_ply - 1, ind_ply, ..., ind_ply + n_contig - 1}

    INPUTS

    - ss: stacking sequences of each panel of the structure (list or arrays)
    - n_contig: maximum of adjacent plies that can have the same fibre orientaion
    """
    viola_contig = [[]]*len(ss)
    for ind_panel, ss_panel in enumerate(ss):
        viola_contig[ind_panel] = np.zeros((ss_panel.size,), dtype=bool)
        before = 0
        for ind_ply in range(ss_panel.size):
            if not is_contig(
                    ss_panel[ind_ply: ind_ply + n_contig + 1], n_contig):
                if before == 0:
                    viola_contig[ind_panel][ind_ply] += 1
                    before = 1
                else:
                    viola_contig[ind_panel][ind_ply] += 2
                    before = 1
            else:
                before = 0
    return viola_contig


def calc_penalty_contig_mp(ss, constraints):
    """
    returns the total number of violations of the contiguity constraint by a
    multipanel structure

    INPUTS

    - ss: stacking sequences of each panel of the structure (list or arrays)
    - constraints: set of constraints

    """
    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule'):
            if constraints.dam_tol_rule in {2, 3}:
                ss = [el[1:-1] for el in ss]
        else:
            if constraints.n_plies_dam_tol == 2:
                ss = [el[1:-1] for el in ss]

    viola_contig = calc_matrix_viol_contig(ss, constraints.n_contig)
    return np.array([np.sum(el) for el in viola_contig]) # sum per panel

def calc_penalty_contig_ss(ss, constraints):
    """
    returns the total number of violations of the contiguity constraint in
    a stacking sequence

    INPUTS

    - ss: panel stacking sequence
    - constraints: set of constraints
    """
    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule'):
            if constraints.dam_tol_rule in {2, 3}:
                ss = np.copy(ss)[1:-1]
        else:
            if constraints.n_plies_dam_tol == 2:
                ss = np.copy(ss)[1:-1]

    return np.sum(calc_matrix_viol_contig([ss], constraints.n_contig))


if __name__ == "__main__":
    print('*** Test for the function is_contig ***\n')
#    print('Inputs:\n')
#    ss = np.array([0, 45, 45, 45, 45, 90, -45])
#    n_contig = 3
#    print(f'ss = {ss}, n_contig = {n_contig}\n')
#    print('outputs:\n')
#    print(is_contig(ss, n_contig))

    print('*** Test for the function calc_contig_vector ***\n')
    print('Inputs:\n')
    ss = np.array([0, 45, 45, 45, 45, 90, 90, 90, 90, 90, -45])
    n_contig = 3
    print(f'ss = {ss}, n_contig = {n_contig}\n')
    print('outputs:\n')
    print(calc_contig_vector(ss, n_contig))

    print('*** Test for the function calc_matrix_viol_contig ***\n')
    print('Inputs:\n')
    ss = [np.array([0, 90, 90]), np.array([0, 0, 0, 0, 0, 0])]
    n_contig = 3
    print(f'ss: {ss}')
    print(f'n_contig: {n_contig}')
    print('outputs:\n')
    viola_contig = calc_matrix_viol_contig(ss, n_contig)
    print(f'viola_contig: {viola_contig}\n')

    print('*** Test for the function calc_penalty_contig_mp ***\n')
    print('Inputs:\n')
    constraints = Constraints(sym=True)
    constraints.dam_tol = True
    constraints.n_contig = 3
    ss = [np.array([0, 90, 90]), np.array([0, 0, 0, 0, 0, 0])]
    print(f'ss: {ss}')
    print('outputs:\n')
    viola_contig = calc_penalty_contig_mp(ss, constraints)
    print(f'Result: {viola_contig}\n')

    print('*** Test for the function calc_penalty_contig_ss ***\n')
    constraints = Constraints(sym=True)
    constraints.dam_tol = True
    constraints.n_contig = 3
    print('Inputs:\n')
    ss = np.array([0, 90, 90, 90, 90])
    print(f'ss: {ss}')
    print(f'n_contig: {constraints.n_contig}')
    print('outputs:\n')
    print(f'Result: {calc_penalty_contig_ss(ss, constraints)}\n')
