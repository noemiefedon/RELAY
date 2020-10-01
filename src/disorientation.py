# -*- coding: utf-8 -*-
"""
Functions related to the disorientation constraint

- calc_n_penalty_diso_ss
    returns the number of disorientation constraint violations in a stacking
    sequence

- calc_penalty_diso
    returns the penalty value associated to the disorientation constraint

- is_diso
    returns True if the difference of angle between the two fibre orientations
    is smaller than or equal to a limit value  or if one of the angles is
    missing, False otherwise

- is_diso_ss
    returns True if the differences of angle between the adajancent angles are
    all smaller than or equal to delta_angle, False otherwise

- calc_n_viola_diso
    returns the number of violations of the disorientation constraint in each
    panel of a laminate design

- calc_penalty_diso_mp_0
    returns the matrix of violations of the disorientation constraint by a
    multipanel structure

    viola_diso[ind_panel, index_ply] = 1 when the plies of indices
    {index_ply, index_ply + 1} violates the constraint, 0 otherwise

- calc_penalty_diso_mp_1
    returns the number of disorientation constraint violations per panel in a
    multipanel structure
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
sys.path.append(r'C:\RELAY')
from src.constraints import Constraints

def calc_penalty_diso(n, viola_diso=np.array(()), diso=False, coeff_diso=1,
                      norm_diso_contig=1):
    """
    returns the penalty value associated to the disorientation and contiguity
    constraints

    INPUTS

    - n: number of stacking sequences tested
    - diso if the disorientation constraint is active
    - viola_diso: matrix of violations of the disorientation constraint
    - coeff_diso: weight of the penalty for the cdisorientation constraint
    - norm_diso_contig: sum of the ply counts of each panel related to the
    current group of plies
    """
    if diso:
        return (coeff_diso*np.sum(viola_diso, axis=1))/norm_diso_contig
    return np.zeros((n,), int)


def is_diso(angle1, angle2, delta_angle):
    """
    returns True if the difference of angle between angle1 and angle2 is
    smaller than or equal to delta_angle or if one angle is missing,
    False otherwise

    The angles are considered in degrees.
    """
    if abs(angle1 - angle2) > delta_angle and \
    abs(180 + angle1 - angle2) > delta_angle and \
    abs(angle1 - angle2 - 180) > delta_angle:
        return False
    return True


def is_diso_ss(ss, delta_angle, dam_tol=False,
               dam_tol_rule=None, n_plies_dam_tol=None):
    """
    returns True if the differences of angle between the adajancent angles are
    all smaller than or equal to delta_angle, False otherwise

    The angles are considered in degrees.
    """
    if dam_tol:
        if (dam_tol_rule is not None and dam_tol_rule in {2, 3})\
        or (n_plies_dam_tol is not None and n_plies_dam_tol == 2):
            for ind_ply in range(ss.size - 3):
                if not is_diso(ss[ind_ply + 1], ss[ind_ply + 2], delta_angle):
                    return False
        if dam_tol_rule is None and n_plies_dam_tol is None:
            raise Exception('Not enough input arugments!')
    else:
        for ind_ply in range(ss.size - 1):
            if not is_diso(ss[ind_ply], ss[ind_ply + 1], delta_angle):
                return False
    return True


def calc_n_viola_diso(
        mother_n_viola_diso, mother_ss_top, child_ss, level, n_panels, n_plies,
        pdl, constraints, middle_ply=False, reflexion=False,
        mother_ss_bot=None, last_level=False, odd_level=False):
    """
    returns the number of violations of the disorientation constraint in each
    panel of a laminate design

    INPUTS

    - mother_n_viola_diso: number of violations of the disorientation rule in
    each panel of the incomplete laminate design
    - mother_ss: panel stacking sequences of the incomplete laminate design
    - child_ss: possible fibre orientation for the new ply
    - level: level in the beam search tree
    - n_panels: number of panels in the laminate structure
    - n_plies: number of plies in the thickest panel
    - pdl: matrix of ply drops
    - constraints: set of constraints
    - mother_ss_top: top stacking sequences of the incomplete laminate design
    - mother_ss_bot: bottom stacking sequence  of the incomplete laminate
    design
    - reflexion = True during design of innermost plies of symmetric laminates
    - middle_ply = True during design of symmetric laminates with odd number of
    plies
    - last_level = True for the last level of the beam search
    - odd_level = True for a level of the beam search corresponding to a ply in
    the top of the laminate
    """
    n_viola_diso = np.copy(mother_n_viola_diso)

    ### last ply in asymmetric laminate
    if not constraints.sym and last_level:
        for ind_panel in range(n_panels):
            if pdl[ind_panel, level - 1] >= 0:
                if not odd_level: # top ply added
                    if not is_diso(child_ss[0],
                                   mother_ss_bot[ind_panel][0],
                                   constraints.delta_angle):
                        n_viola_diso[ind_panel] += 1

                else:
                    if not is_diso(child_ss[0],
                                   mother_ss_top[ind_panel][-1],
                                   constraints.delta_angle):
                        n_viola_diso[ind_panel] += 1

            else:
                if not is_diso(mother_ss_top[ind_panel][-1],
                               mother_ss_bot[ind_panel][0],
                               constraints.delta_angle):
                    n_viola_diso[ind_panel] += 1

    ### top
    if constraints.sym or not odd_level:
        for ind_panel in range(n_panels):
            if pdl[ind_panel, level - 1] >= 0:

                if mother_ss_top[ind_panel].size: # previous ply in the group
                    if not is_diso(child_ss[0],
                                   mother_ss_top[ind_panel][-1],
                                   constraints.delta_angle):
                        if constraints.sym:
                            n_viola_diso[ind_panel] += 2
                        else:
                            n_viola_diso[ind_panel] += 1

            # if a midlle ply is designed simultaneously
            if child_ss.size == 2 and pdl[ind_panel, level - 1] < 0 \
            and pdl[ind_panel, level] >= 0:
                    if not is_diso(mother_ss_top[ind_panel][-1],
                                   child_ss[1],
                                   constraints.delta_angle):
                        if constraints.sym:
                            n_viola_diso[ind_panel] += 2
                        else:
                            n_viola_diso[ind_panel] += 1
    else: ### bottom
        for ind_panel in range(n_panels):
            if pdl[ind_panel, level - 1] >= 0:
                if mother_ss_bot[ind_panel].size: # previous ply in the group
                    if not is_diso(child_ss[0],
                                    mother_ss_bot[ind_panel][0],
                                    constraints.delta_angle):
                        n_viola_diso[ind_panel] += 1
    return n_viola_diso


def calc_penalty_diso_mp_0(ss, delta_angle):
    """
    returns the matrix of violations of the disorientation constraint by a
    multipanel structure

    - viola_diso[ind_panel, index_ply] = 1 when the plies of indices
    {index_ply, index_ply + 1} violates the constraint, 0 otherwise

    INPUTS

    - ss: stacking sequences of each panel of the structure (list or arrays)
    - delta_angle: maximum variation of the fibre angle for adjacent plies
    """
    viola_diso = [[]]*len(ss)
    for ind_panel in range(len(ss)):
        viola_diso[ind_panel] = np.zeros((
            ss[ind_panel].size -1,), dtype=bool)
        for index_ply in range(ss[ind_panel].size - 1):
            viola_diso[ind_panel][index_ply] = not is_diso(
                ss[ind_panel][index_ply],
                ss[ind_panel][index_ply + 1],
                delta_angle)
    return viola_diso


def calc_penalty_diso_mp_1(ss, constraints):
    """
    returns the number of disorientation constraint violations per panel in a
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
    viola_diso = calc_penalty_diso_mp_0(ss, constraints.delta_angle)
    return np.array([np.sum(el) for el in viola_diso]) # sum per panel

def calc_n_penalty_diso_ss(stack, constraints):
    """
    returns the number of disorientation constraint violations in a stacking
    sequence

    INPUTS

    - stack: stacking sequence
    - constraints: set of constraints
    """
    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule'):
            if constraints.dam_tol_rule in {2, 3}:
                stack = np.copy(stack)[1:-1]
        else:
            if constraints.n_plies_dam_tol == 2:
                stack = np.copy(stack)[1:-1]
    return np.sum(calc_penalty_diso_mp_0([stack], constraints.delta_angle))

if __name__ == "__main__":
    print('*** Test for the function is_diso ***\n')
    print('Inputs:\n')
    angle1 = 90
    angle2 = np.array((90))
    print('outputs:\n')
    print(is_diso(angle1, angle2, delta_angle=45))

    print('*** Test for the function is_diso_ss ***\n')
    print('Inputs:\n')
    ss = np.array([45, 0, -45, 0, 45, 90, 45, 45, 0, -45, 0, 45])
    print('outputs:\n')
    print(is_diso_ss(ss, delta_angle=45, dam_tol=True, dam_tol_rule=2))

    print('*** Test for the function calc_n_viola_diso ***\n')
    print('Inputs:\n')
    mother_n_viola_diso = np.array([0, 0])
    mother_ss = [np.array([0, 0]), np.array([90])]
    child_ss = np.array([[0], [45], [90], [-45]])
    level = 3
    n_panels = 2
    n_plies = 4
    pdl = np.array([[0, 1, 2, 3], [0, -1, 2, 3]])
    ss_before = [np.array(()), np.array(())]
    ss_after = [np.array(()), np.array(())]
    bottom = False
    constraints = Constraints(
        sym=True, dam_tol=True, dam_tol_rule=2, delta_angle=45)
    print(f'mother_ss: {mother_ss}')
    print(f'child_ss: {child_ss}')
    print(f'mother_n_viola_diso: {mother_n_viola_diso}\n')
    print('outputs:\n')
    viola_diso = calc_n_viola_diso(
        mother_n_viola_diso, mother_ss, child_ss, level, n_panels, n_plies, pdl,
        constraints)
    print(f'viola_diso: {viola_diso}\n')

    print('*** Test for the function calc_penalty_diso_mp_0 ***\n')
    print('Inputs:\n')
    ss = [np.array([0, 90, 90]), np.array([0, 0, 0, 45, -45, 0])]
    delta_angle = 45
    print(f'ss = {ss}, delta_angle = {delta_angle}\n')
    print('outputs:\n')
    print(calc_penalty_diso_mp_0(ss, delta_angle))

    print('*** Test for the function calc_penalty_diso_mp_1 ***\n')
    print('Inputs:\n')
    constraints = Constraints(sym=True)
    constraints.dam_tol = False
    constraints.delta_angle = 45
    ss = [np.array([-45, 45, 90, 0, 0, 45, -45]),
          np.array([-45, 45, 90, -45, 0, 0, 45, -45]),
          np.array([-45, 45, 45, -45])]
    print(f'ss = {ss}')
    print('outputs:\n')
    print(calc_penalty_diso_mp_1(ss, constraints))

    print('*** Test for the function calc_n_penalty_diso_ss ***\n')
    constraints = Constraints(sym=True)
    constraints.dam_tol = False
    constraints.delta_angle = 45
    print('Inputs:\n')
    ss = np.array([0, 90, 90, 0, 90, 90])
    print(f'ss: {ss}')
    print(f'delta_angle: {constraints.delta_angle}')
    print('outputs:\n')
    print(f'Result: {calc_n_penalty_diso_ss(ss, constraints)}\n')
