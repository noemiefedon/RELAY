# -*- coding: utf-8 -*-
"""
repair usefull during repair for disorientation and contiguity

- smallest_row
    returns the list with the smallest successive values in a list of lists

- initialise_repair_diso_contig
    prepares the reapir for disorientation and/or contiguity
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_ss
from src.repair_tools import calc_delta_angle_2

def smallest_row(permut):
    """
    returns the list with the smallest successive values in a list of lists
    """
    if not permut:
        return None # return unrepaired stacking sequence

    permut = np.array(permut)
#    print('permut', permut)

    for ind in range(permut.shape[1]):
#        print('ind', ind)
        permut = permut[permut[:,ind].argsort(),]
#        print('permut', permut)

        n_to_keep = 1
        for ind2 in range(1, permut.shape[0]):
            if permut[ind2, ind] == permut[0, ind]:
                n_to_keep += 1
            else:
                break
#        print('n_to_keep', n_to_keep)

        permut = permut[:n_to_keep]
#        print('permut', permut)

        if permut.shape[0] == 1:
            break

    return permut[0]


def initialise_repair_diso_contig(
        ss, constraints, from_inside, n_D1):
    """
    prepares the repair for disorientation and/or contiguity

    INPUTS

    ss: stacking sequence (array)
    constraints: design and manufacturing constraints
    from_inside = True, from inwards repair
    n_D1: number of plies in the last permutation

    OUTPUTS

    ind_1: indices of plies to be determined in the repaired laminate
    ind_2: indices of plies to be selected in the initial laminate
    """
    if constraints.sym and not from_inside:

        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                # two outer plies unchanged for sure
                stack = np.array([ss[0], ss[1]], dtype='int16')
                ind_1 = range(1, ss.size // 2 - n_D1 - 1)
                ind_2 = range(2, ss.size // 2)
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                # two outer plies unchanged for sure
                stack = np.array([ss[0], ss[1]], dtype='int16')
                ind_1 = range(1, ss.size // 2 - n_D1 - 1)
                ind_2 = range(2, ss.size // 2)
            else:
                # only 1 outer ply unchanged for sure
                stack = np.array([ss[0]], dtype='int16')
                ind_1 = range(0, ss.size // 2 - n_D1 - 1)
                ind_2 = range(1, ss.size // 2)
        else:
            # only 1 outer ply unchanged for sure
            stack = np.array([ss[0]], dtype='int16')
            ind_1 = range(0, ss.size // 2 - n_D1 - 1)
            ind_2 = range(1, ss.size // 2)

        return stack, np.array(ind_1), np.array(ind_2)

    if not from_inside: # and not constraints.sym
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                # two outer plies unchanged for sure
                stack_top = np.array([ss[0], ss[1]], dtype='int16')
                stack_bottom = np.array([ss[-2], ss[-1]], dtype='int16')
                ind_beg = 1
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                # two outer plies unchanged for sure
                stack_top = np.array([ss[0], ss[1]], dtype='int16')
                stack_bottom = np.array([ss[-2], ss[-1]], dtype='int16')
                ind_beg = 1
            else:
                # only 1 outer ply unchanged for sure
                stack_top = np.array([ss[0]], dtype='int16')
                stack_bottom = np.array([ss[-1]], dtype='int16')
                ind_beg = 0
        else:
            # only 1 outer ply unchanged for sure
            stack_top = np.array([ss[0]], dtype='int16')
            stack_bottom = np.array([ss[-1]], dtype='int16')
            ind_beg = 0

        indices = np.arange(ss.size)
        beginning = np.copy(indices[ind_beg:indices.size // 2])
        ending = np.copy(indices[indices.size // 2:ss.size - ind_beg][::-1])
#        print('indices', indices)
#        print('beginning', beginning)
#        print('ending', ending)
#        print('indices[::2]', indices[::2])
#        print('indices[1::2]', indices[1::2])
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                indices = np.zeros((ss.size - 2,), int)
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                indices = np.zeros((ss.size - 2,), int)
            else:
                indices = np.zeros((ss.size,), int)
        else:
            indices = np.zeros((ss.size,), int)
        indices[::2] = ending
        indices[1::2] = beginning
#        print('indices', indices)

        ind_2 = indices[2:]
        if indices.size - n_D1 - 2 > 0:
            ind_1 = indices[:indices.size - n_D1 - 2]
        else:
            ind_1 = np.array((), int)

        return stack_top, stack_bottom, np.array(ind_1), np.array(ind_2)

    if constraints.sym: # and from_inside

        if ss.size % 2 == 0:
            stack = np.array([
                ss[ss.size // 2],
                ss[ss.size // 2]], dtype='int16')
        else:
            stack = np.array([ss[ss.size // 2]], dtype='int16')

        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                # two outer plies unchanged for sure
                stack_dam_tol = np.array([ss[0], ss[1]], dtype='int16')
                ind_1 = range(2 + n_D1,
                              ss.size // 2 - 1 + ss.size % 2)[::-1]
                ind_2 = range(2, ss.size // 2 - 1 + ss.size % 2)[::-1]
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                # two outer plies unchanged for sure
                stack_dam_tol = np.array([ss[0], ss[1]], dtype='int16')
                ind_1 = range(2 + n_D1,
                              ss.size // 2 - 1 + ss.size % 2)[::-1]
                ind_2 = range(2, ss.size // 2 - 1 + ss.size % 2)[::-1]
            else:
                # one outer ply unchanged for sure
                stack_dam_tol = np.array([ss[0]], dtype='int16')
                ind_1 = range(1 + n_D1,
                              ss.size // 2 - 1 + ss.size % 2)[::-1]
                ind_2 = range(1, ss.size // 2 - 1 + ss.size % 2)[::-1]
        else:
            # no outer plies unchanged for sure
            stack_dam_tol = np.array((), dtype='int16')
            ind_1 = range(0 + n_D1,
                          ss.size // 2 - 1 + ss.size % 2)[::-1]
            ind_2 = range(0, ss.size // 2 - 1 + ss.size % 2)[::-1]

        return stack_dam_tol, stack, np.array(ind_1), np.array(ind_2)

    # if not constraints.sym and from_inside

    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule') \
        and constraints.dam_tol_rule in {2, 3}:
            # two outer plies unchanged for sure
            stack_top = np.array([ss[0], ss[1]], dtype='int16')
            stack_bottom = np.array([ss[-2], ss[-1]], dtype='int16')
            indices = np.arange(2, ss.size - 2)
        elif not hasattr(constraints, 'dam_tol_rule') and \
        constraints.n_plies_dam_tol == 2:
            # two outer plies unchanged for sure
            stack_top = np.array([ss[0], ss[1]], dtype='int16')
            stack_bottom = np.array([ss[-2], ss[-1]], dtype='int16')
            indices = np.arange(2, ss.size - 2)
        else:
            # one outer ply unchanged for sure
            stack_top = np.array([ss[0]], dtype='int16')
            stack_bottom = np.array([ss[-1]], dtype='int16')
            indices = np.arange(1, ss.size - 1)
    else:
        # no outer plies unchanged for sure
        stack_top = np.array((), dtype='int16')
        stack_bottom = np.array((), dtype='int16')
        indices = np.arange(0, ss.size)

    beginning = np.copy(indices[0:indices.size // 2][::-1])
    ending = np.copy(indices[indices.size // 2:ss.size])
    stack = np.array([ss[ending[0]]], int)
#    print('indices', indices)
#    print('beginning', beginning)
#    print('ending', ending)
#    print('indices[::2]', indices[::2])
#    print('indices[1::2]', indices[1::2])
    indices = np.zeros((beginning.size + ending.size,), int)
    indices[::2] = ending
    indices[1::2] = beginning
#    print('indices', indices)

    ind_2 = indices[1:]
    if indices.size - n_D1 - 1 > 0:
        ind_1 = indices[:indices.size - n_D1 - 1]
    else:
        ind_1 = np.array((), int)

    return stack_top, stack, stack_bottom, np.array(ind_1), np.array(ind_2)


def modify_indices(ind_2, ind_ply_1, bottom):
    """
    modifies the indices of plies to be selected in the initial laminate during
    the inwards repair for disorientation and contiguity for asymmetric
    laminates to keep the laminate unchanged if possible

    INPUTS
    """
    ind_2_mod = np.copy(ind_2)
#    print('@@@@')
#    print('ind_ply_1', ind_ply_1)
#    print('ind_2', ind_2)
#    print('bottom', bottom)
#    if (bottom and ind_2[1] == ind_ply_1 - 1) \
#    or (not bottom and ind_2[1] == ind_ply_1 + 1):
#        ind_2_mod[0], ind_2_mod[1] = ind_2_mod[1], ind_2_mod[0]
    return ind_2_mod

def order_plies_to_test(adjacent_plies, queue, constraints):
    """
    list the plies from the queue in order of preference for them to be
    satcked next to the plies in the adjacent_plies list

    preference for:
        - plies at the same fibre orientations than the adjacent plies, to
        simplify the reapir for disorientation and/or contiguity later
        - plies with the least difference in fibre orientation.
        - positive diffrence in fibre orientation chosen if ambiguity
    """
    if len(adjacent_plies) == 0:
        return np.unique(queue)

    if len(adjacent_plies) == 1:
        ordered_plies = []
        angle1 = adjacent_plies[0]
        diff_angle = np.zeros((constraints.n_set_of_angles,), int)

        for ind_angle, angle2 in enumerate(constraints.set_of_angles):
            diff_angle[ind_angle] = calc_delta_angle_2(angle1, angle2)
#            print('angle2', angle2, 'diff_angle', diff_angle[ind_angle])
        abs_diff_angle = abs(diff_angle)

#        print('diff_angle', diff_angle)
#        print('abs_diff_angle', abs_diff_angle)

        while (abs_diff_angle - 666 != 0).any():
            indices = np.argwhere(abs_diff_angle == min(abs_diff_angle))
#            print('abs_diff_angle', abs_diff_angle)
#            print('indices', indices)
            if indices.size == 2:
                if diff_angle[indices[0]] > 0:
                    ordered_plies.append(
                        constraints.set_of_angles[indices[0][0]])
                    ordered_plies.append(
                        constraints.set_of_angles[indices[1][0]])
                else:
                    ordered_plies.append(
                        constraints.set_of_angles[indices[1][0]])
                    ordered_plies.append(
                        constraints.set_of_angles[indices[0][0]])
                abs_diff_angle[indices[0]] = 666
                abs_diff_angle[indices[1]] = 666
            else:
                ordered_plies.append(
                    constraints.set_of_angles[indices[0][0]])
                abs_diff_angle[indices[0]] = 666
#            print('ordered_plies', ordered_plies)
#        print('ordered_plies', ordered_plies)

        return [angle for angle in ordered_plies if angle in queue]


    if len(adjacent_plies) == 2:
        raise Exception('Case not considered yet')


if __name__ == "__main__":
    print('\n*** Test for the function order_plies_to_test ***')
#    constraints = Constraints(
#        sym=False,
#        bal=True,
#        dam_tol=True,
#        diso=True,
#        contig=True,
#        delta_angle=45,
#        n_contig=3,
#        set_of_angles=[0, 45, -45, 30, -30, 60, -60, 90])
#    adjacent_plies = [0]
#    queue = [-30, 0, 60, -45, 60, -30, 60]
#    print(order_plies_to_test(adjacent_plies, queue, constraints))
#
#    print('\n*** Test for the function smallest_row ***')
#    permut = [[0, 2, 1], [2, 1, 0], [1, 2, 0], [0, 1, 2]]
#    print('best_permut', smallest_row(permut))
#    print('expectecd best_permut', [0, 1, 2])


    print('\n*** Test for the function initialise_repair_diso_contig ***')
    ss = np.arange(10)
    ss = np.hstack((ss, np.flip(ss)))
    print_ss(ss)
    constraints = Constraints(
        sym=False,
        dam_tol=True,
        dam_tol_rule=3,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=3)
    from_inside = True
    n_D1 = 3

    if constraints.sym and not from_inside:
        stack, ind_1, ind_2 = initialise_repair_diso_contig(
            ss, constraints, from_inside, n_D1)
        print('stack', end='')
        print_ss(stack, 40)
    elif not constraints.sym and not from_inside:
        stack_top, stack_bottom, ind_1, ind_2 \
        = initialise_repair_diso_contig(
            ss, constraints, from_inside, n_D1)
        print('stack_top', end='')
        print_ss(stack_top, 40)
        print('stack_bottom', end='')
        print_ss(stack_bottom, 40)
    elif constraints.sym and from_inside:
        stack_dam_tol, stack, ind_1, ind_2 \
        = initialise_repair_diso_contig(
            ss, constraints, from_inside, n_D1)
        print('stack_dam_tol', end='')
        print_ss( stack_dam_tol, 40)
        print('stack', end='')
        print_ss(stack, 40)
    elif not constraints.sym and from_inside:
        stack_top, stack, stack_bottom, ind_1, ind_2 \
        = initialise_repair_diso_contig(
            ss, constraints, from_inside, n_D1)
        print('stack_top', end='')
        print_ss(stack_top, 40)
        print('stack_bottom', end='')
        print_ss(stack_bottom, 40)
        print('stack', end='')
        print_ss(stack, 40)

    print('ind_1', list(ind_1))
    print('ind_2', list(ind_2))
