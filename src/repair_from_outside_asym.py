# -*- coding: utf-8 -*-
"""
repair for disorientation and contiguity

- repair_diso_contig_from_outside_asym
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for asymmetric laminates with swaps performed inwardly
    in the laminates
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import itertools
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.disorientation import is_diso, is_diso_ss
from src.contiguity import is_contig
from src.pretty_print import print_ss
from src.repair_diso_contig_tools import initialise_repair_diso_contig
from src.repair_diso_contig_tools import smallest_row
from src.repair_diso_contig_tools import order_plies_to_test
from src.repair_diso_contig_tools import modify_indices

def repair_diso_contig_from_outside_asym(
        ss, ply_queue, constraints, n_D1):
    """
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for asymmetric laminates with swaps performed inwardly
    in the laminates

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    - n_D1: number of plies in the last permutation
    """
    stack_top, stack_bottom, ind_1, ind_2 = initialise_repair_diso_contig(
        ss, constraints, from_inside=False,
        n_D1=n_D1)
#    print('ind_1', ind_1)
#    print('ind_2', ind_2)

    ind_step = -1
    for ind_step, ind_ply_1 in enumerate(ind_1):
        ind_2_mod = modify_indices(ind_2, ind_ply_1, bottom=ind_step % 2 == 0)
#        print('ind_ply_1', ind_ply_1)
#        print('ind_2_mod', ind_2_mod)
        angle_found = False

        for ind_ply_2 in ind_2_mod:
#            print('    ind_ply_2', ind_ply_2)

            if ss[ind_ply_2] != 666:

                if ind_step % 2 == 0:
                    stack_test = np.hstack((ss[ind_ply_2], stack_bottom))
#                    print('stack_test', stack_test)
                    if constraints.diso and not is_diso(
                            stack_test[0], stack_test[1],
                            constraints.delta_angle):
                        continue
                    if constraints.contig and not is_contig(
                            stack_test[:2*constraints.n_contig],
                            constraints.n_contig):
                        continue
                    stack_bottom = stack_test
#                    print('stack_bottom', end='')
#                    print_ss(stack_bottom)
                else:
                    stack_test = np.hstack((stack_top, ss[ind_ply_2]))
#                    print('stack_test', stack_test)
                    if constraints.diso and not is_diso(
                            stack_test[-1], stack_test[-2],
                            constraints.delta_angle):
                        continue
                    if constraints.contig and not is_contig(
                            stack_test[-2*constraints.n_contig:],
                            constraints.n_contig):
                        continue
                    stack_top = stack_test
#                    print('stack_top', end='')
#                    print_ss(stack_top)
                ind_2 = np.delete(ind_2, np.argwhere(ind_2 == ind_ply_2))
                angle_found = True
                break

            else: # let's pick a ply from the queue

                if ind_step % 2 == 0:
                    angles_to_test = order_plies_to_test(
                        [stack_bottom[0]], ply_queue, constraints)
                    for angle in angles_to_test:
                        stack_test = np.hstack((angle, stack_bottom))
#                        print('stack_test', stack_test)
                        if constraints.diso and not is_diso(
                                stack_test[0], stack_test[1],
                                constraints.delta_angle):
                            continue
                        if constraints.contig and not is_contig(
                                stack_test[:2*constraints.n_contig],
                                constraints.n_contig):
                            continue
                        stack_bottom = stack_test
#                        print('stack_bottom', end='')
#                        print_ss(stack_bottom)
                        angle_found = True
                        ind_2 = np.delete(
                            ind_2, np.argwhere(ind_2 == ind_ply_2)[0])
                        ply_queue.remove(angle)
                        break
                else:
                    angles_to_test = order_plies_to_test(
                        [stack_top[-1]], ply_queue, constraints)
                    for angle in angles_to_test:
                        stack_test = np.hstack((stack_top, angle))
#                        print('stack_test', stack_test)
                        if constraints.diso and not is_diso(
                                stack_test[-1], stack_test[-2],
                                constraints.delta_angle):
                            continue
                        if constraints.contig and not is_contig(
                                stack_test[-2*constraints.n_contig:],
                                constraints.n_contig):
                            continue
                        stack_top = stack_test
#                        print('stack_top', end='')
#                        print_ss(stack_top)
                        angle_found = True
                        ind_2 = np.delete(
                            ind_2, np.argwhere(ind_2 == ind_ply_2)[0])
                        ply_queue.remove(angle)
                        break

                if angle_found:
                    break

        if not angle_found:
            # return semi-repaired stacking sequence
            for ind in ind_2:
                if ss[ind] != 666:
                    stack_top = np.hstack((stack_top, ss[ind]))
                else:
                    stack_top = np.hstack((stack_top, ply_queue.pop(0)))
            stack = np.hstack((stack_top, stack_bottom))
            return stack, False

    if stack_bottom.size + stack_top.size \
    < ss.size - n_D1:
        # return semi-repaired stacking sequence
        raise Exception('This should not happen anymore')

#    print('stack_bottom', end='')
#    print_ss(stack_bottom)
#    print('stack_top', end='')
#    print_ss(stack_top)

    # plies at the centre all designed simultaneously
    good_permut = []
    real_n_D1 = ind_2.size
    remaining_plies = np.empty((real_n_D1,), int)
    for counter, ind in enumerate(ind_2):
        if ss[ind] != 666:
            remaining_plies[counter] = ss[ind]
        else:
            remaining_plies[counter] = ply_queue.pop(0)
#    print('remaining_plies', remaining_plies)


    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule') \
        and constraints.dam_tol_rule in {2, 3}:
            mini_bottom = max(1, stack_bottom.size - real_n_D1 - 1)
            mini_top = min(stack_bottom.size - 1, real_n_D1 + 1)
        if not hasattr(constraints, 'dam_tol_rule') and \
        constraints.n_plies_dam_tol == 2:
            mini_bottom = max(1, stack_bottom.size - real_n_D1 - 1)
            mini_top = min(stack_bottom.size - 1, real_n_D1 + 1)
        else:
            mini_bottom = max(0, stack_bottom.size - real_n_D1 - 1)
            mini_top = min(stack_bottom.size, real_n_D1 + 1)
    else:
        mini_bottom = max(0, stack_bottom.size - real_n_D1 - 1)
        mini_top = min(stack_bottom.size, real_n_D1 + 1)

    ind_2 = np.sort(ind_2)
    for elem in itertools.permutations(
        range(real_n_D1)):
#        print('elem', elem)
        ss_group = remaining_plies[np.array(elem, int)]
        stack_test = np.hstack((
            stack_top[mini_bottom:], ss_group, stack_bottom[:mini_top]))
#        print('stack_test', stack_test)
        if constraints.diso and not is_diso_ss(
                stack_test, constraints.delta_angle, dam_tol=False):
            continue
        if constraints.contig and not is_contig(
                stack_test, constraints.n_contig_c):
            continue
        good_permut.append(np.array(elem, int))
#        print('stack_test', stack_test, elem)
#    print('good_permut', good_permut)

    good_permut = smallest_row(good_permut)
#    print('good_permut', good_permut)

    if good_permut is None:
        # return semi-repaired stacking sequence
        stack = np.hstack((stack_top, remaining_plies, stack_bottom))
        return stack, False

#    remaining_plies = remaining_plies[good_permut]
    stack = np.hstack((stack_top, remaining_plies[good_permut], stack_bottom))

    if stack.size != ss.size:
        raise Exception('This should not happen')
    return stack, True

if __name__ == "__main__":
    print('\n*** Test for repair_diso_contig_from_outside_asym ***')
    ss_ini = np.array([-45, 0, 666, 666, 666, 90, 45, 90], int)
    ply_queue = [-45, 0, 45]
    constraints = Constraints(
        sym=False,
        dam_tol=False,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=4,
        set_of_angles=[0, 45, -45, 90, 30, -30, 60, -60, 15, -15, 75, -75])
    print_ss(ss_ini, 200)
    n_D1 = 6
    print('ss_ini.size', ss_ini.size)
    print('ply_queue', ply_queue)
    ss, completed = repair_diso_contig_from_outside_asym(
        ss_ini, ply_queue, constraints, n_D1)
    print('Repair successful?', completed)
    print_ss(ss, 200)
