# -*- coding: utf-8 -*-
"""
repair for disorientation and contiguity

- repair_diso_contig_from_outside_sym
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for symmetric laminates with swaps performed inwardly
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


def repair_diso_contig_from_outside_sym(
        ss, ply_queue, constraints, n_D1):
    """
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for symmetric laminates with swaps performed inwardly
    in the laminates

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    - n_D1: number of plies in the last permutation
    """
#    print(constraints)
#    print('n_D1', n_D1)
#    print_ss(ss)
#    print('ply_queue', ply_queue)

    stack, ind_1, ind_2 = initialise_repair_diso_contig(
        ss, constraints, from_inside=False, n_D1=n_D1)
#    print('ind_1', ind_1)
#    print('ind_2', ind_2)

    for ind_ply_1 in ind_1:
#        print('ind_ply_1', ind_ply_1)
        angle_found = False

        for ind_ply_2 in ind_2:
#            print('    ind_ply_2', ind_ply_2)

            if ss[ind_ply_2] != 666:
                stack_test = np.hstack((stack, ss[ind_ply_2]))
#                print('stack_test', stack_test)
                if constraints.diso and not is_diso(
                        stack_test[-1], stack_test[-2],
                        constraints.delta_angle):
                    continue

                if constraints.dam_tol:
                    if hasattr(constraints, 'dam_tol_rule') \
                    and constraints.dam_tol_rule in {2, 3}:
                        mini = max(1,
                                   stack_test.size - constraints.n_contig - 1)
                    if not hasattr(constraints, 'dam_tol_rule') and \
                    constraints.n_plies_dam_tol == 2:
                        mini = max(1,
                                   stack_test.size - constraints.n_contig - 1)
                    else:
                        mini = max(0,
                                   stack_test.size - constraints.n_contig - 1)
                else:
                    mini = max(0, stack_test.size - constraints.n_contig - 1)

                if constraints.contig and not is_contig(
                        stack_test[mini:], constraints.n_contig):
                    continue
                stack = stack_test
                ind_2 = np.delete(ind_2, np.argwhere(ind_2 == ind_ply_2))
                angle_found = True
#                print_ss(stack, 13)
                break

            else: # let's pick a ply from the queue
                angles_to_test = order_plies_to_test(
                    [stack[-1]], ply_queue, constraints)
                for angle in angles_to_test:
                    stack_test = np.hstack((stack, angle))
#                    print('stack_test', stack_test)
                    if constraints.diso and not is_diso(
                            stack_test[-1], stack_test[-2],
                            constraints.delta_angle):
                        continue

                    if constraints.dam_tol:
                        if hasattr(constraints, 'dam_tol_rule') \
                        and constraints.dam_tol_rule in {2, 3}:
                            mini = max(1,
                                       stack_test.size - constraints.n_contig - 1)
                        if not hasattr(constraints, 'dam_tol_rule') and \
                        constraints.n_plies_dam_tol == 2:
                            mini = max(1,
                                       stack_test.size - constraints.n_contig - 1)
                        else:
                            mini = max(0,
                                       stack_test.size - constraints.n_contig - 1)
                    else:
                        mini = max(0, stack_test.size - constraints.n_contig - 1)

                    if constraints.contig and not is_contig(
                            stack_test[mini:], constraints.n_contig):
                        continue
                    stack = stack_test
                    ind_2 = np.delete(
                        ind_2, np.argwhere(ind_2 == ind_ply_2)[0])
                    ply_queue.remove(angle)
#                    print_ss(stack, 13)
                    angle_found = True
                    break

                if angle_found:
                    break

        if not angle_found:
            # return semi-repaired stacking sequence
            for ind in ind_2:
                if ss[ind] != 666:
                    stack = np.hstack((stack, ss[ind]))
                else:
                    stack = np.hstack((stack, ply_queue.pop(0)))
            if ss.size % 2:
                stack = np.hstack((stack, ss[ss.size // 2], np.flip(stack)))
            else:
                stack = np.hstack((stack, np.flip(stack)))
            return stack, False

    if stack.size < ss.size // 2 - n_D1:
        # return semi-repaired stacking sequence
        raise Exception('This should not happen anymore')

#    print_ss(stack)

    # plies at the centre all designed simultaneously
    good_permut = []
    real_n_D1 = ind_2.size
    remaining_plies = np.empty((real_n_D1,), int)
    for counter, ind in enumerate(ind_2):
        if ss[ind] != 666:
            remaining_plies[counter] = ss[ind]
        else:
            remaining_plies[counter] = ply_queue.pop(0)
    print('remaining_plies', remaining_plies)

    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule') \
        and constraints.dam_tol_rule in {2, 3}:
            mini = max(1, stack.size - real_n_D1 - 1)
        if not hasattr(constraints, 'dam_tol_rule') and \
        constraints.n_plies_dam_tol == 2:
            mini = max(1, stack.size - real_n_D1 - 1)
        else:
            mini = max(0, stack.size - real_n_D1 - 1)
    else:
        mini = max(0, stack.size - real_n_D1 - 1)

    for elem in itertools.permutations(
        range(real_n_D1)):
        ss_group = remaining_plies[np.array(elem, int)]

        stack_test = np.hstack((stack[mini:], ss_group))
        if ss.size % 2:
            stack_test = np.hstack((
                stack_test, ss[ss.size // 2],
                np.flip(stack_test[-real_n_D1:])))
        else:
            stack_test = np.hstack((
                stack_test,
                np.flip(stack_test[-real_n_D1:])))
        if constraints.diso and not is_diso_ss(
                stack_test, constraints.delta_angle, dam_tol=False):
            continue
        if constraints.contig and not is_contig(
                stack_test, constraints.n_contig_c):
            continue
        good_permut.append(np.array(elem, int))
    print('good_permut', good_permut)

    good_permut = smallest_row(good_permut)
    if good_permut is None:
        # return semi-repaired stacking sequence
        stack = np.hstack((stack, remaining_plies))
        if ss.size % 2:
            stack = np.hstack((stack, ss[ss.size // 2]))
        stack = np.hstack((stack, np.flip(stack)))
        if ss.size % 2:
            stack = np.delete(stack, np.s_[ss.size // 2])
        return stack, False

    stack = np.hstack((stack, remaining_plies[good_permut]))
    if ss.size % 2:
        stack = np.hstack((stack, ss[ss.size // 2], np.flip(stack)))
    else:
        stack = np.hstack((stack, np.flip(stack)))
    if stack.size != ss.size:
        raise Exception('This should not happen')
    return stack, True

if __name__ == "__main__":
    print('\n*** Test for repair_diso_contig_from_outside_sym ***')
    constraints = Constraints(
        sym=True,
        bal=False,
        ipo=False,
        dam_tol=False,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=5,
        set_of_angles=[0, 45, -45, 90, 30, -30, 60, -60, 75, -75, 15, -15])

    ss = np.array([0, 45, -45, 666, 666, 666, 666], int)
    ss = np.hstack((ss, np.flip(ss)))
    ply_queue = [-45, 90, 90, 90]

    print_ss(ss, 40)
    n_D1 = 3
    print('ss.size', ss.size)
    ss, completed = repair_diso_contig_from_outside_sym(
        ss, ply_queue, constraints, n_D1)
    print('Repair successful?', completed)
    print_ss(ss, 40)

    if not is_diso_ss(ss, constraints.delta_angle, dam_tol=False):
        raise Exception('Disorientation rule not satisfied')

    if not is_contig(ss, constraints.n_contig):
        raise Exception('Contiguity rule not satisfied')
