# -*- coding: utf-8 -*-
"""
repair for disorientation and contiguity

- repair_diso_contig_from_inside_sym
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for symmetric laminates with swaps from the inside out
    of the laminates
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

def repair_diso_contig_from_inside_sym(
        ss, ply_queue, constraints, n_D1):
    """
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for symmetric laminates with swaps from the inside out
    of the laminates

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    - n_D1: number of plies in the last permutation
    """
    stack_dam_tol, stack, ind_1, ind_2 = \
    initialise_repair_diso_contig(
        ss, constraints, from_inside=True,
        n_D1=n_D1)
#    print('ind_1', ind_1)
#    print('ind_2', ind_2)

    for ind_ply_1 in ind_1:
#        print('ind_ply_1', ind_ply_1)
        angle_found = False

        for ind_ply_2 in ind_2:
#            print('    ind_ply_2', ind_ply_2)

            if ss[ind_ply_2] != 666:
                stack_test = np.hstack((ss[ind_ply_2], stack, ss[ind_ply_2]))
#                print('stack_test', stack_test)
                if constraints.diso and not is_diso(
                        stack_test[0], stack_test[1],
                        constraints.delta_angle):
                    continue
                if constraints.contig and not is_contig(
                        stack_test[:2*constraints.n_contig],
                        constraints.n_contig):
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
                    stack_test = np.hstack((angle, stack, angle))
#                    print('stack_test', stack_test)
                    if constraints.diso and not is_diso(
                            stack_test[0], stack_test[1],
                            constraints.delta_angle):
                        continue
                    if constraints.contig and not is_contig(
                            stack_test[:2*constraints.n_contig],
                            constraints.n_contig):
                        continue
                    stack = stack_test
                    ind_2 = np.delete(
                        ind_2, np.argwhere(ind_2 == ind_ply_2)[0])
                    ply_queue.remove(angle)
                    angle_found = True
                    print_ss(stack, 13)
                    break

                if angle_found:
                    break

        if not angle_found:
            # return semi-repaired stacking sequence
            for ind_ply_2 in ind_2:
                if ss[ind_ply_2] != 666:
                    angle = ss[ind_ply_2]
                else:
                    angle = ply_queue.pop(0)
                stack = np.hstack((angle, stack, angle))
            stack = np.hstack((stack_dam_tol, stack, np.flip(stack_dam_tol)))
            return stack, False

    if stack.size + 2*stack_dam_tol.size \
    < ss.size - 2*n_D1:
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
#    print('remaining_plies', remaining_plies)


    for elem in itertools.permutations(
        range(real_n_D1)):
        ss_group = remaining_plies[np.array(elem, int)]
        if constraints.dam_tol:
            stack_test = np.hstack((stack_dam_tol[-1], ss_group))
        else:
            stack_test = ss_group
        stack_test = np.hstack((stack_test, stack, np.flip(stack_test)))
        stack_test = stack_test[:2 + 2 * constraints.n_contig \
                                + real_n_D1]
#        print('stack_test', stack_test)
        if constraints.diso and not is_diso_ss(
                stack_test, constraints.delta_angle, dam_tol=False):
            continue
        if constraints.contig and not is_contig(
                stack_test, constraints.n_contig_c):
            continue
        good_permut.append(np.array(elem, int))
#    print('good_permut', good_permut)

    good_permut = smallest_row(good_permut)
    if good_permut is None:
        # return semi-repaired stacking sequence
        stack = np.hstack((np.flip(remaining_plies), stack, remaining_plies))
        stack = np.hstack((stack_dam_tol , stack, np.flip(stack_dam_tol)))
        return stack, False

    stack = np.hstack((
        stack_dam_tol,
        remaining_plies[good_permut],
        stack,
        np.flip(remaining_plies[good_permut]),
        np.flip(stack_dam_tol)))
    if stack.size != ss.size:
        raise Exception('This should not happen')
    return stack, True

if __name__ == "__main__":
    print('\n*** Test for repair_diso_contig_from_inside_sym ***')
    ss = np.array([0, 90, 45, 45, 666, 666, 666], int)
    ss = np.hstack((ss, 0, np.flip(ss)))
    constraints = Constraints(
        sym=True,
        dam_tol=True,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=3,
        set_of_angles=[0, 45, -45, 90])
    print_ss(ss, 40)
    n_D1 = 4
    ply_queue = [-45, -45, 0]
    print('ss.size', ss.size)
    ss, completed = repair_diso_contig_from_inside_sym(
        ss, ply_queue, constraints, n_D1)
    print('Repair successful?', completed)
    print_ss(ss, 40)
