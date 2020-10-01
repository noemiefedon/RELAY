# -*- coding: utf-8 -*-
"""
repair for disorientation and contiguity

- repair_diso_contig_from_inside_asym
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for asymmetric laminates with swaps from the inside out
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

def repair_diso_contig_from_inside_asym(
        ss, ply_queue, constraints, n_D1):
    """
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule for asymmetric laminates with swaps from the inside out
    of the laminates

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    - n_D1: number of plies in the last permutation
    """
    stack_top, stack, stack_bottom, ind_1, ind_2 = \
    initialise_repair_diso_contig(
        ss, constraints, True, n_D1)
#    print('ind_1', ind_1)
#    print('ind_2', ind_2)

    ind_step = -1
    for ind_step, ind_ply_1 in enumerate(ind_1):
#        print('ind_ply_1', ind_ply_1)

        for ind_ply_2 in ind_2:
#            print('    ind_ply_2', ind_ply_2)
            angle_found = False

            if ss[ind_ply_2] != 666:

                if ind_step % 2 == 0:
                    stack_test = np.hstack((ss[ind_ply_2], stack))
    #                print('stack_test', stack_test)
                    if constraints.diso and not is_diso(
                            stack_test[0], stack_test[1],
                            constraints.delta_angle):
                        continue
                    if constraints.contig and not is_contig(
                            stack_test[:2*constraints.n_contig],
                            constraints.n_contig):
                        continue
                else:
                    stack_test = np.hstack((stack, ss[ind_ply_2]))
    #                print('stack_test', stack_test)
                    if constraints.diso and not is_diso(
                            stack_test[-1], stack_test[-2],
                            constraints.delta_angle):
                        continue
                    if constraints.contig and not is_contig(
                            stack_test[-2*constraints.n_contig:],
                            constraints.n_contig):
                        continue
                stack = stack_test
                ind_2 = np.delete(ind_2, np.argwhere(ind_2 == ind_ply_2))
                angle_found = True
    #            print_ss(stack, 13)
                break

            else: # let's pick a ply from the queue

                if ind_step % 2 == 0:
                    angles_to_test = order_plies_to_test(
                    [stack[0]], ply_queue, constraints)
                    for angle in angles_to_test:
                        stack_test = np.hstack((angle, stack))
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
                        ind_2 = np.delete(
                            ind_2, np.argwhere(ind_2 == ind_ply_2)[0])
                        angle_found = True
                        ply_queue.remove(angle)
            #            print_ss(stack, 13)
                        break
                else:
                    angles_to_test = order_plies_to_test(
                        [stack[-1]], ply_queue, constraints)
                    for angle in angles_to_test:
                        stack_test = np.hstack((stack, angle))
        #                print('stack_test', stack_test)
                        if constraints.diso and not is_diso(
                                stack_test[-1], stack_test[-2],
                                constraints.delta_angle):
                            continue
                        if constraints.contig and not is_contig(
                                stack_test[-2*constraints.n_contig:],
                                constraints.n_contig):
                            continue
                        stack = stack_test
                        ind_2 = np.delete(
                            ind_2, np.argwhere(ind_2 == ind_ply_2)[0])
                        angle_found = True
                        ply_queue.remove(angle)
            #            print_ss(stack, 13)
                        break

                if angle_found:
                    break

        if not angle_found:
            # return semi-repaired stacking sequence
            for ind in ind_2:
                ind_step += 1
                if ind_step % 2 == 0:
                    if ss[ind] != 666:
                        stack = np.hstack((ss[ind], stack))
                    else:
                        stack = np.hstack((ply_queue.pop(0), stack))
                else:
                    if ss[ind] != 666:
                        stack = np.hstack((stack, ss[ind]))
                    else:
                        stack = np.hstack((stack, ply_queue.pop(0)))
            stack = np.hstack((stack_top, stack, stack_bottom))
            if stack.size != ss.size:
                raise Exception('This should not happen')
            return stack, False

    if stack_bottom.size + stack_top.size  + stack.size \
    < ss.size - n_D1:
        # return semi-repaired stacking sequence
        raise Exception('This should not happen anymore')

#    print('stack', stack)
#    print('stack_top', stack_top)
#    print('stack_bottom', stack_bottom)

    # plies at the centre all designed simultaneously
    ind_2_before = np.copy(ind_2)
    ply_queue_ini =ply_queue[:]
    ind_2 = np.sort(ind_2)
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
#        print('elem', elem)
        ss_group = remaining_plies[np.array(elem, int)]
        if constraints.dam_tol:
            stack_test = np.hstack((
                stack_top[-1],
                ss_group[:real_n_D1//2],
                stack,
                ss_group[real_n_D1//2:],
                stack_bottom[0]))
        else:
            stack_test = np.hstack((
                ss_group[:real_n_D1//2],
                stack,
                ss_group[real_n_D1//2:]))
#        print('stack_test', stack_test)
        stack_test_1 = stack_test[
            :3 + real_n_D1 + constraints.n_contig]
        if constraints.diso and not is_diso_ss(
                stack_test_1, constraints.delta_angle, dam_tol=False):
            continue
        if constraints.contig and not is_contig(
                stack_test_1, constraints.n_contig_c):
            continue
        number = stack_test.size - 4 - real_n_D1 \
        - constraints.n_contig
        if number > 0:
            stack_test_2 = stack_test[number:]
            if constraints.diso and not is_diso_ss(
                    stack_test_2, constraints.delta_angle, dam_tol=False):
                continue
            if constraints.contig and not is_contig(
                    stack_test_2, constraints.n_contig_c):
                continue
        good_permut.append(np.array(elem, int))
#    print('good_permut', good_permut)

    good_permut = smallest_row(good_permut)
    if good_permut is None:
        # return semi-repaired stacking sequence
        for ind in ind_2_before:
            ind_step += 1
            if ind_step % 2 == 0:
                if ss[ind] != 666:
                    stack = np.hstack((ss[ind], stack))
                else:
                    stack = np.hstack((ply_queue_ini.pop(0), stack))
            else:
                if ss[ind] != 666:
                    stack = np.hstack((stack, ss[ind]))
                else:
                    stack = np.hstack((stack, ply_queue_ini.pop(0)))
        stack = np.hstack((stack_top, stack, stack_bottom))
        if stack.size != ss.size:
            raise Exception('This should not happen')
        return stack, False

    ss_group = remaining_plies[good_permut]
    stack = np.hstack((
        stack_top,
        ss_group[:real_n_D1//2],
        stack,
        ss_group[real_n_D1//2:],
        stack_bottom))
    if stack.size != ss.size:
        raise Exception('This should not happen')
    return stack, True

if __name__ == "__main__":
    print('\n*** Test for repair_diso_contig_from_inside_asym ***')
    ss_ini = np.array([ 90, -45, -45,   0,   0, -45, -45,   0,
                   45,  90,  90,  90,  90, -45,   0,  45,
                   0, -45,  90, -45,  90,  45,  90,  90,  90,
                   90, -45, -45, -45,   0,   0, -45, -45,  90,
                   -45, -45,   0, -45, -45,
                   0,   0,   0,  45,  45,   0,  45,  45,   0,  45,   0], int)
    constraints = Constraints(
        sym=False,
        dam_tol=True,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=4,
        set_of_angles=[0, 45, -45, 90])
    print_ss(ss_ini, 200)
    n_D1 = 4
    ply_queue = []
    print('ss_ini.size', ss_ini.size)
    ss, completed = repair_diso_contig_from_inside_asym(
        ss_ini, ply_queue, constraints, n_D1)
    print('Repair successful?', completed)
    print_ss(ss - ss_ini, 200)