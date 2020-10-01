# -*- coding: utf-8 -*-
"""
repair for 10% rule and balance

- repair_10_bal
    repairs a laminate regarding the 10% rule and balance

- calc_mini_10:
    returns the minimum number of plies in the 0/90/+45/-45 directions to
    satisfy the 10% rule

- calc_current_10_2:
    returns the current number of plies in the 0/90/+45/-45 directions

- is_equal
    returns True if the set of partial stacking sequence + ply queue matches
    the initial stacking sequence
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import math as ma
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_ss
from src.repair_tools import RepairError


def repair_10_bal(ss_ini, mini_10, constraints):
    """
    repairs a laminate regarding the 10% rule and balance
    """

#    if not constraints.ipo and constraints.percent_45_135 != 0:
#        raise Exception("""
#Repair for 10% rule not implemented for laminates with no balance requirements
#and a limit percentage for the ply orientated in the combined +-45 direction!
#""")

#    print('initial')
#    print_ss(ss_ini)

    ss = np.copy(ss_ini)
    if not constraints.ipo and not constraints.rule_10_percent:
        return ss, []

    if constraints.sym and ss.size % 2 and ss[ss.size // 2] not in {0, 90}:
        ss[ss.size // 2] = 0

    if constraints.sym:
        ind_plies = np.array(range(0, ss.size // 2))
    else:
        ind_plies = np.arange(ss.size)
        beginning = np.copy(ind_plies[0:ind_plies.size // 2])
        ending = np.copy(ind_plies[ind_plies.size // 2:ss.size][::-1])
        ind_plies = np.zeros((ss.size,), int)
        ind_plies[::2] = ending
        ind_plies[1::2] = beginning
#    print('ind_plies', list(ind_plies))

    ply_queue = []
    if constraints.rule_10_percent:
        for elem in range(ma.ceil(mini_10[0])):
            ply_queue.append(0)
        for elem in range(ma.ceil(mini_10[1])):
            ply_queue.append(90)
        for elem in range(ma.ceil(mini_10[2])):
            ply_queue.append(45)
        for elem in range(ma.ceil(mini_10[3])):
            ply_queue.append(-45)
#    print('initial ply queue', ply_queue)

    if constraints.rule_10_percent and constraints.percent_45_135:
        missing_extra_45_135 = ma.ceil(mini_10[4]) \
        - ma.ceil(mini_10[2]) - ma.ceil(mini_10[3])
    else:
        missing_extra_45_135 = 0

    counter_remaining_plies = ind_plies.size

    if constraints.rule_10_percent:
        pass

    change = False
    for counter, ind_ply in enumerate(ind_plies):
#        print()
#        print('ind_ply', ind_ply)
#        print('new_angle', ss[ind_ply])
#        print('ply_queue', ply_queue)
#        print_ss(ss)

        ply_queue_before = ply_queue[:]
        new_angle = ss[ind_ply]
#        print('ind_ply', ind_ply, 'new_angle', new_angle)
        counter_remaining_plies -= 1

        if new_angle in ply_queue:
            ply_queue.remove(new_angle)
        else:
            if constraints.ipo and not new_angle in (0, 90):
                ply_queue.append(-new_angle)
            if not constraints.ipo and new_angle in (45, -45):
                missing_extra_45_135 = max(0, missing_extra_45_135 - 1)

#        print('ply_queue', ply_queue, len(ply_queue))
        if counter_remaining_plies < len(ply_queue) + missing_extra_45_135:
            change = True
            last_ply = counter
            ply_queue = ply_queue_before
            for ind_ply in ind_plies[counter:]:
                ss[ind_ply] = 666
                if constraints.sym:
                    ss[ss.size - ind_ply - 1] = 666
            break
#        print_ss(ss)
#        print('ply_queue', ply_queue)

#    print('last_ply', last_ply)
#    print('ply_queue', ply_queue, len(ply_queue))

#    print_ss(ss)

    for ind in range(missing_extra_45_135):
        if ind % 2:
            ply_queue.append(45)
        else:
            ply_queue.append(-45)

    if change and last_ply + len(ply_queue) != len(ind_plies):
#        ply_queue_2 = ply_queue[:]
#        ply_queue_2.append(90)
#        if (constraints.sym \
#            and np.isclose(np.sort(np.array(2*ply_queue_2)),
#                           np.sort(ss_ini[ss == 666])).all()) \
#        or (not constraints.sym \
#            and np.isclose(np.sort(np.array(ply_queue_2)),
#                           np.sort(ss_ini[ss == 666])).all()):
#            ply_queue.append(90)
#        else:
        ply_queue.append(0)

    return ss, ply_queue


def calc_mini_10(constraints, n_plies):
    """
    returns the minimum number of plies in the 0/90/+45/-45/+-45 directions to
    satisfy the 10% rule (array)

    INPUTS

        ss: stacking sequence (array)
        constraints: constraints (instance of the class Constraints)
    """
    mini_10 = np.zeros((5,), float)
    mini_10[0] = ma.ceil(constraints.percent_0 * n_plies)
    mini_10[1] = ma.ceil(constraints.percent_90 * n_plies)
    mini_10[2] = ma.ceil(constraints.percent_45 * n_plies)
    mini_10[3] = ma.ceil(constraints.percent_135 * n_plies)
    mini_10[4] = ma.ceil(constraints.percent_45_135 * n_plies)
    if constraints.ipo:
        mini_10[2] = max(mini_10[2], mini_10[3])
        if mini_10[4] % 2:
            mini_10[4] += 1
            mini_10[4] = max(mini_10[4], 2 * mini_10[2])
        mini_10[2] = max(mini_10[2], mini_10[4] // 2)
        mini_10[3] = mini_10[2]

    if constraints.sym:
        mini_10 /= 2
        # middle ply can only be oriented at 0 or 90 degrees
        if n_plies % 2:
            mini_10[2:] = np.ceil(mini_10[2:])
        else:
            mini_10 = np.ceil(mini_10)

        if constraints.ipo:
            if mini_10[4] % 2:
                mini_10[4] += 1
                mini_10[4] = max(mini_10[4], 2 * mini_10[2])
            mini_10[2] = max(mini_10[2], mini_10[4] // 2)
            mini_10[3] = mini_10[2]

    return mini_10

def calc_current_10_2(ss, sym):
    """
    returns the current number of plies in the 0/90/+45/-45 directions

    INPUTS

        ss: stacking sequence (array)
        sym: True for symmetric laminates (boolean)
    """
    current_10 = np.zeros((4,), float)
    if sym:
        lenn = ss.size // 2
        current_10[0] = sum(ss[:lenn] == 0)
        current_10[1] = sum(ss[:lenn] == 90)
        current_10[2] = sum(ss[:lenn] == 45)
        current_10[3] = sum(ss[:lenn] == -45)
        current_10[4] = current_10[2] + current_10[3]
        if ss.size % 2:
            if ss[lenn] == 0:
                current_10[0] += 1/2
            elif ss[lenn] == 90:
                current_10[1] += 1/2
            else:
                raise RepairError("""
This should not happen, plies at the midle surface at another fibre orientation
than 0 or 90 deg""")
    else:
        current_10[0] = sum(ss == 0)
        current_10[1] = sum(ss == 90)
        current_10[2] = sum(ss == 45)
        current_10[3] = sum(ss == -45)
        current_10[4] = current_10[2] + current_10[3]
    return current_10

def is_equal(ss, ply_queue, ss_ini, sym):
    """
    returns True if the set of partial stacking sequence + ply queue matches
    the initial stacking sequence
    """
#    print('@')
#    print_ss(ss_ini, 200)
#    print_ss(ss, 200)
#    print(ply_queue)

    if not np.isclose(ss[ss != 666], ss_ini[ss != 666] ).all():
        return False

#    print('np.sort(np.array(ply_queue))', np.sort(np.array(ply_queue)))
#    print('np.sort(ss_ini[ss == 666])', np.sort(ss_ini[ss == 666]))
#
    if sym:
        if not np.isclose(np.sort(np.array(2*ply_queue)),
                          np.sort(ss_ini[ss == 666])).all():
            return False
    else:
        if not np.isclose(np.sort(np.array(ply_queue)),
                          np.sort(ss_ini[ss == 666])).all():
            return False
    return True

if __name__ == "__main__":
    constraints = Constraints(
        sym=True,
        bal=False,
        dam_tol=False,
        rule_10_percent=True,
        ipo=True,
        n_contig=4,
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10,
        set_of_angles=[0, 45, -45, 90])


#    print('\n\n*** Test for the function is_equal ***')
#    ss_ini = np.array([
#        -45, 45, 60, 15, -15, 60, 30, 45, 0])
#    ss = np.array([
#        -45, 45, 60, 15, -15, 60, 30, 45, 0])
#    ply_queue = []
#    print(is_equal(ss, ply_queue, ss_ini, constraints.sym))

    print('\n\n*** Test for the function calc_mini_10 ***')
#    mini_10 = calc_mini_10(constraints, 40)
#    print('\nmini_10', mini_10)
#
#    n_45 = ma.ceil(mini_10[2])
#    n_135 = ma.ceil(mini_10[3])
#    n_45_135 = ma.ceil(mini_10[4])
#
#    if constraints.rule_10_percent and constraints.percent_45_135:
#        missing_extra_45_135 = ma.ceil(mini_10[4]) \
#        - ma.ceil(mini_10[2]) - ma.ceil(mini_10[3])
#    else:
#        missing_extra_45_135 = 0
#
#    print('n_45', n_45)
#    print('n_135', n_135)
#    print('n_45_135', n_45_135)
#    print('missing_extra_45_135', missing_extra_45_135)
#
#    print('\n\n*** Test for the function calc_current_10_2 ***')
#    ss = np.array([
#        45, 45, 90, 45, 45, 90, -45, 0, 0, -45, 90,
#        -45, 0, 0, -45, 90, 45, 45, 90, 45, 45], int)
#    print('\nInitial stacking sequence')
#    print_ss(ss, 40)
#    current_10 = calc_current_10_2(ss, sym=constraints.sym)
#    print('\ncurrent_10', current_10)
#
    print('\n*** Test for the function repair_10_bal***')
    ss_ini = np.array([60, 45, 60, 0], int)
    ss_ini = np.array([60, 45, 60, 15, -15, 30, 45], int)
    ss_ini = np.array([-45, 45, -45, 0, 0, 0, -45, 90, 45, 45], int)

    if constraints.sym:
#        ss_ini = np.hstack((ss_ini, 0, np.flip(ss_ini)))
        ss_ini = np.hstack((ss_ini, np.flip(ss_ini)))

#    ss_ini = "60 15 60 -15 90 60"
#    ss_ini = np.array(ss_ini.split()).astype(int)

    print('\nInitial stacking sequence')
    print_ss(ss_ini, 2000)
    print('ss_ini.zize', ss_ini.size)
    mini_10 = calc_mini_10(constraints, ss_ini.size)
    print('mini_10', mini_10)
    ss, ply_queue = repair_10_bal(ss_ini, mini_10, constraints)
    print('\nSolution stacking sequence')
    print_ss(ss, 2000)
    print('ply_queue', ply_queue)
