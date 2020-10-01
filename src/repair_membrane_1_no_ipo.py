# -*- coding: utf-8 -*-
"""
- calc_objA_options
    calculates the possible in-plane objective function values achievable by
    modifying one fibre orientation

- repair_membrane_1_no_ipo:
    repair for membrane properties only accounting for one panel when the
    laminate does not have to remain balanced

- calc_lampamA_options_3
    calculates the elementary changes of in-plane lamination parameters
    when modifying a fibre orientation

- calc_objA_options_3
    calculates the possible in-plane objective function values achievable by
    modifying one fibre orientation

- calc_ind_plies
    makes:
        - a list of all ply indices which can be modified during the refinement
        for membrane properties, sorted by starting with innermost plies
        - a list of the ply indices in each fibre direction which can be
        modified during the refinement for membrane properties, sorted by
        starting with innermost plies

- calc_lampamA_ply_queue
    calculates in-plane lamination parameters based on a partially retrieved
    stacking sequence and a list of plies for the innermost plies whose
    positions are left to be determined
    """
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import math as ma
import numpy as np

sys.path.append(r'C:\RELAY')
from src.sorting import sortAccording
from src.constraints import Constraints
from src.pretty_print import print_ss
from src.lampam_functions import calc_lampam
from src.repair_10_bal import calc_mini_10
from src.repair_10_bal import calc_current_10_2
from src.repair_tools import RepairError

def repair_membrane_1_no_ipo(
        ss_ini, ply_queue_ini, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints):
    """
    repair for membrane properties only accounting for one panel when the
    laminate does not have to remain balanced

    modifies the stacking sequence to converge towards the in-plane target
    lamination parameters. The modifications preserves the satisfaction to the
    10% rule, to the balance requirements and to the damage tolerance
    constraints.

    This step is required because step 1 and 2 have a significant effect on the
    in-plane lamination parameters.

    The fibre orientations are modified one by one.

    INPUTS

    - ss_ini: partially retrieved stacking sequence
    - ply_queue_ini: queue of plies for innermost plies
    - mini_10: number of plies required for the 10 % rule in the 0/90/45/-45
    fibre directions
    - in_plane_coeffs: coefficients in the in-plane objective function
    - p_A: coefficient for the proportion
        of the laminate thickness that can be modified during the repair
        for membrane properties
    - lampam_target: lamination parameter targets
    - constraints: design and manufacturing constraints
    - p_A: coefficient for the
    proportion of the laminate thickness that can be modified during the repair
    for membrane properties
    """
    n_plies = ss_ini.size

    ss = np.copy(ss_ini)
    ply_queue = ply_queue_ini[:]

    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)
    objA = sum(in_plane_coeffs * ((lampamA - lampam_target[0:4]) ** 2))
#    print('objA', objA)

    ss_list = [np.copy(ss)]
    ply_queue_list = [ply_queue[:]]
    lampamA_list = [lampamA]
    objA_list = [objA]

    excess_10 = calc_excess_10(ss, ply_queue, mini_10, constraints.sym)

    indices_1, indices_per_angle = calc_ind_plies(
        ss, n_plies, ply_queue, constraints, p_A)
    indices_to_sort = list(indices_1)
    indices_to_sort.insert(0, -1)
#    print('indices_1', list(indices_1))
#    print('indices_per_angle', list(indices_per_angle))
#    print('indices_to_sort', indices_to_sort)

    lampamA_options = calc_lampamA_options_3(n_plies, constraints)
    objA_options = calc_objA_options_3(
        lampamA, lampamA_options, lampam_target, constraints, in_plane_coeffs)
#    print('objA_options', objA_options)

    while np.min(objA_options) + 1e-20 < objA and objA > 1e-10:
        # attempts at modifying a couple of angled plies
        ind_angle1, ind_angle2 = np.unravel_index(
            np.argmin(objA_options, axis=None), objA_options.shape)
        angle1 = constraints.set_of_angles[ind_angle1]
        angle2 = constraints.set_of_angles[ind_angle2]
#        print('test  angle1', angle1, 'to angle2', angle2)
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle', indices_per_angle)

        # if no ply to be deleted
        if len(indices_per_angle[ind_angle1]) < 1:
            objA_options[ind_angle1, ind_angle2] = 1e10
            continue

        # attention to not break the 10% rule
        if angle1 == 0:
            if excess_10[0] < 1:
                objA_options[ind_angle1, ind_angle2] = 1e10
                continue
            excess_10[0] -= 1
        elif angle1 == 90:
            # attention to not break the 10% rule
            if excess_10[1] < 1:
                objA_options[ind_angle1, ind_angle2] = 1e10
                continue
            excess_10[1] -= 1
        elif angle1 == 45:
            # attention to not break the 10% rule
            if excess_10[2] < 1:
                objA_options[ind_angle1, ind_angle2] = 1e10
                continue
            excess_10[2] -= 1
        elif angle1 == -45:
            # attention to not break the 10% rule
            if excess_10[3] < 1:
                objA_options[ind_angle1, ind_angle2] = 1e10
                continue
            excess_10[3] -= 1

#        print(angle1, ' plies changed into ', angle2, 'plies')
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle[ind_angle1]', indices_per_angle[ind_angle1])
#        print('indices_per_angle[ind_angle2]', indices_per_angle[ind_angle2])

        if angle2 == 0:
            excess_10[0] += 1
        elif angle2 == 90:
            excess_10[1] += 1
        elif angle2 == 45:
            excess_10[2] += 1
        elif angle2 == -45:
            excess_10[3] += 1

        lampamA += lampamA_options[ind_angle2] - lampamA_options[ind_angle1]
        objA = objA_options[ind_angle1, ind_angle2]

        # modification of the stacking sequence
        ind_ply_1 = indices_per_angle[ind_angle1].pop(0)
#        print('ind_ply_1', ind_ply_1)

        if ind_ply_1 == 6666: # ply from the queue
            ply_queue.remove(angle1)
            ply_queue.append(angle2)
        else:
            ss[ind_ply_1] = angle2
            if constraints.sym:
                ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

        ss_list.insert(0, np.copy(ss))
        ply_queue_list.insert(0, ply_queue[:])
        lampamA_list.insert(0, np.copy(lampamA))
        objA_list.insert(0, objA)

        indices_per_angle[ind_angle2].append(ind_ply_1)
        if constraints.sym:
            indices_per_angle[ind_angle2].sort(reverse=True)
        else:
            sortAccording(indices_per_angle[ind_angle2], indices_to_sort)
            indices_per_angle[ind_angle2].reverse()

#        print('indices_per_angle', indices_per_angle)
#        print('objA', objA)
        if objA < 1e-10:
            break

        objA_options = calc_objA_options_3(
            lampamA, lampamA_options, lampam_target, constraints,
            in_plane_coeffs)
#        print('objA_options', objA_options)

    return ss_list, ply_queue_list, lampamA_list, objA_list


def calc_ind_plies(ss, n_plies, ply_queue, constraints, p_A):
    """
    makes:
        - a list of all ply indices which can be modified during the refinement
        for membrane properties, sorted by starting with innermost plies
        - a list of the ply indices in each fibre direction which can be
        modified during the refinement for membrane properties, sorted by
        starting with innermost plies

    Notes:
        - al lplies from the queue of plies are included.
        - middle plies of symmetric laminates are not included.
        - the rest of the plies are included only if they are part of the
        inner part of laminate representing p_A %
        of the overall laminate thickness.
    """
    ind_min = ma.floor(
        (1 - p_A/100)*(n_plies/2))
#    print('ind_min', ind_min)
    if constraints.sym:
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                indices_1 = range(max(2, ind_min), n_plies // 2)[::-1]
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                indices_1 = range(max(2, ind_min), n_plies // 2)[::-1]
            else:
                indices_1 = range(max(1, ind_min), n_plies // 2)[::-1]
        else:
            indices_1 = range(ind_min, n_plies // 2)[::-1]
    else:
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                ind_1 = list(range(max(2, ind_min), n_plies // 2)[::-1])
                ind_2 = list(range(
                    n_plies // 2, min(n_plies - ind_min, n_plies - 2)))
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                ind_1 = list(range(max(2, ind_min), n_plies // 2)[::-1])
                ind_2 = list(range(
                    n_plies // 2, min(n_plies - ind_min, n_plies - 2)))
            else:
                ind_1 = list(range(max(1, ind_min), n_plies // 2)[::-1])
                ind_2 = list(range(
                    n_plies // 2, min(n_plies - ind_min, n_plies - 1)))
        else:
            ind_1 = list(range(ind_min, n_plies // 2)[::-1])
            ind_2 = list(range(n_plies // 2, n_plies - ind_min))
        #print('ind_1', ind_1, 'ind_2', ind_2)
        indices_1 = np.zeros((len(ind_1) + len(ind_2),), 'int16')
        indices_1[::2] = ind_2
        indices_1[1::2] = ind_1
#    print('indices_1', list(indices_1))

    indices_per_angle = []
    for ind_angle in range(constraints.n_set_of_angles):
        indices_per_angle.append([])

    for ind_ply_1 in indices_1:
        if ss[ind_ply_1] != 666:
            indices_per_angle[
                constraints.ind_angles_dict[ss[ind_ply_1]]].append(ind_ply_1)

    for angle in ply_queue:
        indices_per_angle[constraints.ind_angles_dict[angle]].insert(0, 6666)

    return indices_1, indices_per_angle


def calc_lampamA_options_3(n_plies, constraints):
    """
    calculates the elementary changes of in-plane lamination parameters
    when modifying a fibre orientation
    """
    lampamA_options = np.empty((constraints.n_set_of_angles, 4), float)
    for ind_angle, angle in enumerate(constraints.set_of_angles):
        lampamA_options[ind_angle] = np.copy(
            constraints.cos_sin[ind_angle]).reshape(4)
    if not constraints.sym:
        lampamA_options *= 1 / n_plies
    else:
        lampamA_options *= 2 / n_plies
    return lampamA_options

def calc_objA_options_3(
        lampamA, lampamA_options, lampam_target, constraints, in_plane_coeffs):
    """
    calculates the possible in-plane objective function values achievable by
    modifying one fibre orientation

    objA_options[ind_pos_angle1, ind_pos_angle2] for angle1 ply changed
    to angle2 plies
    """
    objA_options = 1e10*np.ones((constraints.n_set_of_angles,
                                 constraints.n_set_of_angles), float)
    for ind_pos_angle1 in range(constraints.n_set_of_angles):
        for ind_pos_angle2 in range(constraints.n_set_of_angles):
            if ind_pos_angle1 == ind_pos_angle2:
                continue
            objA_options[ind_pos_angle1, ind_pos_angle2] = sum(
                in_plane_coeffs *((
                        lampamA \
                        - lampamA_options[ind_pos_angle1] \
                        + lampamA_options[ind_pos_angle2] \
                        - lampam_target[0:4])**2))
    return objA_options

def calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints):
    """
    calculates in-plane lamination parameters based on a partially retrieved
    stacking sequence and a list of plies for the innermost plies whose
    positions are left to be determined
    """
    cos_sin = np.zeros((4,), float)

    if not constraints.sym:
        for angle in ss:
            if angle != 666:
                cos_sin += constraints.cos_sin[
                    constraints.ind_angles_dict[angle]].reshape((4, ))

        for angle in ply_queue:
            cos_sin += constraints.cos_sin[
                constraints.ind_angles_dict[angle]].reshape((4, ))

        return (1 / n_plies) * cos_sin

    for angle in ss[:np.size(ss) // 2]:
        if angle != 666:
            cos_sin += constraints.cos_sin[
                constraints.ind_angles_dict[angle]].reshape((4, ))

    for angle in ply_queue:
        cos_sin += constraints.cos_sin[
            constraints.ind_angles_dict[angle]].reshape((4, ))

    if np.size(ss) % 2:
        cos_sin += 0.5 * constraints.cos_sin[
            constraints.ind_angles_dict[ss[n_plies // 2]]].reshape((4,))

    return (2 / n_plies) * cos_sin

def calc_excess_10(ss, ply_queue, mini_10, sym):
    """
returns the current number of plies in the 0/90/+45/-45 directions

    INPUTS

        ss: stacking sequence (array)
        sym: True for symmetric laminates (boolean)
    """
    ply_queue = np.array(ply_queue)
    current_10 = np.zeros((5,), float)
    if sym:
        lenn = ss.size // 2
        current_10[0] = sum(ss[:lenn] == 0) + sum(ply_queue == 0)
        current_10[1] = sum(ss[:lenn] == 90) + sum(ply_queue == 90)
        current_10[2] = sum(ss[:lenn] == 45) + sum(ply_queue == 45)
        current_10[3] = sum(ss[:lenn] == -45) + sum(ply_queue == -45)
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
        current_10[0] = sum(ss == 0) + sum(ply_queue == 0)
        current_10[1] = sum(ss == 90) + sum(ply_queue == 90)
        current_10[2] = sum(ss == 45) + sum(ply_queue == 45)
        current_10[3] = sum(ss == -45) + sum(ply_queue == -45)
        current_10[4] = current_10[2] + current_10[3]
    return current_10 - mini_10

if __name__ == "__main__":

    print('\n\n*** Test for the function calc_excess_10 ***')
    constraints = Constraints(
        sym=True,
        ipo=True,
        dam_tol=False,
        rule_10_percent=True,
        percent_0=10,
        percent_45=10,
        percent_90=10,
        percent_135=10,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    ss = np.array([0, 45, 666, 666, 666, 666, 666, 666, 666,
                   666, 666, 666, 666, 666, 666, 666, 45, 0], int)
    ply_queue = [90, 90, -45, 90, 90, 45, 0]
    mini_10 = calc_mini_10(constraints, ss.size)
    print('\nInitial stacking sequence')
    print_ss(ss, 40)
    excess_10 = calc_excess_10(ss, ply_queue, mini_10, sym=constraints.sym)
    print('\nexcess_10', excess_10)

    print('\n*** Test for the function calc_ind_plies ***')
    constraints = Constraints(
        sym=True,
        ipo=True,
        dam_tol=False,
        rule_10_percent=True,
        percent_0=10,
        percent_45=5,
        percent_90=10,
        percent_135=5,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    p_A = 50
    ss = np.array([
      45,  90,  45, 0, -45, 0, 666, 666], int)
    ss = np.hstack((ss, np.flip(ss)))
    ply_queue = [90, -45]
    n_plies = ss.size
    indices_1, indices_per_angle = calc_ind_plies(
        ss, n_plies, ply_queue, constraints, p_A)
    print('indices_1', list(indices_1))
    print('indices_per_angle', indices_per_angle)

    print('\n*** Test for the function calc_lampamA_ply_queue ***')
    constraints = Constraints(
        sym=False,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    ply_queue = [90, 90, 90, -45, -45, -45, 90, 90, 90, -45, -45, -45]
    n_plies = ss.size
    lampamA = calc_lampamA_ply_queue(
        ss, n_plies, ply_queue, constraints)

    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   90, 90, 90, -45, -45, -45, 90, 90, 90, -45, -45, -45,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    lampamA_check = calc_lampam(ss, constraints)[0:4]
    print('lampamA', lampamA)
    print('lampamA_check', lampamA_check)

    print()
    constraints = Constraints(
        sym=True,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    ply_queue = [90, 90, 90, -45, -45, -45]
    n_plies = ss.size
    lampamA = calc_lampamA_ply_queue(
        ss, n_plies, ply_queue, constraints)

    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   90, 90, 90, -45, -45, -45, 90, 90, 90, -45, -45, -45,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    lampamA_check = calc_lampam(ss, constraints)[0:4]
    print('lampamA', lampamA)
    print('lampamA_check', lampamA_check)


    print('\n*** Test for the function repair_membrane_1_no_ipo ***')
    constraints = Constraints(
        sym=True,
        ipo=False,
        dam_tol=False,
        rule_10_percent=True,
        percent_0=10,
        percent_45=10,
        percent_90=10,
        percent_135=10,
        set_of_angles=[0, 45, -45, 90])
#    set_of_angles=[0, 45, 30, -30, 60, -60, -45, 90])
    p_A = 100
    in_plane_coeffs = np.array([1, 1, 0, 0])
    ss_target = np.array([0], int)
    print('\nTarget stacking sequence')
    print_ss(ss_target, 40)
    lampam_target = calc_lampam(ss_target)
    ss = np.array([0, 45, 666, 666, 666, 666, 666, 666, 666,
                   666, 666, 666, 666, 666, 666, 666, 45, 0], int)
    ply_queue = [90, 90, -45, 90, 90, 45, 0]
#    ss = np.array([60, 45, 0, 0, 30, 0, -30,
#                   -30, 0, 30, 0, 0, 45, 60], int)
#    ply_queue = []
    print('\nInitial stacking sequence')
    print_ss(ss, 40)
    mini_10 = calc_mini_10(constraints, ss.size)
    ss_list, ply_queue_list, lampamA_list, objA_list = repair_membrane_1_no_ipo(
        ss, ply_queue, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints)
    print('\nSolution stacking sequences')
    for index in range(len(ss_list)):
        print_ss(ss_list[index], 20)
        print(ply_queue_list[index], 20)
