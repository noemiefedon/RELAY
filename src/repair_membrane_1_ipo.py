# -*- coding: utf-8 -*-
"""
- calc_lampamA_options_1
    calculates the elementary changes of in-plane lamination parameters
    when modifying the fibre orientations of couples of angled plies

- calc_objA_options_1
    calculates the possible in-plane objective function values achievable by
    modifying the fibre orientations of couples of angled plies

- repair_membrane_1_ipo:
    repair for membrane properties only accounting for one panel when the
    laminate must remain balanced
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_ss
from src.lampam_functions import calc_lampam
from src.repair_10_bal import calc_mini_10
from src.repair_membrane_1_no_ipo import calc_excess_10
from src.repair_membrane_1_no_ipo import calc_ind_plies
from src.repair_membrane_1_no_ipo import calc_lampamA_ply_queue
from src.ten_percent_rule import is_ten_percent_rule
from src.sorting import sortAccording

def repair_membrane_1_ipo(
        ss_ini, ply_queue_ini, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints):
    """
    modifies the stacking sequence to converge towards the in-plane target
    lamination parameters. The modifications preserves the satisfaction to the
    10% rule, to the balance requirements and to the damage tolerance
    constraints.

    This step is required because step 1 and 2 have a significant effect on the
    in-plane lamination parameters.

    First the fibre orientations of couples of angled plies situated in the
    middle of the laminate are modified to ensure convergence for the
    in-plane lamination parameters.
    Then the fibre orientations of 0 and 90 deg plies may be modified to the
    other orientation.

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
    - in_plane_coeffs: coefficients in the in-plane objective function
    - p_A: coefficient for the proportion
        of the laminate thickness that can be modified during the repair
        for membrane properties
    """
    n_plies = ss_ini.size

    ss = np.copy(ss_ini)
    ply_queue = ply_queue_ini[:]

#    print()
#    print('beginning repair membrane')
#    print_ss(ss)
#    print(ply_queue)
#    print('lampam_target', lampam_target[0:4])
#    print('mini_10', mini_10)
#    print('in_plane_coeffs', in_plane_coeffs)
#    print('p_A',
#          p_A)
#    print(constraints)
#    print()

    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)
    objA = sum(in_plane_coeffs * ((lampamA - lampam_target[0:4]) ** 2))
#    print()
#    print('lampamA', lampamA)
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

    lampamA_options = calc_lampamA_options_1(n_plies, constraints)
    objA_options = calc_objA_options_1(
        lampamA, lampamA_options, lampam_target, constraints, in_plane_coeffs)


    while np.min(objA_options) + 1e-20 < objA  and objA > 1e-10:
#        print('objA', objA)
#        print('objA_options', objA_options)

        # attempts at modifying a couple of angled plies
        ind_pos_angle1, ind_pos_angle2 = np.unravel_index(
            np.argmin(objA_options, axis=None), objA_options.shape)
        angle1 = constraints.pos_angles[ind_pos_angle1]
        angle2 = constraints.pos_angles[ind_pos_angle2]
#        print('angle1', angle1, 'angle2', angle2)
        ind_angle1 = constraints.ind_angles_dict[angle1]
        ind_angle1_minus = constraints.ind_angles_dict[-angle1]
        ind_angle2 = constraints.ind_angles_dict[angle2]
        ind_angle2_minus = constraints.ind_angles_dict[-angle2]
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)

#        print('indices_per_angle', indices_per_angle)

        # if plies +-theta exist in the middle of the laminate
        if angle1 == 0:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 2:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
            # attention to not break the 10% rule
            if excess_10[0] < 2:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
            excess_10[0] -= 2
        elif angle1 == 90:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 2:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
            # attention to not break the 10% rule
            if excess_10[1] < 2:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
            excess_10[1] -= 2
        elif angle1 == 45:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 1 \
            or len(indices_per_angle[ind_angle1_minus]) < 1:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
            # attention to not break the 10% rule
            if excess_10[2] < 1 or excess_10[3] < 1:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
            excess_10[2] -= 1
            excess_10[3] -= 1
        else:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 1 \
            or len(indices_per_angle[ind_angle1_minus]) < 1:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue

#        print('objA_options after clean', objA_options)
#        print('+-', angle1, ' plies changed into +-', angle2, 'plies')
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle[ind_angle1]', indices_per_angle[ind_angle1])
#        print('indices_per_angle[ind_angle2]', indices_per_angle[ind_angle2])

        if angle2 == 0:
            excess_10[0] += 2
        elif angle2 == 90:
            excess_10[1] += 2
        elif angle2 == 45:
            excess_10[2] += 1
            excess_10[3] += 1

        lampamA += lampamA_options[ind_pos_angle2] - lampamA_options[
            ind_pos_angle1]
        objA = objA_options[ind_pos_angle1, ind_pos_angle2]
#        print()
#        print('lampamA', lampamA)
#        print('objA', objA)

        # modification of the stacking sequence
        ind_ply_1 = indices_per_angle[ind_angle1].pop(0)
        ind_ply_2 = indices_per_angle[ind_angle1_minus].pop(0)
#        print('ind_ply_1', ind_ply_1)
#        print('ind_ply_2', ind_ply_2)
#        print('ply_queue', ply_queue)

        if ind_ply_1 == 6666: # ply from the queue
            ply_queue.remove(angle1)
            ply_queue.append(angle2)
        else:
            ss[ind_ply_1] = angle2
            if constraints.sym:
                ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

        if ind_ply_2 == 6666: # ply from the queue
            if angle1 == 90:
                ply_queue.remove(90)
            else:
                ply_queue.remove(-angle1)

            if angle2 == 90:
                ply_queue.append(90)
            else:
                ply_queue.append(-angle2)
        else:
            if angle2 != 90:
                ss[ind_ply_2] = -angle2
            else:
                ss[ind_ply_2] = 90
            if constraints.sym:
                ss[ss.size - ind_ply_2 - 1] = ss[ind_ply_2]

#        lampamA_check = calc_lampamA_ply_queue(
#            ss, ss.size, ply_queue, constraints)

        ss_list.insert(0, np.copy(ss))
        ply_queue_list.insert(0, ply_queue[:])
        lampamA_list.insert(0, np.copy(lampamA))
        objA_list.insert(0, objA)

        indices_per_angle[ind_angle2].append(ind_ply_1)
        indices_per_angle[ind_angle2_minus].append(ind_ply_2)
        if constraints.sym:
            indices_per_angle[ind_angle2].sort(reverse=True)
            indices_per_angle[ind_angle2_minus].sort(reverse=True)
        else:
            sortAccording(indices_per_angle[ind_angle2], indices_to_sort)
            sortAccording(indices_per_angle[ind_angle2_minus], indices_to_sort)
            indices_per_angle[ind_angle2].reverse()
            indices_per_angle[ind_angle2_minus].reverse()

#        print('indices_per_angle', indices_per_angle)
#        print('objA', objA)
        if objA < 1e-10:
            break

        objA_options = calc_objA_options_1(
            lampamA, lampamA_options, lampam_target, constraints,
            in_plane_coeffs)

#    print('objA', objA)

    # attempt at changing a 0 deg ply into a 90 deg ply
    ind_0 = np.where(constraints.pos_angles == 0)[0][0]
    ind_90 = np.where(constraints.pos_angles == 90)[0][0]

    if indices_per_angle[constraints.ind_angles_dict[0]] \
    and excess_10[0] > 0.5:
        obj_0_to_90 = sum(in_plane_coeffs*((
            (lampamA_options[ind_90] - lampamA_options[ind_0])/2 \
            + lampamA - lampam_target[0:4])**2))
#        print('obj_0_to_90', obj_0_to_90)

        if obj_0_to_90 + 1e-20 < objA:
#            print('excess_10[0]', excess_10[0])
#            print('0 deg ply changed to 90 deg ply')
            objA = obj_0_to_90
            lampamA += (lampamA_options[ind_90] - lampamA_options[ind_0])/2
            ind_ply_1 = indices_per_angle[constraints.index0].pop(0)
            if ind_ply_1 == 6666: # ply from the queue
                ply_queue.remove(0)
                ply_queue.append(90)
            else:
                ss[ind_ply_1] = 90
                if constraints.sym:
                    ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

            ss_list.insert(0, np.copy(ss))
            ply_queue_list.insert(0, ply_queue[:])
            lampamA_list.insert(0, np.copy(lampamA))
            objA_list.insert(0, objA)
#            print()
#            print('lampamA', lampamA)
#            print('objA', objA)

            return ss_list, ply_queue_list, lampamA_list, objA_list

    # attempt at changing a 90 deg ply into a 0 deg ply
    if indices_per_angle[constraints.ind_angles_dict[90]] \
    and excess_10[1] > 0.5:

        obj_90_to_0 = sum(in_plane_coeffs * ((
            (lampamA_options[ind_0] - lampamA_options[ind_90])/2 \
            + lampamA - lampam_target[0:4])**2))
#        print('obj_90_to_0', obj_90_to_0)

        if obj_90_to_0 + 1e-20 < objA:
#            print('90 deg ply changed to 0 deg ply')
            objA = obj_90_to_0
            lampamA += (lampamA_options[ind_0] - lampamA_options[ind_90])/2
            ind_ply_1 = indices_per_angle[constraints.index90].pop(0)
            if ind_ply_1 == 6666: # ply from the queue
                ply_queue.remove(90)
                ply_queue.append(0)
            else:
                ss[ind_ply_1] = 0
                if constraints.sym:
                    ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

            ss_list.insert(0, np.copy(ss))
            ply_queue_list.insert(0, ply_queue[:])
            lampamA_list.insert(0, np.copy(lampamA))
            objA_list.insert(0, objA)
#            print()
#            print('lampamA', lampamA)
#            print('objA', objA)

    return ss_list, ply_queue_list, lampamA_list, objA_list


def calc_lampamA_options_1(n_plies, constraints):
    """
    calculates the elementary changes of in-plane lamination parameters
    when modifying the fibre orientations of couples of angled plies
    """
    lampamA_options = np.empty((constraints.pos_angles.size, 4), float)
    for ind_angle, angle in enumerate(constraints.pos_angles):
        lampamA_options[ind_angle] = np.copy(constraints.cos_sin[
            constraints.ind_angles_dict[angle]]).reshape(4)
        lampamA_options[ind_angle] += constraints.cos_sin[
            constraints.ind_angles_dict[-angle]].reshape(4)
    if not constraints.sym:
        lampamA_options *= (1 / n_plies)
    else:
        lampamA_options *= (2 / n_plies)
    return lampamA_options


def calc_objA_options_1(
        lampamA, lampamA_options, lampam_target, constraints, in_plane_coeffs):
    """
    calculates the possible in-plane objective function values achievable by
    modifying the fibre orientations of couples of angled plies
    objA_options[ind_pos_angle1, ind_pos_angle2] for +-angle1 plies changed
    to +-angle2 plies

    The plies to be modified are chosen among the innermost layers to
    reduce the out-of-plane modifications
    """
    objA_options = 1e10*np.ones((constraints.pos_angles.size,
                                 constraints.pos_angles.size), float)
    for ind_pos_angle1 in range(constraints.pos_angles.size):
        for ind_pos_angle2 in range(constraints.pos_angles.size):
            if ind_pos_angle1 == ind_pos_angle2:
                continue
            objA_options[ind_pos_angle1, ind_pos_angle2] = sum(
                in_plane_coeffs * ((
                    lampamA \
                    - lampamA_options[ind_pos_angle1] \
                    + lampamA_options[ind_pos_angle2] \
                    - lampam_target[0:4]) ** 2))
    return objA_options

if __name__ == "__main__":

    print('\n*** Test for the function repair_membrane_1_ipo ***')
    constraints = Constraints(
        sym=True,
        bal=True,
        ipo=True,
        dam_tol=False,
        rule_10_percent=True,
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10,
        set_of_angles=[0, 45, -45, 90])
#    set_of_angles=[0, 45, 30, -30, 60, -60, -45, 90])
    p_A = 70
    in_plane_coeffs = np.array([0.5, 0.5, 0, 0])
    ss_target = np.array([0], int)
    lampam_target = calc_lampam(ss_target)

    ss = np.array([45, 90, 45, 45, 0, -45, -45, 0, 666, 666], int)

    ss = np.hstack((ss, np.flip(ss)))
    ply_queue = [90, -45]

    print('\nInitial stacking sequence')
    print_ss(ss, 200)
    print('ply_queue', ply_queue)
    mini_10 = calc_mini_10(constraints, ss.size)
    print('mini_10', mini_10)

    ss_list, ply_queue_list, lampamA_list, objA_list = repair_membrane_1_ipo(
        ss, ply_queue, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints)

#    print('\nSolution stacking sequence')
#    print_ss(ss_list[0], 200)
#    print('ply_queue', ply_queue_list[0])

#    lampamA = lampamA_list[0]
#    lampamA_check = calc_lampamA_ply_queue(ss, ss.size, ply_queue, constraints)
#    if not (abs(lampamA_check - lampamA) < 1e-10).all():
#        print()
#        print('lampamA      ', lampamA)
#        print('lampamA_check', lampamA_check)
#        print('lampam_target', lampam_target)

#    print('\nSolution stacking sequence 0')
#    print_ss(ss_list[0], 200)
#    print(ply_queue_list[0], 200)

    if not is_ten_percent_rule(
            constraints, stack=ss_list[0], ply_queue=ply_queue_list[0]):
        raise Exception('10% rule not satisfied membrane')

    for ind in range(len(ss_list)):
        print()
        print('ss_list[ind]', ss_list[ind])
        print('ply_queue_list[ind]', ply_queue_list[ind])
        if not is_ten_percent_rule(
                constraints, stack=ss_list[ind],
                ply_queue=ply_queue_list[ind]):
            raise Exception('10% rule not satisfied membrane')
