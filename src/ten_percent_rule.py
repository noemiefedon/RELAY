# -*- coding: utf-8 -*-
"""
Application of the 10% rule for the design of a laminate

- display_ply_counts
    displays the ply counts in each fibre direction for a laminate lay-up

- is_ten_percent_rule
    returns True for a panel stacking sequence satisfying the 10% rule,
    False otherwise

- calc_penalty_10_ss and calc_penalty_10_pc
    returns the stacking sequence penalty for 10% rule

- ten_percent_rule
    returns only the stacking sequences that satisfy the 10% rule when added to
    plies for which the ply orientations have been previously determined

- ten_percent_rule_lampam
    returns only the stacking sequences for which the lamiantion peremeters
    satisfy the 10% rule with restrictions on lamination parameters when added
    to a stack of plies

- calc_n_plies_per_angle
    returns the ply counts in each fibre direction
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_list_ss

def display_ply_counts(stack, constraints):
    '''
    displays the ply counts in each fibre direction for a laminate lay-up

    INPUTS

    - stack (array): stacking sequence
    - constraints (instance of the class Constraints): set of design guidelines
    '''
    print('number of 0 plies: ', sum(stack == 0))
    for angle in constraints.angles_bal[:, 0]:
        print('number of +' + str(angle) + ' plies: ', sum(stack == angle))
        print('number of -' + str(angle) + ' plies: ', sum(stack == -angle))
    print('number of 90 plies: ', sum(stack == 90))
    return 0


def is_ten_percent_rule(
        constraints, stack=None, ply_queue=None, n_plies_per_angle=None,
        equality_45_135=False, equality_0_90=False):
    """
    checks the satisfaction to the 10% rule

    Args:
        constraints (instance of the class Constraints): set of design
            guidelines
        stack (numpy array): partial stacking sequence with the angle 666 used
            for unknown ply fibre orientations
        ply_queue (list): ply fibre orientations for the remaining plies in the
            stacking sequence
        n_plies_per_angle (numpy array): ply counts in each fibre orientation
        equality_45_135 (boolean): True if +45/-45 plies are not differentiated
        equality_0_90 (boolean): True if 0/90 plies are not differentiated

    Returns:
        boolean: True if the stacking sequence 'stack' satisfies the 10% rule,
        False otherwise.

    Examples:
        >>> constraints=Constraints(rule_10_percent=True, percent_0=50)
        >>> is_ten_percent_rule(constraints, stack=np.array([0, 45, 90], int))
        False
    """
    if n_plies_per_angle is not None:
        if constraints.percent_tot > 0:
            n_total = sum(n_plies_per_angle)
            percent_0 = n_plies_per_angle[constraints.index0] / n_total
            percent_45 = n_plies_per_angle[constraints.index45] / n_total
            percent_90 = n_plies_per_angle[constraints.index90] / n_total
            percent_135 = n_plies_per_angle[constraints.index135] / n_total
            if percent_0 < constraints.percent_0 \
            or percent_45 < constraints.percent_45 \
            or percent_90 < constraints.percent_90 \
            or percent_135 < constraints.percent_135 \
            or percent_45 + percent_135 < constraints.percent_45_135:
                return False
        return True

    if isinstance(stack, list):
        n_total = len(stack)
    else:
        n_total = stack.size

    if ply_queue is not None:

        if constraints.sym:
            percent_0 = 2 * (
                sum(stack[:stack.size // 2] == 0) + ply_queue.count(0))
            percent_45 = 2 * (
                sum(stack[:stack.size // 2] == 45) + ply_queue.count(45))
            percent_90 = 2 * (
                sum(stack[:stack.size // 2] == 90) + ply_queue.count(90))
            percent_135 = 2 * (
                sum(stack[:stack.size // 2] == -45) + ply_queue.count(-45))
            if stack.size % 2:
                mid_ply_angle = stack[stack.size % 2]
                if mid_ply_angle == 0:
                    percent_0 += 1
                if mid_ply_angle == 90:
                    percent_90 += 1
                if mid_ply_angle == 45:
                    percent_45 += 1
                if mid_ply_angle == -45:
                    percent_135 += 1
        else:
            percent_0 = sum(stack == 0) + ply_queue.count(0)
            percent_45 = sum(stack == 45) + ply_queue.count(45)
            percent_90 = sum(stack == 90) + ply_queue.count(90)
            percent_135 = sum(stack == -45) + ply_queue.count(-45)
    else:
        percent_0 = sum(stack == 0)
        percent_45 = sum(stack == 45)
        percent_90 = sum(stack == 90)
        percent_135 = sum(stack == -45)

    percent_0 /= n_total
    percent_45 /= n_total
    percent_90 /= n_total
    percent_135 /= n_total

#    print(percent_0, constraints.percent_0)
#    print(percent_90, constraints.percent_90)
#    print(percent_45, constraints.percent_45)
#    print(percent_135, constraints.percent_135)
#    print(percent_45 + percent_135, constraints.percent_45_135)

    if not equality_0_90 and (percent_0 + 1e-15 < constraints.percent_0\
                              or percent_90 + 1e-15 < constraints.percent_90):
        return False

    if equality_0_90 and (percent_0 + percent_90 + 1e-15 \
                          < constraints.percent_0 + constraints.percent_90):
        return False

    if not equality_45_135 and (percent_45 + 1e-15 < constraints.percent_45 \
                                or percent_135 + 1e-15 < constraints.percent_135):
        return False

    if equality_45_135 and (percent_45 + percent_135 + 1e-15 < \
                            constraints.percent_45 + constraints.percent_135):
        return False

    if percent_45 + percent_135 + 1e-15 < constraints.percent_45_135:
        return False

    return True


def calc_penalty_10_ss(ss, constraints):
    """
    returns the stacking sequence penalty for 10% rule

    INPUTS

    - ss: stacking sequences
    """
    if constraints.percent_tot > 0:
        ss = np.array(ss)
        n_total = ss.size
        percent_0 = np.sum(ss == 0)/n_total
        percent_45 = np.sum(ss == 45)/n_total
        percent_90 = np.sum(ss == 90)/n_total
        percent_135 = np.sum(ss == -45)/n_total
        return (max(0, constraints.percent_0 - percent_0)
                + max(0, constraints.percent_45 - percent_45)
                + max(0, constraints.percent_90 - percent_90)
                + max(0, constraints.percent_135 - percent_135)
                + max(0, constraints.percent_45_135 - percent_45-percent_135))\
        / constraints.percent_tot
    return 0


def calc_penalty_10_pc(n_plies_per_angle, constraints, cummul_areas=1):
    """
    returns the penalty for 10% rule based on n_plies_per_angle
    """
    if constraints.percent_tot > 0:
        if (isinstance(n_plies_per_angle, np.ndarray) \
            and n_plies_per_angle.ndim == 1)\
        or isinstance(n_plies_per_angle, list):
            n_total = sum(n_plies_per_angle)
            if n_total:
                percent_0 = n_plies_per_angle[constraints.index0]/n_total
                percent_45 = n_plies_per_angle[constraints.index45]/n_total
                percent_90 = n_plies_per_angle[constraints.index90]/n_total
                percent_135 = n_plies_per_angle[constraints.index135]/n_total
                return cummul_areas * (
                    max(0, constraints.percent_0 - percent_0)
                    + max(0, constraints.percent_45 - percent_45)
                    + max(0, constraints.percent_90 - percent_90)
                    + max(0, constraints.percent_135 - percent_135)
                    + max(0, constraints.percent_45_135 \
                          - percent_45 - percent_135))/constraints.percent_tot

#        print('n_plies_per_angle.shape', n_plies_per_angle.shape)
        penalties = np.zeros((n_plies_per_angle.shape[0],))
        for ind_ss in range(n_plies_per_angle.shape[0]):
            n_total = sum(n_plies_per_angle[ind_ss])
            if n_total:
                percent_0 = n_plies_per_angle[
                    ind_ss][constraints.index0]/n_total
                percent_45 = n_plies_per_angle[
                    ind_ss][constraints.index45]/n_total
                percent_90 = n_plies_per_angle[
                    ind_ss][constraints.index90]/n_total
                percent_135 = n_plies_per_angle[
                    ind_ss][constraints.index135]/n_total
                penalties[ind_ss] = (
                    max(0, constraints.percent_0 - percent_0)
                    + max(0, constraints.percent_45 - percent_45)
                    + max(0, constraints.percent_90 - percent_90)
                    + max(0, constraints.percent_135 - percent_135)
                    + max(0, constraints.percent_45_135 \
                          - percent_45 - percent_135))/constraints.percent_tot
    #            print('penalties.shape', penalties.shape)
        return cummul_areas * penalties
    return 0


def ten_percent_rule(
        n_plies_per_angle_ini, constraints, parameters,
        middle_ply, angle, angle2=None):
    '''
    returns only the stacking sequences that satisfy the 10% rule when added to
    plies for which the ply orientations have been previously determined

    ONLY USED BY LAYLA, NOT BELLA

    OUTPUTS

    - angle: the selected sublaminate stacking sequences
    - angle2: the selected sublaminate stacking sequences if a
    second sublaminate is given as input for angle2

    INPUTS
     - n_plies_per_angle_ini: number of initials plies per fibre
    orientation in the same order indicated in constraints.set_of_angles
     - constraints.rule_10_percent = True implies the 10% rule is active
     - middle_ply = 0 if there is no ply overlapping the mid-surface for the
    sublaminate considered here, otherwise, middle_ply is equal to the position
    number of this ply.
    - angle: the first sublaminate stacking sequences
    - angle:2 matrix storing the second sublaminate stacking sequences
    '''
    if angle.ndim == 1:
        angle = angle.reshape((1, angle.size))
    if constraints.rule_10_percent and parameters.penalty_10_pc_switch:
        # For S or U laminates
        if angle2 is None:
            size_ss = angle.shape[0]
            for ind_ss in np.arange(size_ss)[::-1]:
                n_plies_per_angle = np.copy(n_plies_per_angle_ini)
                for ind_ply in range(angle.shape[1]):
                    index = constraints.ind_angles_dict[angle[ind_ss, ind_ply]]
                    n_plies_per_angle[index] += 1
                # Remove the extra half contribution for potential middle ply
                if middle_ply != 0:
                    index = constraints.ind_angles_dict[angle[ind_ss, -1]]
                    n_plies_per_angle[index] -= 1/2
                n_total = sum(n_plies_per_angle)
                percent_0 = n_plies_per_angle[constraints.index0]/n_total
                percent_45 = n_plies_per_angle[constraints.index45]/n_total
                percent_90 = n_plies_per_angle[constraints.index90]/n_total
                percent_135 = n_plies_per_angle[constraints.index135]/n_total
#                print('percent_0', percent_0)
#                print('percent_90', percent_90)
#                print('percent_45', percent_45)
#                print('percent_135', percent_135)
                if percent_0 + 1e-15 < constraints.percent_0 \
                or percent_45 + 1e-15 < constraints.percent_45 \
                or percent_90 + 1e-15 < constraints.percent_90 \
                or percent_135 + 1e-15 < constraints.percent_135 \
                or percent_45 + percent_135 + 1e-15 < constraints.percent_45_135:
                    angle = np.delete(angle, np.s_[ind_ss], axis=0)
            return angle, angle2

        # For SB or B laminates
        size_ss = angle.shape[0]
        for ind_ss in np.arange(size_ss)[::-1]:
            n_plies_per_angle = np.copy(n_plies_per_angle_ini)
            for ind_ply in range(angle.shape[1]):
                index = constraints.ind_angles_dict[angle[ind_ss, ind_ply]]
                n_plies_per_angle[index] += 1
            for ind_ply in range(angle2.shape[1]):
                index = constraints.ind_angles_dict[angle2[ind_ss, ind_ply]]
                n_plies_per_angle[index] += 1
            if middle_ply != 0:
                index = constraints.ind_angles_dict[angle[ind_ss, -1]]
                n_plies_per_angle[index] -= 1/2
            n_total = sum(n_plies_per_angle)
            percent_0 = n_plies_per_angle[constraints.index0]/n_total
            percent_45 = n_plies_per_angle[constraints.index45]/n_total
            percent_90 = n_plies_per_angle[constraints.index90]/n_total
            percent_135 = n_plies_per_angle[constraints.index135]/n_total
            if percent_0 + 1e-15 < constraints.percent_0 \
            or percent_45 + 1e-15 < constraints.percent_45 \
            or percent_90 + 1e-15 < constraints.percent_90 \
            or percent_135 + 1e-15 < constraints.percent_135 \
            or percent_45 + percent_135 + 1e-15 < constraints.percent_45_135:
                angle = np.delete(angle, np.s_[ind_ss], axis=0)
                angle2 = np.delete(angle2, np.s_[ind_ss], axis=0)
    return angle, angle2


def calc_n_plies_per_angle(
        n_plies_per_angle_ini, constraints, middle_ply, angle, angle2=None):
    '''
    returns the ply counts in each fibre direction

    INPUTS
     - angle: the sublaminate stacking sequences
    - angle2: the second sublaminate stacking sequences
    - n_plies_per_angle_ini: number of initials plies per fibre
    orientation in the same order indicated in constraints.set_of_angles
     - constraints.rule_10_percent = True implies the 10% rule is active
    '''
    if angle.ndim == 1:
        angle = angle.reshape((1, angle.size))
    size_ss = angle.shape[0]
    n_plies_per_angle_tab = np.matlib.repmat(n_plies_per_angle_ini, size_ss, 1)
    if angle2 is None: # For S or U laminates
        for ind_ss in np.arange(size_ss)[::-1]:
            n_plies_per_angle = np.copy(n_plies_per_angle_ini)
            for ind_ply in range(angle.shape[1]):
                index = constraints.ind_angles_dict[angle[ind_ss, ind_ply]]
                n_plies_per_angle[index] += 1
            if middle_ply != 0:
                index = constraints.ind_angles_dict[angle[ind_ss, -1]]
                n_plies_per_angle[index] -= 1/2
            n_plies_per_angle_tab[ind_ss] = n_plies_per_angle
    else:
        for ind_ss in np.arange(size_ss)[::-1]:
            n_plies_per_angle = np.copy(n_plies_per_angle_ini)
            for ind_ply in range(angle.shape[1]):
                index = constraints.ind_angles_dict[angle[ind_ss, ind_ply]]
                n_plies_per_angle[index] += 1
            for ind_ply in range(angle2.shape[1]):
                index = constraints.ind_angles_dict[angle2[ind_ss, ind_ply]]
                n_plies_per_angle[index] += 1
            if middle_ply != 0:
                index = constraints.ind_angles_dict[angle[ind_ss, -1]]
                n_plies_per_angle[index] -= 1/2
            n_plies_per_angle_tab[ind_ss] = n_plies_per_angle
    return n_plies_per_angle_tab


def ten_percent_rule_lampam(lampam_tab_tab, constraints, angle, angle2=None):
    '''
    returns only the stacking sequences for which the lamiantion peremeters
    satisfy the 10% rule with restrictions on lamination parameters when added
    to a stack of plies

    OUTPUTS

    - lampam_tab_tab: the lamination parameters satisfing
    the modified 10% rule
    - angle: the matrix returned is equal to the input matrix with some rows
    deleted. These rows deleted are the same as the rows deleted in
    lampam_tab_tab
    - angle2: the matrix returned is equal to the input matrix angle2 with
    some rows deleted which are of same indices as the rows deleted in
    lampam_tab_tab

    INPUTS

    - lampam_tab_tab: the set of lamination parameters
    - constraints.rule_10_percent = true implies the 10% rule is active
    - angle: the stacking sequences related to the lamination
    parameters stored in lampam_tab_tab
    - angle2: another matrix storing stacking sequences
    '''
    if angle.ndim == 1:
        angle = angle.reshape((1, angle.size))

    if not constraints.rule_10_percent or parameters.penalty_10_pc_switch:
        return lampam_tab_tab, angle, angle2

    if angle2 is None:
        sizee = lampam_tab_tab.shape[0]
        for ind in np.arange(sizee)[::-1]:
            if lampam_tab_tab[ind, 2] + 1e-15 \
            < 2*lampam_tab_tab[ind, 0]+4*constraints.percent_90 - 1:
                lampam_tab_tab = np.delete(lampam_tab_tab, np.s_[ind], axis=0)
                angle = np.delete(angle, np.s_[ind], axis=0)

            elif lampam_tab_tab[ind, 2] + 1e-15\
            < -2*lampam_tab_tab[ind, 0]+4*constraints.percent_0 - 1:
                lampam_tab_tab = np.delete(lampam_tab_tab, np.s_[ind], axis=0)
                angle = np.delete(angle, np.s_[ind], axis=0)

            elif lampam_tab_tab[ind, 2] > 1 - 2*constraints.percent_45_135 + 1e-15:
                lampam_tab_tab = np.delete(lampam_tab_tab, np.s_[ind], axis=0)
                angle = np.delete(angle, np.s_[ind], axis=0)

    else:
        sizee = lampam_tab_tab.shape[0]
        for ind in np.arange(sizee)[::-1]:
            if lampam_tab_tab[ind, 2] + 1e-15 \
            < 2*lampam_tab_tab[ind, 0]+4*constraints.percent_90 -1:
                lampam_tab_tab = np.delete(lampam_tab_tab, np.s_[ind], axis=0)
                angle = np.delete(angle, np.s_[ind], axis=0)
                angle2 = np.delete(angle2, np.s_[ind], axis=0)

            elif lampam_tab_tab[ind, 2] + 1e-15\
            < -2*lampam_tab_tab[ind, 0]+4*constraints.percent_0 -1:
                lampam_tab_tab = np.delete(lampam_tab_tab, np.s_[ind], axis=0)
                angle = np.delete(angle, np.s_[ind], axis=0)
                angle2 = np.delete(angle2, np.s_[ind], axis=0)

            elif lampam_tab_tab[ind, 2] > 1 -2*constraints.percent_45_135 + 1e-15:
                lampam_tab_tab = np.delete(lampam_tab_tab, np.s_[ind], axis=0)
                angle = np.delete(angle, np.s_[ind], axis=0)
                angle2 = np.delete(angle2, np.s_[ind], axis=0)

    return lampam_tab_tab, angle, angle2


if __name__ == "__main__":
    constraints = Constraints(sym=True)
    constraints.rule_10_percent = True
    constraints.set_percentages(
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10)

    print('*** Test for the function display_ply_counts ***\n')
    display_ply_counts(stack=np.array([
       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  45,  90,  90,  90,  90,  90,  90, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45,  90,  90,  90,  90,  90,  90,  45,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], int), constraints=constraints)

    print('*** Test for the function is_ten_percent_rule ***\n')
    stack = np.array([
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  45,  90,  90,  90,  90,  90,  90, -45, -45, -45, -45, -45, -45, -45, -45, -45, -45,  90,  90,  90,  90,  90,  90,  45,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0], int)
    print(is_ten_percent_rule(constraints, stack=stack), '\n')

#    n_plies_per_angle = np.array([8., 1., 1., 1.])
#    print(is_ten_percent_rule(
#        constraints, n_plies_per_angle=n_plies_per_angle), '\n')

#    print('\n*** Test for the function calc_penalty_10_ss ***\n')
#    ss = [0, 45, 45, 90, -45]
#    print(calc_penalty_10_ss(ss, constraints))
#
#    print('\n*** Test for the function calc_penalty_10_pc ***\n')
#    n_plies_per_angle = [0, 0, 0, 0, 0, 10]
#    print(calc_penalty_10_pc(n_plies_per_angle, constraints))
#
#    print('\n*** Test for the function ten_percent_rule ***\n')
#    middle_ply = 0
#    n_plies_per_angle = np.array([0., 0., 0., 0.])
#    ss = np.array([[45, -45, 0, 45, 90],
#                   [45, 90, 45, 45, 45],
#                   [0, 45, 45, 45, 90],
#                   [90, 45, 45, 45, 0],
#                   [45, 45, 45, 90, 45],
#                   [45, 45, 0, -45, 45]])
#    print('Input stacking sequences:\n')
#    print_list_ss(ss)
#    test, _ = ten_percent_rule(n_plies_per_angle, constraints, parameters,
#                               middle_ply, ss)
#    if test.shape[0]:
#        print('Stacking sequences satisfying the rule:\n')
#        print_list_ss(test)
#    else:
#        print('No stacking sequence satisfy the rule\n')
#
#    print('\n*** Test for the function calc_n_plies_per_angle ***\n')
#    parameters.penalty_10_pc_switch = True
#    middle_ply = 0
#    n_plies_per_angle = np.array([0., 0., 0., 0.])
#    ss = np.array([[45, -45, 0, 45, 90],
#                   [45, 90, 45, 45, 45],
#                   [0, 45, 45, 45, 90],
#                   [90, 45, 45, 45, 0],
#                   [45, 45, 45, 90, 45],
#                   [45, 45, 0, -45, 45]])
#    print('Input stacking sequences:\n')
#    print_list_ss(ss)
#    n_plies_per_angle = calc_n_plies_per_angle(
#        n_plies_per_angle, constraints, middle_ply, ss)
#    print('ply counts', n_plies_per_angle)
#
#    print('\n*** Test for the function ten_percent_rule_lampam ***\n')
#    print('Input stacking sequences:\n')
#    parameters.penalty_10_pc_switch = False
#    middle_ply = 2
#    ss = np.array([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]])
#    lampam_tab_tab = np.array([
#        [0.15, 0.25, -0.5, 0.],
#        [0.15, 0.35, -0.5, 0.],
#        [0.05, 0.25, -0.5, 0.],
#        [0.05, 0.35, -0.5, 0.],
#        [0.2, 0.3, -0.4, 0.]])
#    print_list_ss(ss)
#    lampam_tab_tab, test, _ = ten_percent_rule_lampam(
#        lampam_tab_tab, constraints, ss)
#    if test.shape[0]:
#        print('Stacking sequences satisfying the rule:\n')
#        print_list_ss(test)
#    else:
#        print('No stacking sequence satisfy the rule\n')
#
#
