# -*- coding: utf-8 -*-
"""
repair tools

- RepairError:
    class for errors occuring when repairing a stacking sequence

- objective_repair:
    calculates the objective function values to chose the plies to change to
    improve the balance level of a laminate or its satisfaction to the 10% rule

- calc_delta_angle:
    calculates the angle difference between two angles given in degrees

- calc_distanciation
    if one ply position as input:
        returns the normalised distance between a ply and the middle surface
    if two ply positions as input:
        returns the difference of distances for the two plies regarding to the
        middle surface
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np

class RepairError(Exception):
    " Errors during the repair of a multi-panel design"
    pass

def objective_repair(
        delta_angle, ind_ply_1, ind_ply_2, n_plies, status_bal=0,
        coeff_bal_status=0.5, coeff_delta_angle=0.5, coeff_position=0.5):
    """
    calculates the objective function values to chose the plies to change to
    improve the balance level of a laminate or its satisfaction to the 10% rule
    (to minimise)

    delta_angle / 180: normalised value for the difference in fibre orientation

    2 * diff_position / n_plies: normalised value expressing the ply proximity
        with the middle surface

    INPUTS
        delta_angle: difference of fibre orientation between two plies
        ind_ply_1: index of the first ply
        ind_ply_2: index of the second ply
        n_plies: number of plies in the stacking sequence
        status_bal: staus of balance for the ply orientation of the second ply
        coeff_bal_status: coefficient for the preference of picking plies in
            excess and not in shortage for regarding the balance constraint
        coeff_delta_angle: coefficient for the preference of picking plies with
            small difference in fibre orientation
        coeff_position: coefficient for the preference of picking plies
            regarding its position in the stack
    """
    return - coeff_bal_status * status_bal \
+ coeff_delta_angle * delta_angle / 180 \
+ coeff_position * calc_distanciation(ind_ply_1, n_plies, ind_ply_2)


def calc_delta_angle(angle1, angle2):
    "calculates the angle difference between two angles given in degrees"
    return min(abs(angle2 - angle1),
               abs(180 + angle2 - angle1),
               abs(angle2 - 180 - angle1))

def calc_delta_angle_2(angle1, angle2):
    "calculates the angle difference between two angles given in degrees"
    my_array = np.array([
            angle2 - angle1, 180 + angle2 - angle1, angle2 - 180 - angle1])
    result = my_array[np.argmin(abs(my_array))]
    if result == -90:
        return 90
    return result


def calc_distanciation(ind_ply_1, n_plies, ind_ply_2=None):
    """
    if no value gien to ind_ply_2 :
        returns the normalised distance between the ply and the middle surface
    otherwise:
        returns the difference of distances for the two plies regarding to the
        middle surface
    """
    if ind_ply_1 < n_plies/2:
        position_1 = ind_ply_1
    else:
        position_1 = n_plies - ind_ply_1 - 1

    if ind_ply_2 is None:
        position_2 = n_plies / 2
    else:
        if ind_ply_2 < n_plies/2:
            position_2 = ind_ply_2
        else:
            position_2 = n_plies - ind_ply_2 - 1
    return 2 * abs(position_1 - position_2) / n_plies

if __name__ == "__main__":
    print('\n*** Test for the function calc_distanciation ***')
    print('calc_distanciation(2, 40)', calc_distanciation(2, 40))
    print('calc_distanciation(18, 40)', calc_distanciation(18, 40))

    print('\n*** Test for the function calc_delta_angle_2 ***')
    print('calc_delta_angle_2(-60, 30)', calc_delta_angle_2(-60, 30))
    print('calc_delta_angle_2(30, -60)', calc_delta_angle_2(30, -60))
    print('calc_delta_angle_2(-45, 30)', calc_delta_angle_2(-45, 30))
    print('calc_delta_angle_2(30, -45)', calc_delta_angle_2(30, -45))
    print('calc_delta_angle_2(-45, 45)', calc_delta_angle_2(-45, 45))
    print('calc_delta_angle_2(75, -75)', calc_delta_angle_2(75, -75))
    print('calc_delta_angle_2(-75, 75)', calc_delta_angle_2(-75, 75))