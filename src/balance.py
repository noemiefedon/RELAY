# -*- coding: utf-8 -*-
"""
Functions related to orthotropy requirements

- calc_penalty_bal
    calculates penalties for the balance constraint

- is_balanced
    checks if a stacking sequence is balanced
    """
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints

def is_balanced(stack, constraints):
    """
    checks if a stacking sequence is balanced

    Args:
        stack (numpy array): stacking sequence
        constraints (instance of the class Constraints): set of design
            guidelines

    Returns:
        boolean: True if the stacking sequence 'stack' is balanced, False
        otherwise.

    Examples:
        >>> is_balanced(np.array([0, 45, 90], int), Constraints())
        False

        >>> is_balanced(np.array([0, 45, 90, -45], int), Constraints())
        True
    """
    my_set = set(np.abs(stack))
    my_set.discard(0)
    my_set.discard(90)
    for angle in my_set:
        if sum(stack == angle) != sum(stack == -angle):
            return False
    return True

def calc_penalty_bal(n_plies_per_angle, constraints, cummul_areas=1):
    """
    calculates penalties for the balance constraint

    Args:
        n_plies_per_angle (numpy array): ply counts for each fibre direction
            for one or more laminates
        constraints (instance of the class Constraints): set of design
            guidelines.
            constraints.angles_bal (numpy array) has three columns for the:
                 - fibre orientations of angle ply theta
                 - indices of +theta in constraints.set_of_angles
                 - indices of -theta in constraints.set_of_angles
        cummul_areas (scalar): sum of the areas of the plies retrieved so far
            + the current plies

    Returns:
        the penalties for the balance constraint for one or several laminates
        based on ply counts given in each fibre direction, weighted by the
        coefficient 'cummul_areas'
    """
    if n_plies_per_angle.ndim == 1: # one laminate
        penalty_ipo_pc = 0
        my_sum = sum(n_plies_per_angle)
        if my_sum:
            for ind_angle in range(constraints.angles_bal.shape[0]):
                penalty_ipo_pc += abs(
                    n_plies_per_angle[constraints.angles_bal[ind_angle][1]] \
                    - n_plies_per_angle[constraints.angles_bal[ind_angle][2]])
            return cummul_areas * penalty_ipo_pc / sum(n_plies_per_angle)
        return 0

    # more than one laminate
    penalty_ipo_pc = np.zeros((n_plies_per_angle.shape[0],), float)
    for ind_ss in range(n_plies_per_angle.shape[0]):
        my_sum = sum(n_plies_per_angle[ind_ss])
        if my_sum:
            for ind_angle in range(constraints.angles_bal.shape[0]):
                penalty_ipo_pc[ind_ss] += abs(
                    n_plies_per_angle[ind_ss][
                        constraints.angles_bal[ind_angle][1]] \
                    - n_plies_per_angle[ind_ss][
                        constraints.angles_bal[ind_angle][2]])
            penalty_ipo_pc[ind_ss] /= my_sum

    return cummul_areas * penalty_ipo_pc

if __name__ == "__main__":

    print('\n*** Test for the functions is_balanced***\n')
    constraints = Constraints(ipo=True, oopo=True)
    print(is_balanced(
        stack = np.array([0, 45, 90, -45, 60, 30, 12, -12, 13]),
        constraints = Constraints()))

    print('\n*** Test for the functions calc_penalty_bal***\n')
    n_plies_per_angle = np.array([0, 0, 1, 0])
    constraints = Constraints(ipo=True, oopo=True)
    print(calc_penalty_bal(n_plies_per_angle, constraints))
    n_plies_per_angle = np.array([[2, 1, 2, 2], [2, 1, 2, 2]])
    constraints = Constraints(ipo=True, oopo=True)
    print(calc_penalty_bal(n_plies_per_angle, constraints))
