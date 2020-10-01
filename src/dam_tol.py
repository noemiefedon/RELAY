# -*- coding: utf-8 -*-
"""
Functions related to the damage tolerance constraint

- is_damtol
    checks if a stacking sequence satisfies the damage tolerance guideline
    """
__version__ = '2.0'
__author__ = 'Noemie Fedon'


def is_dam_tol(stack, constraints):
    """
    checks if a stacking sequence satisfies the damage tolerance guideline

    Args:
        stack (numpy array): stacking sequence
        constraints (instance of the class Constraints): set of design
            guidelines

    Returns:
        boolean: True if the stacking sequence 'stack' satisfies the damage
            tolerance guideline, False otherwise.

    Examples:
        >>> is_damtol(np.array([0, 45, 90], int), Constraints(damtol=True))
        False

        >>> is_balanced(np.array([45, -45, 0, 45, -45], int), Constraints(damtol=True))
        True
    """
    if not constraints.dam_tol:
        return True

    my_set = set([45, -45])

    if hasattr(constraints, 'dam_tol_rule'):

        if constraints.dam_tol_rule == 1:
            if not stack[0] in my_set or not stack[-1] in my_set:
                return False
        elif constraints.dam_tol_rule == 2:
            if not stack[0] in my_set or not stack[1] in my_set \
            or not stack[-1] in my_set or not stack[-2] in my_set \
            or not stack[0] == - stack[1] or not stack[-1] == - stack[-2]:
                return False
        elif constraints.dam_tol_rule == 3:
            if not stack[0] in my_set or not stack[1] in my_set \
            or not stack[-1] in my_set or not stack[-2] in my_set:
                return False
        else:
            raise Exception('Not coded')

    else:
        if constraints.n_plies_dam_tol == 1:
            if not stack[0] in my_set or not stack[-1] in my_set:
                return False
        elif constraints.n_plies_dam_tol == 2:
            if not stack[0] in my_set or not stack[1] in my_set \
            or not stack[-1] in my_set or not stack[-2] in my_set \
            or not stack[0] == - stack[1] or not stack[-1] == - stack[-2]:
                return False
        else:
            raise Exception('Not coded')

    return True
