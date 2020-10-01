# -*- coding: utf-8 -*-
"""
Functions calculating lamination parameters
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def calc_lampam_from_delta_lp_matrix(stack, constraints, delta_lampams):
    """
    returns the lamination parameters of a laminate

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    - delta_lampams: ply partial lamination parameters
    """
    lampam = np.zeros((12,), float)
    for ind_ply in range(delta_lampams.shape[0]):
        lampam += delta_lampams[
            ind_ply, constraints.ind_angles_dict[stack[ind_ply]]]
    return lampam
