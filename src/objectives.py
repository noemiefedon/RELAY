# -*- coding: utf-8 -*-
"""
objective functions
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def objectives(
        lampam,
        targets,
        lampam_weightings,
        constraints,
        parameters,
        mat_prop=None):
    '''
    returns the objective function evaluated for a lamination
    parameter

    INPUTS

    - lampam: lamiation parameters
    - targets.lampam: target lamiation parameters
    - lampam_weightings: set of coefficients to be applied to the square terms
    of the norm 2, specific for each step of the algorithm
    - constraints: design and manufacturing constraints
    - parameters: optimiser parameters
    - mat_prop: material properties of the laminae
    '''
    if lampam.ndim == 1:

        if parameters.type_obj_func == 2:
            return sum(lampam_weightings*(lampam - targets.lampam)**2)
        elif parameters.type_obj_func == 1:
            return sum(lampam_weightings*abs(lampam - targets.lampam))
        else:
            raise Exception('TYpe of objective function not recognised')

    result = np.zeros((lampam.shape[0],), float)

    for ind in range(lampam.shape[0]):
        result[ind] = objectives(
            lampam=lampam[ind],
            targets=targets,
            lampam_weightings=lampam_weightings,
            constraints=constraints,
            parameters=parameters,
            mat_prop=mat_prop)
    return result
