# -*- coding: utf-8 -*-
"""
Repair strategy

- repair_ss:
    attempts at repairing a stacking sequence for the following constraints:
        - damage tolerance
        - contiguity
        - disorientation
        - 10% rule
        - balance
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
import numpy.matlib

sys.path.append(r'C:\RELAY')
from src.parameters import Parameters
from src.constraints import Constraints
from src.objectives import objectives
from src.pretty_print import print_ss, print_list_ss
from src.repair_10_bal import repair_10_bal
from src.repair_10_bal import calc_mini_10
from src.repair_membrane import repair_membrane
from src.repair_flexural import repair_flexural
from src.repair_diso_contig import repair_diso_contig_list
from src.one_stack import check_ss_manufacturability
from src.lampam_functions import calc_lampam

def repair_ss(
        ss, constraints, parameters, lampam_target, obj_no_constraints=None,
        count_obj=False):
    """
    repairs stacking sequences to meet design and manufacturing guidelines
    and evaluates the performance of the repaired stacking sequence

    The repair process is deterministic and attempts at conducting minimal
    modification of the original stacking sequence with a preference for
    modifying outer plies that have the least influence on out-of-plane
    properties.

    step 1: repair for the 10% rule and balance
    step 2: refinement for in-plane lamination parameter convergence
    step 3: repair for disorientation and contiguity
    step 4: refinement for out-of-plane lamination parameter convergence
    (step 5: attribute a poor objective function value to unrepaired layups)

    OUTPUTS
    -

    INPUTS
    - ss: stacking sequence of the laminate
    - lampam_target: lamination parameter targets
    - constraints: instance of the class Constraints
    - parameters: instance of the class Parameters
    - count_obj: flag to count the number of objective function calls
    (- obj_no_constraints: objective function value of the initial stacking
     sequence with no consideration of design and manufacturing constraints)
    """
    ss_ini = np.copy(ss)
    mini_10 = calc_mini_10(constraints, ss.size)
#    print('before repair')
#    print_ss(ss_ini)
    #--------------------------------------------------------------------------
    # step 1 / repair for the 10% rule and balance
    #--------------------------------------------------------------------------
    ss, ply_queue = repair_10_bal(ss, mini_10, constraints)
#    print('after repair 10 and balance')
#    print_ss(ss)
#    print(ply_queue)
    #--------------------------------------------------------------------------
    # step 2 / improvement of the in-plane lamination parameter convergence
    #--------------------------------------------------------------------------
    ss_list, ply_queue_list, _ = repair_membrane(
        ss=ss,
        ply_queue=ply_queue,
        mini_10=mini_10,
        in_plane_coeffs=parameters.weighting_finalA,
        parameters=parameters,
        constraints=constraints,
        lampam_target=lampam_target)
#    print('after repair for membrane properties')
#    for ind in range(len(ss_list)):
#        print('ind', ind)
#        print('ss_list[ind]', ss_list[ind])
#        print('ply_queue_list[ind]', ply_queue_list[ind])
#        if not is_ten_percent_rule(constraints, stack=ss_list[ind],
#                                   ply_queue=ply_queue_list[ind]):
#            print('lampam_target', lampam_target[0:4])
#            raise Exception('10% rule not satisfied membrane')
#    print('ss_list[0]')
#    print_ss(ss_list[0])
#    print('ply_queue_list[0]', ply_queue_list[0])
    #--------------------------------------------------------------------------
    # step 3 / repair for disorientation and contiguity
    #--------------------------------------------------------------------------
    ss, completed_inward, completed_outward, ind = repair_diso_contig_list(
        ss_list, ply_queue_list, constraints,
        parameters.n_D1)
#    print('completed_inward, completed_outward, ind',
#          completed_inward, completed_outward, ind)
    if not completed_outward:
#        print('unsuccessful repair for disorientation and/or contiguity')
        if obj_no_constraints is None:
            if count_obj:
                return ss_ini, False, 0
            else:
                return ss_ini, False
        if count_obj:
            return ss_ini, False, 1e10, 0
        else:
            return ss_ini, False, 1e10
#    print('successful repair for disorientation and/or contiguity')
#    print_ss(ss)
    #--------------------------------------------------------------------------
    # step 4 / improvement of the out-of-plane lamination parameter convergence
    #--------------------------------------------------------------------------
    ss = repair_flexural(
        ss=ss,
        out_of_plane_coeffs=parameters.weighting_finalD,
        lampam_target=lampam_target,
        constraints=constraints,
        parameters=parameters,
        count_obj=count_obj)
    if count_obj:
        ss, n_obj_func_D_calls = ss
#    print('    after repair')
#    print_ss(ss)
#    print('lampam_target', lampam_target)

    if obj_no_constraints is None:
        if count_obj:
            return ss, True, n_obj_func_D_calls
        else:
            return ss, True
    #--------------------------------------------------------------------------
    # step 5 /
    #--------------------------------------------------------------------------
    obj_no_constraints = objectives(
        lampam=calc_lampam(ss, constraints),
        lampam_target=lampam_target,
        lampam_weightings=parameters.lampam_weightings_final,
        constraints=constraints,
        parameters=parameters)
    if count_obj:
        return ss, True, 1e10, n_obj_func_D_calls
    else:
        return ss, True, 1e10

if __name__ == "__main__":
    print('\n*** Test for the function repair_ss ***')
    constraints = Constraints(
        sym=True,
        bal=True,
        ipo=True,
        dam_tol=False,
        rule_10_percent=True,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=5,
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10,
        set_of_angles=[0, 45, -45, 90])
    ss = np.array([45,  90,  45,  90, -45, -45, -45, -45,  90,  45,  45,  45,
                   90, -45,   0,   0,   0, -45,  90,  45,  45,  45,  90, -45,
                   -45, -45, -45,  90,  45,  90,  45], int)
    ss_target = 60*np.ones((1,), dtype=int)
    lampam_target = calc_lampam(ss_target)
    #==========================================================================
    # Optimiser Parameters
    #==========================================================================
    ### Techniques to enforce the constraints
    # repair to improve the convergence towards the in-plane lamination parameter
    # targets
    repair_membrane_switch = True
    # repair to improve the convergence towards the out-of-plane lamination
    # parameter targets
    repair_flexural_switch = True
    # balanced laminate scheme
    balanced_scheme = False

    # coefficient for the proportion of the laminate thickness that can be modified
    # during the refinement for membrane properties in the repair process
    p_A = 80
    # number of plies in the last permutation during repair for disorientation
    # and/or contiguity
    n_D1 = 6
    # number of ply shifts tested at each step of the re-designing process during
    # refinement for flexural properties
    n_D2 = 10
    # number of times are redesigned during the refinement of flexural properties
    n_D3 = 2

    # Lamination parameters to be considered in the multi-objective functions
    optimisation_type = 'D'
    set_of_angles = np.array([-45, 0, 45, 90], int)
    if optimisation_type == 'A':
        if set_of_angles is np.array([-45, 0, 45, 90], int):
            lampam_to_be_optimised = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        else:
            lampam_to_be_optimised = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    if optimisation_type == 'D':
        if set_of_angles is np.array([-45, 0, 45, 90], int):
            lampam_to_be_optimised = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0])
        else:
            lampam_to_be_optimised = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
    if optimisation_type == 'AD':
        if set_of_angles is np.array([-45, 0, 45, 90], int):
            lampam_to_be_optimised = np.array([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0])
        else:
            lampam_to_be_optimised = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])

    # Lamination parameters sensitivities from the first-lebel optimiser
    first_level_sensitivities = np.ones((12,), float)

    parameters = Parameters(
        constraints=constraints,
        p_A=p_A,
        n_D1=n_D1,
        n_D2=n_D2,
        n_D3=n_D3,
        first_level_sensitivities=first_level_sensitivities,
        lampam_to_be_optimised=lampam_to_be_optimised,
        repair_membrane_switch=repair_membrane_switch,
        repair_flexural_switch=repair_flexural_switch)

    ss, completed, n_obj_func_D_calls = repair_ss(
        ss, constraints, parameters, lampam_target, count_obj=True)
    print('Repair successful?', completed)
    print_ss(ss, 20)
    print('n_obj_func_D_calls', n_obj_func_D_calls)
    check_ss_manufacturability(ss, constraints)
