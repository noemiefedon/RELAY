# -*- coding: utf-8 -*-
"""
repair for flexural properties
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.internal_diso_contig import internal_diso_contig
from src.pretty_print import print_lampam, print_ss, print_list_ss
from src.lp_functions_2 import calc_lampamD
from src.repair_flexural_tools import find_list_blocks
from src.repair_flexural_tools import calc_delta_lampamD_swap_1


def repair_flexural_1(ss_ini, out_of_plane_coeffs, lampam_target, constraints,
                      count_obj=False):
    """
    repair for flexural properties only accounting for one panel by
    swapping block of plies with the same ply counts to improve the
    out-of-plane lamination parameters' convergence.

    This repair step is required because despite the preference of modifying
    inner plies rather than outer plies during the repair strategy, the
    out-of-plane lamination parameters matching with their targets might have
    been reduced by the first steps of the repair process. The design
    guidelines remain satisfied by the laminates during this repair step.

    INPUTS

    - ss_ini: stacking sequence (array)
    - out_of_plane_coeffs: coefficients in the out-of-plane objective function
    - lampam_target: lamination parameter targets
    - constraints: design and manufacturing constraints
    - count_obj: flag to count the number of objective function calls
    """
    n_obj_func_D_calls = 0

    ss = np.copy(ss_ini)
    n_plies = ss.size
    list_blocks = find_list_blocks(ss, n_plies, constraints)
    if constraints.sym and ss.size % 2:
        midply = ss[ss.size // 2]
    else:
        midply = None

#    print('list_blocks', list_blocks)
#    print_ss(ss)

    status_dam_tol = np.zeros((ss.size,), bool)
    if constraints.dam_tol:
        if hasattr(constraints, 'dam_tol_rule') \
        and constraints.dam_tol_rule in {2, 3}:
            status_dam_tol[0:2] = 1
            status_dam_tol[-2:] = 1
        elif not hasattr(constraints, 'dam_tol_rule') and \
        constraints.n_plies_dam_tol == 2:
            status_dam_tol[0:2] = 1
            status_dam_tol[-2:] = 1
        elif hasattr(constraints, 'dam_tol_rule') \
        and constraints.dam_tol_rule == 1:
            status_dam_tol[0] = 1
            status_dam_tol[-1] = 1
        elif not hasattr(constraints, 'dam_tol_rule') and \
        constraints.n_plies_dam_tol == 1:
            status_dam_tol[0] = 1
            status_dam_tol[-1] = 1

    # current out-of-plane objective function value
    lampamD = calc_lampamD(ss, constraints)
    objD = sum(out_of_plane_coeffs*((lampamD - lampam_target[8:12])**2))
    n_obj_func_D_calls += 1
#    print('objD', objD)

    # swapping the outer plies in priority as they are more important for
    # out-of-plane properties
    for index_block_1 in range(len(list_blocks[:-1])):
        block_1 = list_blocks[index_block_1]
        if status_dam_tol[block_1.first_ply_pos]:
            continue
#        if block_1.distance_middle < 0.5:
#            break
        ind_possible_swap = []
        for index_block_2, block_2 \
        in enumerate(list_blocks[index_block_1 + 1:]):
            if status_dam_tol[list_blocks[index_block_2].first_ply_pos]:
                continue
            if block_1.n_block_of_plies == block_2.n_block_of_plies \
            and block_2.angle in block_1.possible_angles \
            and block_1.angle in block_2.possible_angles :
                ind_possible_swap.append(index_block_1 + index_block_2 + 1)
#        print('index_block_1', index_block_1)
#        print('ind_possible_swap', ind_possible_swap)
        if ind_possible_swap: # possible swaps
            new_objectives_D = 1e10*np.ones((len(ind_possible_swap),))
            new_lampamDs = np.zeros((len(ind_possible_swap), 4))
            index = 0
            for index_block_2 in ind_possible_swap:
                # change all the stacking sequence
                new_lampamDs[index] = lampamD + calc_delta_lampamD_swap_1(
                    angle_first=block_1.angle,
                    angle_second=list_blocks[index_block_2].angle,
                    n_first_ply=block_1.first_ply_pos,
                    n_second_ply=list_blocks[index_block_2].first_ply_pos,
                    n_plies_group=block_1.n_block_of_plies,
                    n_plies=n_plies,
                    constraints=constraints)
                new_objectives_D[index] = sum(out_of_plane_coeffs*(
                    (new_lampamDs[index] - lampam_target[8:12])**2))
                n_obj_func_D_calls += 1
                index += 1
#            print('new_objectives_D', new_objectives_D)

            if min(new_objectives_D) + 1e-20 < objD:
#                print('SWAP')

                ind_min = np.argmin(new_objectives_D)
                objD = new_objectives_D[ind_min]
                lampamD = new_lampamDs[ind_min]

                block_2 = list_blocks[ind_possible_swap[ind_min]]
#                print('block_1', block_1.ID, block_1.angle)
#                print('block_2', block_2.ID, block_2.angle)

                ss[block_1.first_ply_pos:
                   block_1.first_ply_pos + block_1.n_block_of_plies] \
                = block_2.angle
                ss[block_2.first_ply_pos:
                   block_2.first_ply_pos + block_2.n_block_of_plies] \
                = block_1.angle
                if constraints.sym:
                    ss[ss.size - block_1.first_ply_pos \
                       - block_1.n_block_of_plies:
                           ss.size - block_1.first_ply_pos] = block_2.angle
                    ss[ss.size - block_2.first_ply_pos \
                       - block_2.n_block_of_plies:
                           ss.size - block_2.first_ply_pos] = block_1.angle

#                print_ss(ss)

                block_2.angle, block_1.angle = block_1.angle, block_2.angle
                block_1.update_possible_angles(
                    constraints, list_blocks, midply)
                block_2.update_possible_angles(
                    constraints, list_blocks, midply)

#                print_ss(ss)
#                test, _ = internal_diso_contig(ss, constraints)
#                if test.size < 1:
#                    raise Exception(
#                        'Contiguity or disorientation rule not satisfied.')
#                print('list_blocks', list_blocks)

            if objD < 1e-10:
                if count_obj:
                    return ss, n_obj_func_D_calls
                return ss

    if count_obj:
        return ss, n_obj_func_D_calls
    return ss

if __name__ == "__main__":

    print('\n*** Test for the function repair_flexural_1 ***')
    constraints = Constraints(
        sym=True,
        ipo=True,
        dam_tol=False,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=4)
    lampam_target = np.array([1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0])
    ss = np.array([90, 90, 45, 0, 0, -45, 90, 45, 0, -45], int)
    if constraints.sym:
        ss = np.hstack((ss, np.flip(ss)))
    print('\nInitial stacking sequence')
    print_ss(ss)
    test, _ = internal_diso_contig(ss, constraints)
    if test.size < 1:
        raise Exception('Contiguity or disorientation rule not satisfied.')

    out_of_plane_coeffs = np.array([1/3, 1/3, 1/3, 0])

    ss = repair_flexural_1(
        ss, out_of_plane_coeffs, lampam_target, constraints)
    print('\nFinal stacking sequence')
    print_ss(ss)
    test, _ = internal_diso_contig(ss, constraints)
    if test.size < 1:
        raise Exception('Contiguity or disorientation rule not satisfied.')