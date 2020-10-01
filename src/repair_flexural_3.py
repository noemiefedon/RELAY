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
from src.pretty_print import print_ss
from src.parameters import Parameters
from src.lampam_functions import calc_lampam

def repair_flexural_3(
        ss_ini, out_of_plane_coeffs, lampam_target, constraints,
        parameters, count_obj=False):
    """
    repair for flexural properties only accounting for one panel by
    inwardly redesigning the laminate to improve the out-of-plane lamination
    parameters' convergence.

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
    - parameters.n_D2: number of ply shifts tested at each step of the process
    - count_obj: flag to count the number of objective function calls
    """
    n_obj_func_D_calls = 0

    ss = np.copy(ss_ini)

    if constraints.sym:
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                ind_2 = np.arange(2, ss.size // 2)
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                ind_2 = np.arange(2, ss.size // 2)
            else:
                ind_2 = np.arange(1, ss.size // 2)
        else:
            ind_2 = np.arange(0, ss.size // 2)
    else:
        indices = np.arange(ss.size)
        beginning = np.copy(indices[0:indices.size // 2])
        ending = np.copy(indices[indices.size // 2:ss.size][::-1])
        ind_2 = np.zeros((ss.size,), int)
        ind_2[::2] = ending
        ind_2[1::2] = beginning
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                ind_2 = ind_2[4:]
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                ind_2 = ind_2[4:]
            else:
                ind_2 = ind_2[2:]
#    print('ind_2', ind_2)

    if constraints.sym:
        n_plies_here = ss.size // 2
        if ss.size % 2:
            index_ply = ss.size // 2
            norm_pos_top = 1 - (2 / ss.size) * (index_ply)
            norm_pos_bot =  1 - (2 / ss.size) * (index_ply + 1)
            sec_mom_areas_mid = (norm_pos_top**3 - norm_pos_bot**3) / 2
            delta_lampam_mid = sec_mom_areas_mid * constraints.cos_sin[
                    constraints.ind_angles_dict[ss[index_ply]]].reshape((4,))
        else:
            delta_lampam_mid = np.zeros((4,), float)
    else:
        n_plies_here = ss.size
        delta_lampam_mid = np.zeros((4,), float)

    cos_sin = np.zeros((n_plies_here, 4), float)
    sec_mom_areas = np.zeros((n_plies_here,), float)

    for index_ply in range(n_plies_here):
        cos_sin[index_ply] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[index_ply]]].reshape((4,))
        norm_pos_top = 1 - (2 / ss.size) * (index_ply)
        norm_pos_bot =  1 - (2 / ss.size) * (index_ply + 1)
        sec_mom_areas[index_ply] = (norm_pos_top**3 - norm_pos_bot**3)

    if not constraints.sym:
        sec_mom_areas /= 2

    lampamD = np.matmul(sec_mom_areas, cos_sin) + delta_lampam_mid
    objD = sum(out_of_plane_coeffs*((lampamD - lampam_target[8:12])**2))
    n_obj_func_D_calls += 1

    best_objD = objD
#    print('sec_mom_areas', sec_mom_areas, sum(sec_mom_areas))
#    print('cos_sin', cos_sin)
#    print('lampamD', lampamD)
#    print_ss(ss)
#    print('objD', objD)

    ind_step = 0
    for index, ind_ply_1 in enumerate(ind_2[:-1]):

#        print('ind_ply_1', ind_ply_1)
#        print('ind_2[index + 1:]', ind_2[index + 1:])

        objD = 1e10*np.ones((ind_2[index + 1:].size,), float)

        for ind_loc, ind_ply_2 \
        in enumerate(ind_2[index + 1:index + 2 + parameters.n_D2]):

#            print('    ind_ply_2', ind_ply_2)

            if not constraints.sym and ind_step % 2 == 0:

                ss_mod = np.copy(ss)
                ss_mod[ind_ply_2:ind_ply_1 + 1] = np.hstack((
                    ss_mod[ind_ply_2 + 1:ind_ply_1 + 1],
                    ss_mod[ind_ply_2]))
#                print('    ss_mod')
#                print_ss(ss_mod)

                # change
                cos_sin[ind_ply_2:ind_ply_1 + 1, :] = np.vstack((
                    cos_sin[ind_ply_2 +1:ind_ply_1 + 1, :],
                    cos_sin[ind_ply_2, :]))

                lampamD = np.matmul(sec_mom_areas, cos_sin) + delta_lampam_mid

                objD[ind_loc] = sum(out_of_plane_coeffs*((
                    lampamD - lampam_target[8:12])**2))
                n_obj_func_D_calls += 1

                # revert change
                cos_sin[ind_ply_2:ind_ply_1 + 1, :] = np.vstack((
                    cos_sin[ind_ply_1, :],
                    cos_sin[ind_ply_2:ind_ply_1, :]))

            else:
                # change
                cos_sin[ind_ply_1:ind_ply_2 + 1, :] = np.vstack((
                    cos_sin[ind_ply_2, :],
                    cos_sin[ind_ply_1:ind_ply_2, :]))

                lampamD = np.matmul(sec_mom_areas, cos_sin) + delta_lampam_mid
                objD[ind_loc] = sum(out_of_plane_coeffs*((
                    lampamD - lampam_target[8:12])**2))
                n_obj_func_D_calls += 1

                # revert change
                cos_sin[ind_ply_1:ind_ply_2 + 1, :] = np.vstack((
                    cos_sin[ind_ply_1 +1:ind_ply_2 + 1, :],
                    cos_sin[ind_ply_1, :]))

#        print('    objDs', objD)
        if min(objD) + 1e-20 < best_objD:

            best_objD = min(objD)
            ind_min = np.argmin(objD)
            ind_ply_2 = ind_2[index + 1:][ind_min]
#            print('    ind_ply_2 for change', ind_ply_2)

            if not constraints.sym and ind_step % 2 == 0:

                ss[ind_ply_2:ind_ply_1 + 1] = np.hstack((
                    ss[ind_ply_2 +1:ind_ply_1 + 1], ss[ind_ply_2]))

                cos_sin[ind_ply_2:ind_ply_1 + 1, :] = np.vstack((
                    cos_sin[ind_ply_2 +1:ind_ply_1 + 1, :],
                    cos_sin[ind_ply_2, :]))
            else:
                ss[ind_ply_1:ind_ply_2 + 1] = np.hstack((
                    ss[ind_ply_2], ss[ind_ply_1:ind_ply_2]))

                cos_sin[ind_ply_1:ind_ply_2 + 1, :] = np.vstack((
                    cos_sin[ind_ply_2, :], cos_sin[ind_ply_1:ind_ply_2, :]))
#        print('for ind_ply_1')
#        print('best_objD', best_objD)
#        print_ss(ss)
#        print()
        ind_step += 1

#    print()
#    print_ss(ss_ini[:ss_ini.size//2])
#    print_ss(ss[:ss_ini.size//2])
#    print('objD', objD)

    if constraints.sym:
        ss[ss.size //2 + ss.size % 2:] = np.flip(ss[:ss.size //2])

    if count_obj:
        return ss, n_obj_func_D_calls
    return ss

if __name__ == "__main__":

    print('\n*** Test for the function repair_flexural_4 ***')
    constraints = Constraints(sym=True, dam_tol=True)
    parameters = Parameters(
        repair_flexural_switch=True, n_D2=2, constraints=constraints)
    ss = np.array([
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0,
        0, 0, 0, 0, 45, 45, 45, 45, 45, 0, 0, 0, 0], int)
    if constraints.sym:
        ss = np.hstack((ss, 0, np.flip(ss)))

    ss = '90  90  90  90  0  0  0  90  0  0  90  0  0  45  0  90  0  0  45  0  0  90  0  0  0  0  0  45  0  0  0  0  0  0  0  0  0  0  45  0'
    ss = np.array(ss.split('  ')).astype(int)
    ss = np.hstack((ss, np.flip(ss)))
    print('\nInitial stacking sequence')
    print_ss(ss, 80)
    print(ss.size)
    ss_target = '90  90  90  90  90  90  45  45  45  45  45  45  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  45  45  45  45  45  45  90  90  90  90  90  90'
    ss_target = np.array(ss_target.split('  ')).astype(int)
    lampam_target = calc_lampam(ss_target)
    out_of_plane_coeffs = np.array([1, 1, 1, 1])
    ss = repair_flexural_3(
        ss, out_of_plane_coeffs, lampam_target, constraints, parameters)
    print('\nFinal stacking sequence')
    print_ss(ss, 30)
