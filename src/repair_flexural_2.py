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
from src.disorientation import is_diso
from src.contiguity import is_contig
from src.lampam_functions import calc_lampam

def repair_flexural_2(
        ss_ini, out_of_plane_coeffs, lampam_target, constraints,
        parameters, count_obj=False):
    """
    repair for flexural properties only accounting for one panel by
    changing ply block sizes to improve the out-of-plane lamination
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
#    print('best_objD', best_objD)
#    print('sec_mom_areas', sec_mom_areas, sum(sec_mom_areas))
#    print('cos_sin', cos_sin)
#    print('lampamD', lampamD)
#    print_ss(ss)
#    print('objD', objD)

    ind_step = 0
    for index, ind_ply_1 in enumerate(ind_2[:-1]):
#        print('ind_ply_1', ind_ply_1)

        objD = 1e10 * np.ones((parameters.n_D2 - 1,), float)

        for ind_loc, ind_ply_2 \
        in enumerate(ind_2[index + 1:index + parameters.n_D2]):

#            print('    ind_ply_2', ind_ply_2)

            if not constraints.sym and ind_step % 2 == 0:

                ss_mod = np.copy(ss)
                ss_mod[ind_ply_2:ind_ply_1 + 1] = np.hstack((
                    ss_mod[ind_ply_2 + 1:ind_ply_1 + 1],
                    ss_mod[ind_ply_2]))
#                print_ss(ss_mod)

                if constraints.diso:
                    if ind_ply_2 > 0 and not is_diso(
                            ss_mod[ind_ply_2], ss_mod[ind_ply_2 - 1],
                            constraints.delta_angle):
                        continue
                    if not is_diso(
                            ss_mod[ind_ply_2], ss_mod[ind_ply_2 + 1],
                            constraints.delta_angle):
                        continue
                    if not is_diso(
                            ss_mod[ind_ply_1], ss_mod[ind_ply_1 - 1],
                            constraints.delta_angle):
                        continue
                    if ind_ply_1 < ss.size - 1 \
                    and not is_diso(
                            ss_mod[ind_ply_1], ss_mod[ind_ply_1 + 1],
                            constraints.delta_angle):
                        continue

                if constraints.contig:
                    mini = max(0, ind_ply_1 - constraints.n_contig)
                    maxi = ind_ply_1 + constraints.n_contig
                    if not is_contig(ss_mod[mini:maxi + 1],
                                     constraints.n_contig):
                        continue
                    mini = max(0, ind_ply_2 - constraints.n_contig)
                    maxi = ind_ply_2 + constraints.n_contig
                    if not is_contig(ss_mod[mini:maxi + 1],
                                     constraints.n_contig):
                        continue


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
                if constraints.sym:
                    ss_mod = np.copy(ss[:ss.size//2])
                    ss_mod[ind_ply_1:ind_ply_2 + 1] = np.hstack((
                        ss_mod[ind_ply_2], ss_mod[ind_ply_1:ind_ply_2]))
#                    print_ss(ss_mod)

                    if constraints.diso:
                        if ind_ply_1 > 0 and not is_diso(
                                ss_mod[ind_ply_1], ss_mod[ind_ply_1 - 1],
                                constraints.delta_angle):
                            continue
                        if not is_diso(
                                ss_mod[ind_ply_1], ss_mod[ind_ply_1 + 1],
                                constraints.delta_angle):
                            continue
                        if not is_diso(
                                ss_mod[ind_ply_2], ss_mod[ind_ply_2 - 1],
                                constraints.delta_angle):
                            continue
                        if ss.size % 2 == 1 and ind_ply_2 == ss.size//2 - 1 \
                        and not is_diso(
                                ss_mod[ind_ply_2], ss[ind_ply_2 + 1],
                                constraints.delta_angle):
                            continue
                        if ind_ply_2 < ss.size//2 - 1\
                        and not is_diso(
                                ss_mod[ind_ply_2], ss_mod[ind_ply_2 + 1],
                                constraints.delta_angle):
                            continue

#                    print('diso ok')

                    if constraints.contig:
                        mini = max(0, ind_ply_1 - constraints.n_contig)
                        maxi = ind_ply_1 + constraints.n_contig
#                        print('mini', mini, 'maxi', maxi)

                        if maxi >= ss.size // 2 - 1:
                            if ss.size % 2 == 1:
                                ss_mod = np.hstack((
                                    ss_mod, ss[ss.size//2], np.flip(ss_mod)))
                            else:
                                ss_mod = np.hstack((ss_mod, np.flip(ss_mod)))
#                            print('SS_mod here', ss_mod)
                            maxi = ss.size - mini - 1
                            if not is_contig(ss_mod[mini:maxi + 1],
                                         constraints.n_contig):
                                continue
                        else:
                            if not is_contig(ss_mod[mini:maxi + 1],
                                         constraints.n_contig):
                                continue

                            mini = max(0, ind_ply_2 - constraints.n_contig)
                            maxi = ind_ply_2 + constraints.n_contig
                            if maxi >= ss.size // 2 -1:
                                if ss.size % 2 == 1:
                                    ss_mod = np.hstack((
                                        ss_mod, ss[ss.size//2],
                                        np.flip(ss_mod)))
                                else:
                                    ss_mod = np.hstack((
                                        ss_mod, np.flip(ss_mod)))
#                                print('SS_mod there', ss_mod)
                                maxi = ss.size - mini - 1
                            if not is_contig(ss_mod[mini:maxi + 1],
                                             constraints.n_contig):
                                continue
#                    print('contig ok')

                else:
                    ss_mod = np.copy(ss)
                    ss_mod[ind_ply_1:ind_ply_2 + 1] = np.hstack((
                        ss_mod[ind_ply_2], ss_mod[ind_ply_1:ind_ply_2]))
#                    print_ss(ss_mod)

                    if constraints.diso:
                        if ind_ply_1 > 0 and not is_diso(
                                ss_mod[ind_ply_1], ss_mod[ind_ply_1 - 1],
                                constraints.delta_angle):
                            continue
                        if not is_diso(
                                ss_mod[ind_ply_1], ss_mod[ind_ply_1 + 1],
                                constraints.delta_angle):
                            continue
                        if not is_diso(
                                ss_mod[ind_ply_2], ss_mod[ind_ply_2 - 1],
                                constraints.delta_angle):
                            continue
                        if ind_ply_2 < ss.size - 1\
                        and not is_diso(
                                ss_mod[ind_ply_2], ss_mod[ind_ply_2 + 1],
                                constraints.delta_angle):
                            continue

#                    print('diso ok')

                    if constraints.contig:
                        mini = max(0, ind_ply_1 - constraints.n_contig)
                        maxi = ind_ply_1 + constraints.n_contig
                        if not is_contig(ss_mod[mini:maxi + 1],
                                         constraints.n_contig):
                            continue
                        mini = max(0, ind_ply_2 - constraints.n_contig)
                        maxi = ind_ply_2 + constraints.n_contig
                        if not is_contig(ss_mod[mini:maxi + 1],
                                         constraints.n_contig):
                            continue

#                    print('contig ok')

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

        if min(objD) + 1e-20 < best_objD:

            best_objD = min(objD)
            ind_min = np.argmin(objD)
            ind_ply_2 = ind_2[index + 1:][ind_min]

#            print('CHANGE')
#            print('objD', objD)
#            print('ind_min', ind_min, 'ind_ply_2', ind_ply_2)

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

#        print()
#        print('after iteration')
#        print_ss(ss)
#        print('best_objD', best_objD)
#        print('objD', objD)

        ind_step += 1

#    print_ss(ss)
#    print('objD', objD)

    if constraints.sym:
        ss[ss.size // 2 + ss.size % 2:] = np.flip(ss[:ss.size // 2])

    if count_obj:
        return ss, n_obj_func_D_calls
    return ss

if __name__ == "__main__":

    print('\n*** Test for the function repair_flexural_2 ***')
    constraints = Constraints(
        sym=True,
        dam_tol=False,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=4,)
    parameters = Parameters(
        repair_flexural_switch=True, n_D2=3, constraints=constraints)
    ss = np.array([45, 0, 45, 90, -45, -45, 0], int)
    if constraints.sym:
        ss = np.hstack((ss, np.flip(ss)))
    print('\nInitial stacking sequence')
    print_ss(ss, 80)
    print(ss.size)
    ss_target = 0*np.ones((1,), dtype=int)
    lampam_target = calc_lampam(ss_target)
    out_of_plane_coeffs = np.array([1/3, 1/3, 1/3, 0])
    ss = repair_flexural_2(
        ss, out_of_plane_coeffs, lampam_target, constraints, parameters)
    print_ss(ss, 30)
