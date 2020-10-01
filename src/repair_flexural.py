# -*- coding: utf-8 -*-
"""
repair for flexural properties

- repair_membrane_1:
    repair for membrane properties only accounting for one panel

- repair_membrane_2:
    repair for membrane properties accounting for all the panels
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.repair_flexural_1 import repair_flexural_1
from src.repair_flexural_2 import repair_flexural_2
from src.repair_flexural_3 import repair_flexural_3
from src.disorientation import is_diso_ss
from src.contiguity import is_contig
from src.constraints import Constraints
from src.parameters import Parameters
from src.pretty_print import print_lampam, print_ss, print_list_ss
from src.lampam_functions import calc_lampam

def repair_flexural(
        ss, constraints,
        out_of_plane_coeffs,
        parameters,
        multipanel=None, lampam_target=None,
        reduced_pdl=None, reduced_sst=None, reduced_lampam=None,
        reduced_ss=None, mat=None, repair_flexural_for_ref_panel_only=True,
        count_obj=False):
    """
    repair for flexural properties

    - count_obj: flag to count the number of objective function calls
    """
    n_obj_func_D_calls = 0

    if multipanel is not None and not repair_flexural_for_ref_panel_only:
        raise('Not coded')

    if multipanel is not None:
        if parameters.repair_flexural_switch \
        and not np.isclose(
            np.array([0, 0, 0, 0], float), out_of_plane_coeffs).all():

            if constraints.diso or constraints.contig:
                for ite in range(parameters.n_D3):
                    # swap blocks of plies
                    ss, n_obj = repair_flexural_1(
                        ss, out_of_plane_coeffs, lampam_target, constraints,
                        count_obj=True)
                    n_obj_func_D_calls += n_obj
                    # changes ply block sizes
                    ss, n_obj = repair_flexural_2(
                        ss, out_of_plane_coeffs, lampam_target, constraints,
                        parameters, count_obj=True)
                    n_obj_func_D_calls += n_obj
            else:
                for ite in range(parameters.n_D3):
                    # re-design laminate
                    ss, n_obj = repair_flexural_3(
                        ss, out_of_plane_coeffs, lampam_target, constraints,
                        parameters, count_obj=True)
                    n_obj_func_D_calls += n_obj
        else:
            reduced_sst = convert_ss_ref_to_reduced_sst(
                ss, reduced_pdl=reduced_pdl, ind_ref=multipanel.ind_ref,
                reduced_ss_before=reduced_ss)
        if count_obj:
            return ss, reduced_sst, n_obj_func_D_calls
        return ss, reduced_sst


    if parameters.repair_flexural_switch \
    and not np.isclose(out_of_plane_coeffs,
                       np.array([0, 0, 0, 0], float)).all():
        if constraints.diso or constraints.contig:

            for ite in range(parameters.n_D3):
                # swap blocks of plies
                ss, n_obj = repair_flexural_1(
                    ss, out_of_plane_coeffs, lampam_target, constraints,
                    count_obj=True)
                n_obj_func_D_calls += n_obj

                if hasattr(constraints, 'dam_tol_rule'):
                    if constraints.diso and not is_diso_ss(
                            ss,
                            delta_angle=constraints.delta_angle,
                            dam_tol=constraints.dam_tol,
                            dam_tol_rule=constraints.dam_tol_rule):
                        raise Exception('Disorientation rule not satisfied')
                else:
                    if constraints.diso and not is_diso_ss(
                            ss,
                            delta_angle=constraints.delta_angle,
                            dam_tol=constraints.dam_tol,
                            n_plies_dam_tol=constraints.n_plies_dam_tol):
                        raise Exception('Disorientation rule not satisfied')

                if constraints.contig and not is_contig(
                        ss, constraints.n_contig):
                    raise Exception('Contiguity rule not satisfie')

                # changes ply block sizes
                ss, n_obj = repair_flexural_2(
                    ss, out_of_plane_coeffs, lampam_target, constraints,
                    parameters, count_obj=True)
                n_obj_func_D_calls += n_obj
        else:
            for ite in range(parameters.n_D3):
                # re-design laminate
                ss, n_obj = repair_flexural_3(
                    ss, out_of_plane_coeffs, lampam_target, constraints,
                    parameters, count_obj=True)
                n_obj_func_D_calls += n_obj

        if count_obj:
            return ss, n_obj_func_D_calls
        return ss

    if count_obj:
        return ss, n_obj_func_D_calls
    return ss

if __name__ == "__main__":

    print('\n*** Test for the function repair_flexural ***')
    constraints = Constraints(
        sym=True,
        bal=True,
        ipo=True,
        dam_tol=False,
        rule_10_percent=True,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=4,
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10,
        set_of_angles=[0, 45, -45, 90, 30, -30, 60, -60, 15, -15, 75, -75])
    parameters = Parameters(
        repair_flexural_switch=True, n_D2=2, constraints=constraints)

    ss = np.array([
            -60, -45,  90, 75, -75, -45,   0], int)
    if constraints.sym:
        ss = np.hstack((ss, np.flip(ss)))
    if not is_diso_ss(ss,
                      delta_angle=constraints.delta_angle,
                      dam_tol=constraints.dam_tol,
                      dam_tol_rule=constraints.dam_tol_rule):
        raise Exception('Disorientation rule not satisfied initially')

    if not is_contig(ss, constraints.n_contig):
        raise Exception('Contiguity rule not satisfied initially')


    print('\nInitial stacking sequence')
    print_ss(ss, 1000)
    ss_target = 45*np.ones((1,), dtype=int)
    lampam_target = calc_lampam(ss_target)
    lampam_target = np.array([
        -0.25662112, -0.01727515, -0.73962959, 0.08359081,
        0.43671968, -0.7901057, 0.98404481, 0.65070345,
        0.97056517, 0.7023994, 0.32539113, -0.35851357])
    out_of_plane_coeffs = np.array([1, 1, 1, 1])
    ss = repair_flexural(
        ss, constraints, out_of_plane_coeffs,
        parameters=parameters,
        lampam_target=lampam_target)
    print('\nFinal stacking sequence')
    print_ss(ss, 1000)

    if not is_diso_ss(ss,
                      delta_angle=constraints.delta_angle,
                      dam_tol=constraints.dam_tol,
                      dam_tol_rule=constraints.dam_tol_rule):
        raise Exception('Disorientation rule not satisfied')

    if not is_contig(ss, constraints.n_contig):
        raise Exception('Contiguity rule not satisfied')

