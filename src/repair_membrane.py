# -*- coding: utf-8 -*-
"""
repair for membrane properties

- repair_membrane
    repairs a laminate to improve its in-plane stiffness properties

- repair_membrane_1:
    repair for membrane properties only accounting for one panel

- repair_membrane_2:
    repair for membrane properties accounting for all the panels
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.repair_membrane_1_ipo import repair_membrane_1_ipo
from src.repair_membrane_1_no_ipo import repair_membrane_1_no_ipo
from src.repair_membrane_1_no_ipo import calc_lampamA_ply_queue

def repair_membrane(
        ss, ply_queue, mini_10, in_plane_coeffs, constraints, parameters,
        multipanel=None, lampam_target=None,
        reduced_ss=None, reduced_pdl=None, mat=None):
    """
    repairs a laminate to improve its in-plane stiffness properties
    """
    if not multipanel is None \
    and not parameters.repair_membrane_for_ref_panel_only:
        if parameters.opti_problem != 'Lamination parameters matching' \
        or not parameters.repair_membrane_switch \
        or np.isclose(np.array([0, 0, 0, 0], float), in_plane_coeffs).all():
            ss_list = [ss] # no in-plane optimisation required
            ply_queue_list = [ply_queue]
            lampamA_list = [
                calc_lampamA_ply_queue(ss, ss.size, ply_queue, constraints)]
        else:
            raise Exception('Not coded')
        return ss_list, ply_queue_list, lampamA_list

    if not multipanel is None:
        if parameters.opti_problem != 'Lamination parameters matching' \
        or not parameters.repair_membrane_switch \
        or np.isclose(np.array([0, 0, 0, 0], float), in_plane_coeffs).all():
            ss_list = [ss] # no in-plane optimisation required
            ply_queue_list = [ply_queue]
            lampamA_list = [
                calc_lampamA_ply_queue(ss, ss.size, ply_queue, constraints)]
        else:
            ss_list, ply_queue_list, lampamA_list, _ = repair_membrane_1(
                ss, ply_queue, mini_10,
                in_plane_coeffs, parameters.p_A,
                lampam_target, constraints)
        return ss_list, ply_queue_list, lampamA_list

    if not parameters.repair_membrane_switch \
    or np.isclose(np.array([0, 0, 0, 0], float), in_plane_coeffs).all():
        ss_list = [ss] # no in-plane optimisation required
        ply_queue_list = [ply_queue]
        lampamA_list = [
            calc_lampamA_ply_queue(ss, ss.size, ply_queue, constraints)]
    else:
        ss_list, ply_queue_list, lampamA_list, _ = repair_membrane_1(
            ss, ply_queue, mini_10,
            in_plane_coeffs, parameters.p_A,
            lampam_target, constraints)
    return ss_list, ply_queue_list, lampamA_list


def repair_membrane_1(
        ss, ply_queue, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints):
    """
    repair for membrane properties only accounting for one panel
    """
    if constraints.ipo:
        return repair_membrane_1_ipo(
            ss, ply_queue, mini_10, in_plane_coeffs,
            p_A, lampam_target, constraints)
    return repair_membrane_1_no_ipo(
        ss, ply_queue, mini_10, in_plane_coeffs,
        p_A, lampam_target, constraints)
