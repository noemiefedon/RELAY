# -*- coding: utf-8 -*-
"""
Function to save laminate design set-up

- save_constraints_LAYLA
    save the design and manufacturing constraints on Sheet [Constraints]

- save_parameters_LAYLA
    saves the optimiser parameters on Sheet [Parameters]

- save_materials
    saves the material properties on Sheet [Materials]
"""
import sys
import numpy as np
import pandas as pd

sys.path.append(r'C:\RELAY')
from src.excel import append_df_to_excel


def save_materials(filename, materials):
    """
    saves the material properties on Sheet [Materials]
    """
    table_mat = pd.DataFrame()
    table_mat.loc[0, 'E11'] = materials.E11
    table_mat.loc[0, 'E22'] = materials.E22
    table_mat.loc[0, 'G12'] = materials.G12
    table_mat.loc[0, 'nu12'] = materials.nu12
    table_mat.loc[0, 'nu21'] = materials.nu21
    table_mat.loc[0, 'areal density'] = materials.density_area
    table_mat.loc[0, 'volumic density'] = materials.density_volume
    table_mat.loc[0, 'ply thickness'] = materials.ply_t

    table_mat.loc[0, 'Q11'] = materials.Q11
    table_mat.loc[0, 'Q12'] = materials.Q12
    table_mat.loc[0, 'Q22'] = materials.Q22
    table_mat.loc[0, 'Q66'] = materials.Q66

    table_mat.loc[0, 'U1'] = materials.U1
    table_mat.loc[0, 'U2'] = materials.U2
    table_mat.loc[0, 'U3'] = materials.U3
    table_mat.loc[0, 'U4'] = materials.U4
    table_mat.loc[0, 'U5'] = materials.U5

    table_mat = table_mat.transpose()

    append_df_to_excel(
        filename, table_mat, 'Materials', index=True, header=False)

def save_constraints_LAYLA(filename, constraints, not_constraints=None):
    """
    saves the design and manufacturing constraints on Sheet [Constraints]
    """
    table_const = pd.DataFrame()
    table_const.loc[0, 'symmetry'] = constraints.sym
    table_const.loc[0, 'balance'] = constraints.bal
    table_const.loc[0, 'in-plane orthotropy'] = constraints.ipo
    table_const.loc[0, 'out-of-plane orthotropy'] = constraints.oopo
    table_const.loc[0, 'damage tolerance'] = constraints.dam_tol
    if hasattr(constraints, 'dam_tol_rule'):
        table_const.loc[0, 'dam_tol_rule'] = constraints.dam_tol_rule
    else:
        table_const.loc[0, 'n_plies_dam_tol'] = constraints.n_plies_dam_tol
    table_const.loc[0, '10% rule'] = constraints.rule_10_percent
    if constraints.rule_10_percent:
        table_const.loc[0, 'percent_0'] = constraints.percent_0 * 100
        table_const.loc[0, 'percent_45'] = constraints.percent_45 * 100
        table_const.loc[0, 'percent_90'] = constraints.percent_90 * 100
        table_const.loc[0, 'percent_-45'] = constraints.percent_135 * 100
        table_const.loc[0, 'percent_+-45'] = constraints.percent_45_135 * 100
    else:
        table_const.loc[0, 'percent_0'] = 0
        table_const.loc[0, 'percent_45'] = 0
        table_const.loc[0, 'percent_90'] = 0
        table_const.loc[0, 'percent_-45'] = 0
        table_const.loc[0, 'percent_+-45'] = 0
    table_const.loc[0, 'diso'] = constraints.diso
    table_const.loc[0, 'delta_angle'] = constraints.delta_angle
    table_const.loc[0, 'contig'] = constraints.contig
    table_const.loc[0, 'n_contig'] = constraints.n_contig_c
    sets = np.array(constraints.set_of_angles, dtype=str)
    table_const.loc[0, 'fibre orientations'] = ' '.join(sets)

    if not_constraints:
        if not_constraints.diso:
            table_const.loc[0, 'No diso !'] = not_constraints.diso
            table_const.loc[0, 'delta_angle (no)'] \
            = not_constraints.delta_angle
        if not_constraints.contig:
            table_const.loc[0, 'No contig !'] = not_constraints.contig
            table_const.loc[0, 'n_contig (no)'] = not_constraints.n_contig_c
        if not_constraints.bal:
            table_const.loc[0, 'No balance !'] \
            = not_constraints.ipo
        if not_constraints.ipo:
            table_const.loc[0, 'No in-plane orthotropy !'] \
            = not_constraints.ipo
        if not_constraints.rule_10_percent:
            table_const.loc[0, 'No 10% rule !'] \
            = not_constraints.rule_10_percent
            table_const.loc[0, 'percent_0 (no)'] \
            = 100 * not_constraints.percent_0
            table_const.loc[0, 'percent_45 (no)'] \
            = 100 * not_constraints.percent_45
            table_const.loc[0, 'percent_90 (no)'] \
            = 100 * not_constraints.percent_90
            table_const.loc[0, 'percent_-45 (no)'] \
            = 100 * not_constraints.percent_135
            table_const.loc[0, 'percent_+-45 (no)'] \
            = 100 * not_constraints.percent_45_135

    table_const = table_const.transpose()

    append_df_to_excel(
        filename, table_const, 'Constraints', index=True, header=False)


def save_parameters_LAYLA_V02(filename, parameters):
    """
    saves the optimiser parameters on Sheet [Parameters]
    """
    table_param = pd.DataFrame()

    table_param.loc[0, 'maximum number of iterations'] \
    = parameters.n_outer_step
    table_param.loc[0, 'global_node_limit'] \
    = parameters.global_node_limit
    table_param.loc[0, 'global_node_limit_p'] \
    = parameters.global_node_limit_p
    table_param.loc[0, 'local_node_limit'] \
    = parameters.local_node_limit
    table_param.loc[0, 'repair_membrane_switch'] \
    = parameters.repair_membrane_switch
    table_param.loc[0, 'repair_flexural_switch'] \
    = parameters.repair_flexural_switch
    table_param.loc[0, 'p_A'] \
    = parameters.p_A
    table_param.loc[0, 'n_D1'] \
    = parameters.n_D1
    table_param.loc[0, 'n_D2'] \
    = parameters.n_D2
    table_param.loc[0, 'n_D3'] \
    = parameters.n_D3
    table_param.loc[0, 'penalty for the 10% rule (ply counts)'] \
    = parameters.penalty_10_pc_switch
    table_param.loc[0, 'penalty for the 10% rule (lampam)'] \
    = parameters.penalty_10_lampam_switch
    table_param.loc[0, 'penalty for balance'] \
    = parameters.penalty_bal_switch
    table_param.loc[0, 'penalty for in-plane orthotropy'] \
    = parameters.penalty_ipo_switch
    table_param.loc[0, 'coeff for in-plane orthotropy/balance penalty'] \
    = parameters.coeff_bal_ipo
    table_param.loc[0, 'coeff for out-of plane orthotropy penalty'] \
    = parameters.coeff_oopo
    table_param.loc[0, 'coeff for the 10% rule penalty'] \
    = parameters.coeff_10
    table_param.loc[0, 'type of objective function'] \
    = 'Norm ' + str(parameters.type_obj_func)

    table_param.loc[0, 'minimum group size'] = parameters.group_size_min
    table_param.loc[0, 'maximum group size'] \
    = ' '.join(np.array(parameters.group_size_max, dtype=str))

    table_param.loc[0, 'first_level_sensitivities[1]'] \
    = parameters.first_level_sensitivities[0]
    table_param.loc[0, 'first_level_sensitivities[2]'] \
    = parameters.first_level_sensitivities[1]
    table_param.loc[0, 'first_level_sensitivities[3]'] \
    = parameters.first_level_sensitivities[2]
    table_param.loc[0, 'first_level_sensitivities[4]'] \
    = parameters.first_level_sensitivities[3]
    table_param.loc[0, 'first_level_sensitivities[5]'] \
    = parameters.first_level_sensitivities[4]
    table_param.loc[0, 'first_level_sensitivities[6]'] \
    = parameters.first_level_sensitivities[5]
    table_param.loc[0, 'first_level_sensitivities[7]'] \
    = parameters.first_level_sensitivities[6]
    table_param.loc[0, 'first_level_sensitivities[8]'] \
    = parameters.first_level_sensitivities[7]
    table_param.loc[0, 'first_level_sensitivities[9]'] \
    = parameters.first_level_sensitivities[8]
    table_param.loc[0, 'first_level_sensitivities[10]'] \
    = parameters.first_level_sensitivities[9]
    table_param.loc[0, 'first_level_sensitivities[11]'] \
    = parameters.first_level_sensitivities[10]
    table_param.loc[0, 'first_level_sensitivities[12]'] \
    = parameters.first_level_sensitivities[11]

    table_param.loc[0, 'lampam_to_be_optimised[1]'] \
    = parameters.lampam_to_be_optimised[0]
    table_param.loc[0, 'lampam_to_be_optimised[2]'] \
    = parameters.lampam_to_be_optimised[1]
    table_param.loc[0, 'lampam_to_be_optimised[3]'] \
    = parameters.lampam_to_be_optimised[2]
    table_param.loc[0, 'lampam_to_be_optimised[4]'] \
    = parameters.lampam_to_be_optimised[3]
    table_param.loc[0, 'lampam_to_be_optimised[5]'] \
    = parameters.lampam_to_be_optimised[4]
    table_param.loc[0, 'lampam_to_be_optimised[6]'] \
    = parameters.lampam_to_be_optimised[5]
    table_param.loc[0, 'lampam_to_be_optimised[7]'] \
    = parameters.lampam_to_be_optimised[6]
    table_param.loc[0, 'lampam_to_be_optimised[8]'] \
    = parameters.lampam_to_be_optimised[7]
    table_param.loc[0, 'lampam_to_be_optimised[9]'] \
    = parameters.lampam_to_be_optimised[8]
    table_param.loc[0, 'lampam_to_be_optimised[10]'] \
    = parameters.lampam_to_be_optimised[9]
    table_param.loc[0, 'lampam_to_be_optimised[11]'] \
    = parameters.lampam_to_be_optimised[10]
    table_param.loc[0, 'lampam_to_be_optimised[12]'] \
    = parameters.lampam_to_be_optimised[11]

    table_param = table_param.transpose()

    append_df_to_excel(
        filename, table_param, 'Parameters', index=True, header=False)
