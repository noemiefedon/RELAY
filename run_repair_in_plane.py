# -*- coding: utf-8 -*-
"""
scripts to test the repair for balance and 10% rule and the refinement for
membrane properties
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np
import pandas as pd
import time
import sys
sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.materials import Material
from src.parameters import Parameters
from src.repair_10_bal import repair_10_bal
from src.repair_10_bal import calc_mini_10
from src.repair_10_bal import is_equal
from src.repair_membrane import repair_membrane
from src.repair_membrane_1_no_ipo import calc_lampamA_ply_queue
from src.ABD import A_from_lampam, B_from_lampam, D_from_lampam, filter_ABD
from src.excel import autofit_column_widths
from src.excel import delete_file
from src.save_set_up import save_constraints_LAYLA
from src.save_set_up import save_materials
from src.pretty_print import print_lampam, print_ss

#==============================================================================
# Input file
#==============================================================================
guidelines = '10-ipo'
n_plies = 150
fibre_angles = 'trad'
fibre_angles = '3060'
fibre_angles = '15'

file_to_open = '/RELAY/pop/'\
 + fibre_angles + '-' + guidelines + '-' + str(n_plies) + 'plies.xlsx'

result_filename = 'in-plane-repair-' + fibre_angles + '-' + guidelines \
+ '-' + str(n_plies) + 'plies-lp-shift-0.2.xlsx'

delete_file(result_filename)
#==============================================================================
# Material properties
#==============================================================================
data = pd.read_excel(
    file_to_open, sheet_name='Materials', index_col=0, header=None)
data = data.transpose()
E11 = data.loc[1, 'E11']
E22 = data.loc[1, 'E22']
nu12 = data.loc[1, 'nu12']
G12 = data.loc[1, 'G12']
ply_t = data.loc[1, 'ply thickness']
mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12, ply_t=ply_t)
#==============================================================================
# Design & manufacturing constraints
#==============================================================================
data = pd.read_excel(
    file_to_open, sheet_name='Constraints', index_col=0, header=None)
data = data.transpose()
#print(data)

sym = data.loc[1, 'symmetry']
bal = True
ipo = True
oopo = data.loc[1, 'out-of-plane orthotropy']
dam_tol = data.loc[1, 'damage tolerance']
rule_10_percent = True
percent_0 = 10
percent_90 = 10
percent_45 = 0
percent_135 = 0
percent_45_135 = 10
diso = data.loc[1, 'diso']
delta_angle = data.loc[1, 'delta_angle']
contig = data.loc[1, 'contig']
n_contig = data.loc[1, 'n_contig']
set_of_angles = np.array(data.loc[1, 'fibre orientations'].split()).astype(int)

constraints = Constraints(
    sym=sym,
    bal=bal,
    ipo=ipo,
    oopo=oopo,
    dam_tol=dam_tol,
    rule_10_percent=rule_10_percent,
    percent_0=percent_0,
    percent_45=percent_45,
    percent_90=percent_90,
    percent_135=percent_135,
    percent_45_135=percent_45_135,
    diso=diso,
    contig=contig,
    n_contig=n_contig,
    delta_angle=delta_angle,
    set_of_angles=set_of_angles)
#==============================================================================
# Parameters
#==============================================================================
# lamination parameter weightings during membrane property refinement
in_plane_coeffs = np.array([1, 1, 0, 0])
# percentage of laminate thickness for plies that can be modified during
# the refinement of membrane properties
p_A = 80
# number of plies in the last permutation during repair for disorientation
# and/or contiguity
n_D1 = 6
# number of ply shifts tested at each step of the re-designing process during
# refinement of flexural properties
n_D2 = 10
# number of times the algorithms 1 and 2 are repeated during the flexural
# property refinement
n_D3 = 2
# lamination parameter weightings during flexural property refinement
out_of_plane_coeffs = np.array([1, 1, 1, 0])

table_param = pd.DataFrame()
table_param.loc[0, 'in_plane_coeffs'] \
= ' '.join(np.array(in_plane_coeffs, dtype=str))
table_param.loc[0, 'out_of_plane_coeffs'] \
= ' '.join(np.array(out_of_plane_coeffs, dtype=str))
table_param.loc[0, 'p_A'] = p_A
table_param.loc[0, 'n_D1'] = n_D1
table_param.loc[0, 'n_D2'] = n_D2
table_param.loc[0, 'n_D3'] = n_D3
table_param = table_param.transpose()

parameters = Parameters(
    constraints=constraints,
    p_A=p_A,
    n_D1=n_D1,
    n_D2=n_D2,
    n_D3=n_D3,
    repair_membrane_switch=True,
    repair_flexural_switch=True)
#==============================================================================
# Tests
#==============================================================================
table_10_bal = pd.DataFrame()
table_membrane = pd.DataFrame()
data = pd.read_excel(file_to_open, sheet_name='stacks', index_col=0)
#print(data)

t_cummul_10_bal = 0
t_cummul_membrane = 0
table_10_bal.loc[0, 'average time repair-10-bal (s)'] = 0
table_membrane.loc[0, 'average time repair-membrane (s)'] = 0

for ind in range(0, 50):
    print('ind', ind)

    #==========================================================================
    # Read inputs
    #==========================================================================
    n_plies = data.loc[ind, 'ply_count']

    lampam_ini = np.empty((12,), float)
    lampam_ini[0] = data.loc[ind, 'lampam[1]']
    lampam_ini[1] = data.loc[ind, 'lampam[2]']
    lampam_ini[2] = data.loc[ind, 'lampam[3]']
    lampam_ini[3] = data.loc[ind, 'lampam[4]']
    lampam_ini[4] = data.loc[ind, 'lampam[5]']
    lampam_ini[5] = data.loc[ind, 'lampam[6]']
    lampam_ini[6] = data.loc[ind, 'lampam[7]']
    lampam_ini[7] = data.loc[ind, 'lampam[8]']
    lampam_ini[8] = data.loc[ind, 'lampam[9]']
    lampam_ini[9] = data.loc[ind, 'lampam[10]']
    lampam_ini[10] = data.loc[ind, 'lampam[11]']
    lampam_ini[11] = data.loc[ind, 'lampam[12]']

    lampam_target = lampam_ini + 0.2

    ss_ini = np.array(data.loc[ind, 'ss'].split()).astype(int)
#    print('ss_ini')
#    print_ss(ss_ini, 200)

    A11_ini = data.loc[ind, 'A11']
    A22_ini = data.loc[ind, 'A22']
    A12_ini = data.loc[ind, 'A12']
    A66_ini = data.loc[ind, 'A66']
    A16_ini = data.loc[ind, 'A16']
    A26_ini = data.loc[ind, 'A26']

    #==========================================================================
    # Repair for balance and 10% rule
    #==========================================================================
    t = time.time()

    mini_10 = calc_mini_10(constraints, ss_ini.size)
    ss, ply_queue = repair_10_bal(ss_ini, mini_10, constraints)
#    print('ss after repair balance/10')
#    print_ss(ss)
#    print(ply_queue)

    elapsed_10_bal = time.time() - t
    t_cummul_10_bal += elapsed_10_bal

    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)

    table_10_bal.loc[ind, 'ply_count'] = n_plies

    table_10_bal.loc[ind, 'time (s)'] = elapsed_10_bal

    table_10_bal.loc[ind, 'no change in ss'] = is_equal(
        ss, ply_queue, ss_ini, constraints.sym)

    table_10_bal.loc[ind, 'f_A ini'] = sum(
        in_plane_coeffs * ((lampam_ini[0:4] - lampam_target[0:4]) ** 2))
    table_10_bal.loc[ind, 'f_A solution'] = sum(
        in_plane_coeffs * ((lampamA - lampam_target[0:4]) ** 2))

    table_10_bal.loc[ind, 'diff lampam 1'] = abs(lampam_ini[0]-lampamA[0])
    table_10_bal.loc[ind, 'diff lampam 2'] = abs(lampam_ini[1]-lampamA[1])
    table_10_bal.loc[ind, 'diff lampam 3'] = abs(lampam_ini[2]-lampamA[2])
    table_10_bal.loc[ind, 'diff lampam 4'] = abs(lampam_ini[3]-lampamA[3])

    table_10_bal.loc[ind, 'lampam[1]'] = lampamA[0]
    table_10_bal.loc[ind, 'lampam[2]'] = lampamA[1]
    table_10_bal.loc[ind, 'lampam[3]'] = lampamA[2]
    table_10_bal.loc[ind, 'lampam[4]'] = lampamA[3]

    table_10_bal.loc[ind, 'lampam_ini[1]'] = lampam_ini[0]
    table_10_bal.loc[ind, 'lampam_ini[2]'] = lampam_ini[1]
    table_10_bal.loc[ind, 'lampam_ini[3]'] = lampam_ini[2]
    table_10_bal.loc[ind, 'lampam_ini[4]'] = lampam_ini[3]

    ss_flatten = ' '.join(np.array(ss, dtype=str))
    table_10_bal.loc[ind, 'ss'] = ss_flatten

    ply_queue_flatten = ' '.join(np.array(ply_queue, dtype=str))
    table_10_bal.loc[ind, 'ply_queue'] = ply_queue_flatten

    ss_ini_flatten = ' '.join(np.array(ss_ini, dtype=str))
    table_10_bal.loc[ind, 'ss_ini'] = ss_ini_flatten

    A = A_from_lampam(lampamA, mat)
    filter_ABD(A=A)
    A11 = A[0, 0]
    A22 = A[1, 1]
    A12 = A[0, 1]
    A66 = A[2, 2]
    A16 = A[0, 2]
    A26 = A[1, 2]

    if A11_ini:
        table_10_bal.loc[ind, 'diff A11 percentage'] \
        = 100 * abs((A11 - A11_ini)/A11_ini)
    else:
        table_10_bal.loc[ind, 'diff A11 percentage'] = 0
    if A22_ini:
        table_10_bal.loc[ind, 'diff A22 percentage'] \
        = 100 * abs((A22 - A22_ini)/A22_ini)
    else:
        table_10_bal.loc[ind, 'diff A22 percentage'] = 0
    if A12_ini:
        table_10_bal.loc[ind, 'diff A12 percentage'] \
        = 100 * abs((A12 - A12_ini)/A12_ini)
    else:
        table_10_bal.loc[ind, 'diff A12 percentage'] = 0
    if A66_ini:
        table_10_bal.loc[ind, 'diff A66 percentage'] \
        = 100 * abs((A66 - A66_ini)/A66_ini)
    else:
        table_10_bal.loc[ind, 'diff A66 percentage'] = 0
    if A16_ini:
        table_10_bal.loc[ind, 'diff A16 percentage'] \
        = 100 * abs((A16 - A16_ini)/A16_ini)
    else:
        table_10_bal.loc[ind, 'diff A16 percentage'] = 0
    if A26_ini:
        table_10_bal.loc[ind, 'diff A26 percentage'] \
        = 100 * abs((A26 - A26_ini)/A26_ini)
    else:
        table_10_bal.loc[ind, 'diff A26 percentage'] = 0

    #==========================================================================
    # Refinement for membrane properties
    #==========================================================================
    t = time.time()

    ss_list, ply_queue_list, lampamA2_list = repair_membrane(
        ss=ss,
        ply_queue=ply_queue,
        mini_10=mini_10,
        in_plane_coeffs=in_plane_coeffs,
        parameters=parameters,
        constraints=constraints,
        lampam_target=lampam_target)

    ss2 = ss_list[0]
    ply_queue2 = ply_queue_list[0]
    lampamA2 = lampamA2_list[0]
#    print('ss after repair membrane')
#    print_ss(ss2)
#    print(ply_queue2)

    elapsed_membrane = time.time() - t
    t_cummul_membrane += elapsed_membrane

    lampamA2_check = calc_lampamA_ply_queue(
        ss2, n_plies, ply_queue2, constraints)

    if not (abs(lampamA2_check - lampamA2) < 1e-10).all():
        raise Exception('This should not happen')

    table_membrane.loc[ind, 'ply_count'] = n_plies

    table_membrane.loc[ind, 'time (s)'] = elapsed_membrane

    table_membrane.loc[ind, 'no change in ss'] \
    = (abs(lampamA - lampamA2) < 1e-10).all()

    f_A_ini = sum(in_plane_coeffs * ((lampamA - lampam_target[0:4]) ** 2))
    table_membrane.loc[ind, 'f_A ini'] = f_A_ini
    f_A_sol = sum(in_plane_coeffs * ((lampamA2 - lampam_target[0:4]) ** 2))
    table_membrane.loc[ind, 'f_A solution'] = f_A_sol

    table_membrane.loc[ind, 'diff lampam 1 solution'] = abs(
        lampamA2[0]-lampam_target[0])
    table_membrane.loc[ind, 'diff lampam 2 solution'] = abs(
        lampamA2[1]-lampam_target[1])
    table_membrane.loc[ind, 'diff lampam 3 solution'] = abs(
        lampamA2[2]-lampam_target[2])
    table_membrane.loc[ind, 'diff lampam 4 solution'] = abs(
        lampamA2[3]-lampam_target[3])

    table_membrane.loc[ind, 'diff lampam 1 before'] = abs(
        lampam_target[0]-lampamA[0])
    table_membrane.loc[ind, 'diff lampam 2 before'] = abs(
        lampam_target[1]-lampamA[1])
    table_membrane.loc[ind, 'diff lampam 3 before'] = abs(
        lampam_target[2]-lampamA[2])
    table_membrane.loc[ind, 'diff lampam 4 before'] = abs(
        lampam_target[3]-lampamA[3])

    table_membrane.loc[ind, 'lampam[1]'] = lampamA2[0]
    table_membrane.loc[ind, 'lampam[2]'] = lampamA2[1]
    table_membrane.loc[ind, 'lampam[3]'] = lampamA2[2]
    table_membrane.loc[ind, 'lampam[4]'] = lampamA2[3]

    table_membrane.loc[ind, 'lampam_target[1]'] = lampam_target[0]
    table_membrane.loc[ind, 'lampam_target[2]'] = lampam_target[1]
    table_membrane.loc[ind, 'lampam_target[3]'] = lampam_target[2]
    table_membrane.loc[ind, 'lampam_target[4]'] = lampam_target[3]

    table_membrane.loc[ind, 'lampam_before[1]'] = lampamA[0]
    table_membrane.loc[ind, 'lampam_before[2]'] = lampamA[1]
    table_membrane.loc[ind, 'lampam_before[3]'] = lampamA[2]
    table_membrane.loc[ind, 'lampam_before[4]'] = lampamA[3]

    ss_flatten = ' '.join(np.array(ss2, dtype=str))
    table_membrane.loc[ind, 'ss'] = ss_flatten

    ply_queue_flatten = ' '.join(np.array(ply_queue2, dtype=str))
    table_membrane.loc[ind, 'ply_queue'] = ply_queue_flatten

    ss_flatten = ' '.join(np.array(ss, dtype=str))
    table_membrane.loc[ind, 'ss_ini'] = ss_flatten

    ply_queue_flatten = ' '.join(np.array(ply_queue, dtype=str))
    table_membrane.loc[ind, 'ply_queue_ini'] = ply_queue_flatten

    A2 = A_from_lampam(lampamA2, mat)
    filter_ABD(A=A)
    A11_2 = A2[0, 0]
    A22_2 = A2[1, 1]
    A12_2 = A2[0, 1]
    A66_2 = A2[2, 2]
    A16_2 = A2[0, 2]
    A26_2 = A2[1, 2]

    if A11_ini:
        table_membrane.loc[ind, 'diff A11 percentage'] \
        = 100 * abs((A11_2 - A11_ini)/A11_ini)
    else:
        table_membrane.loc[ind, 'diff A11 percentage'] = 0
    if A22_ini:
        table_membrane.loc[ind, 'diff A22 percentage'] \
        = 100 * abs((A22_2 - A22_ini)/A22_ini)
    else:
        table_membrane.loc[ind, 'diff A22 percentage'] = 0
    if A12_ini:
        table_membrane.loc[ind, 'diff A12 percentage'] \
        = 100 * abs((A12_2 - A12_ini)/A12_ini)
    else:
        table_membrane.loc[ind, 'diff A12 percentage'] = 0
    if A66_ini:
        table_membrane.loc[ind, 'diff A66 percentage'] \
        = 100 * abs((A66_2 - A66_ini)/A66_ini)
    else:
        table_membrane.loc[ind, 'diff A66 percentage'] = 0
    if A16_ini:
        table_membrane.loc[ind, 'diff A16 percentage'] \
        = 100 * abs((A16_2 - A16_ini)/A16_ini)
    else:
        table_membrane.loc[ind, 'diff A16 percentage'] = 0
    if A26_ini:
        table_membrane.loc[ind, 'diff A26 percentage'] \
        = 100 * abs((A26_2 - A26_ini)/A26_ini)
    else:
        table_membrane.loc[ind, 'diff A26 percentage'] = 0


table_10_bal.loc[0, 'average time repair-10-bal (s)'] = t_cummul_10_bal / 50
table_membrane.loc[0, 'average time repair-membrane (s)'] \
= t_cummul_membrane / 50

writer = pd.ExcelWriter(result_filename)
table_10_bal.to_excel(writer, 'Results-10-bal')
table_membrane.to_excel(writer, 'Results-membrane')
table_param.to_excel(writer, 'Parameters', index=True, header=False)
writer.save()
save_constraints_LAYLA(result_filename, constraints)
save_materials(result_filename, mat)
autofit_column_widths(result_filename)
