# -*- coding: utf-8 -*-
"""
repair for disorientation and contiguity

- repair_diso_contig_list
    successive attempts at at repairing stacking sequences from a list for the
    disorientation and contiguity rules until a successful repair or exhaustion
    of the list

- repair_diso_contig
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_ss
from src.repair_from_inside_sym import repair_diso_contig_from_inside_sym
from src.repair_from_inside_asym import repair_diso_contig_from_inside_asym
from src.repair_from_outside_sym import repair_diso_contig_from_outside_sym
from src.repair_from_outside_asym import repair_diso_contig_from_outside_asym

def repair_diso_contig_list(
        ss_list, ply_queue_list, constraints, n_D1):
    """
    successive attempts at at repairing stacking sequences from a list for the
    disorientation and contiguity rules until a successful repair or exhaustion
    of the list

    - n_D1: number of plies in the last permutation
    """

    for ind, (ss, ply_queue) \
    in enumerate(zip(ss_list, ply_queue_list)):
        ss, completed_inward, completed_outward = repair_diso_contig(
            ss, ply_queue, constraints, n_D1)
        if completed_inward or completed_outward:
            return ss, completed_inward, completed_outward, ind
        if ind == 0:
            ss0 = ss

    return ss0, False, False, 0


def repair_diso_contig(
        ss, ply_queue, constraints, n_D1):
    """
    attempts at repairing the stacking sequence to satisfy the disorientation
    and contiguity rule

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    - from_inside: True if the swaps are performed from inside out of the
    laminate, shifting the modifications outside of the laminate
    - n_D1: number of plies in the last permutation
    """
#    print_ss(ss)
#    print('ply_queue', ply_queue)
#    print('constraints.sym', constraints.sym)
#    print('n_D1', n_D1)

    if constraints.sym:
        ss, inward = repair_diso_contig_from_outside_sym(
            ss, ply_queue, constraints, n_D1)
        if not inward:
            ss, outward = repair_diso_contig_from_inside_sym(
                ss, ply_queue, constraints, n_D1)
            return ss, inward, outward
        return ss, inward, True

    ss, inward = repair_diso_contig_from_outside_asym(
        ss, ply_queue, constraints, n_D1)
    if not inward:
        ss, outward = repair_diso_contig_from_inside_asym(
            ss, ply_queue, constraints, n_D1)
        return ss, inward, outward
    return ss, inward, True


if __name__ == "__main__":
    print('\n*** Test for the function repair_diso_contig ***')
    constraints = Constraints(
        sym=True,
        ipo=False,
        dam_tol=False,
        rule_10_percent=False,
        diso=True,
        contig=True,
        delta_angle=45,
        n_contig=4,
        set_of_angles=[0, 45, -45, 90])
    ss_ini = np.array([0, 45, -45, 666, 666, 666, 666], int)
    ss_ini = np.hstack((ss_ini, np.flip(ss_ini)))
    print_ss(ss_ini, 40)
    print('ss_ini.size', ss_ini.size)
    n_D1 = 6
    ply_queue = [-45, 90, 90, 90]
    ss, inward, outward = repair_diso_contig(
        ss_ini, ply_queue, constraints, n_D1)
    print('Inward repair successful?', inward)
    print('Outward repair successful?', outward)
    print_ss(ss, 40)
