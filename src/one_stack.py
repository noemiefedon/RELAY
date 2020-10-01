# -*- coding: utf-8 -*-
"""
Functions to check a design manufacturability

- check_ss_manufacturability
    checks the manufacturability of a stacking sequence list
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
sys.path.append(r'C:\RELAY')
from src.contiguity import is_contig
from src.disorientation import is_diso_ss
from src.balance import is_balanced
from src.dam_tol import is_dam_tol
from src.ten_percent_rule import is_ten_percent_rule
from src.lp_functions_2 import calc_lampamA
from src.constraints import Constraints
from src.pretty_print import print_ss

def check_ss_manufacturability(
        ss, constraints, no_ipo_check=False, no_bal_check=False,
        equality_45_135=False, equality_0_90=False, n_plies=None):
    """
    checks the manufacturability of a stacking sequence list
    """
    if n_plies is not None and ss.size != n_plies:
        raise Exception("Wrong number of plies")

    if constraints.dam_tol:
        if not is_dam_tol(ss, constraints):
            print_ss(ss)
            raise Exception("Damage tolerance constraint not satisfied")

    if not no_bal_check and constraints.bal:
        if not is_balanced(ss, constraints):
            raise Exception("Balance constraint not satisfied")

    if not no_ipo_check and constraints.ipo:
        lampamA = calc_lampamA(ss, constraints)
        if (abs(lampamA[2:4]) > 1e-10).any():
            print_ss(ss)
            print('lampamA', lampamA)
#            print('ipo')
            raise Exception("In plane orthotropy constraint not satisfied")


    if constraints.diso:
        if hasattr(constraints, 'dam_tol_rule'):
            if not is_diso_ss(ss, constraints.delta_angle,
                              constraints.dam_tol, constraints.dam_tol_rule):
                raise Exception("Disorientation constraint not satisfied")
        else:
            if not is_diso_ss(ss, constraints.delta_angle,
                              constraints.dam_tol, constraints.n_plies_dam_tol):
                raise Exception("Disorientation constraint not satisfied")

    if constraints.contig:
        if not is_contig(ss, constraints.n_contig):
            raise Exception("Contiguity constraint not satisfied")

    if constraints.rule_10_percent:
        if not is_ten_percent_rule(
                constraints, stack=ss,
                equality_45_135=equality_45_135,
                equality_0_90=equality_0_90):
            raise Exception("10% rule not satisfied")

    return 0



if __name__ == "__main__":

    print('\n*** Test for the function check_ss_manufacturability ***')
    constraints = Constraints(
        sym=True,
        bal=True,
        ipo=True,
        oopo=False,
        dam_tol=False,
        rule_10_percent=True,
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10,
        diso=True,
        contig=True,
        n_contig=5,
        delta_angle=45,
        set_of_angles=np.array([0, 45, -45, 90]))
    ss = np.array([  0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0,  45,   0,   0,   0,  45,   0,   0,   0,   0,  45,  90,  90,  90,  90,  90, -45,  90,  90,  90,  90,  45,   0,   0,   0,   0,   0,  45,   0,   0,   0,   0, -45,  90,  90,  90,  90,  90, -45,  90,  90,  90,  90,  90,  45,   0,   0,  45,  90,  90,  90,  90,  45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0,  45,  90,  90,  90, -45,   0,   0,  45,  45,   0,   0, -45,  90,  90,  90,  45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0,  45,  90,  90,  90,  90,  45,   0,   0,  45,  90,  90,  90,  90,  90, -45,  90,  90,  90,  90,  90, -45,   0,   0,   0,   0,  45,   0,   0,   0,   0,   0,  45,  90,  90,  90,  90, -45,  90,  90,  90,  90,  90,  45,   0,   0,   0,   0,  45,   0,   0,   0,  45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0, -45,   0,   0,   0,   0,   0], int)
    check_ss_manufacturability(ss, constraints)
