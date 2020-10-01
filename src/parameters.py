# -*- coding: utf-8 -*-
"""
Class for the parameters of the optimiser
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\LAYLA_and_BELLA')
sys.path.append(r'C:\RELAY')
from src.constraints import Constraints

class Parameters(object):
    " An object for storing the optimiser parameters "
    def __init__(self,
                 constraints,
                 n_outer_step=1,
                 coeff_10=0,
                 coeff_bal_ipo=1,
                 coeff_oopo=1,
                 group_size_min=6,
                 group_size_max=np.array([8]),
                 first_level_sensitivities=np.ones((12,), float),
                 lampam_to_be_optimised=np.ones((12,), float),
                 global_branching_limit=100,
                 global_branching_limit_p=10,
                 local_branching_limit=1000,
                 repair_membrane_switch=True,
                 repair_flexural_switch=True,
                 p_A=100,
                 n_D1=6,
                 n_D2=10,
                 n_D3=1,
                 penalty_10_lampam_switch=False,
                 penalty_10_pc_switch=False,
                 penalty_bal_switch=False,
                 penalty_ipo_switch=False,
                 type_obj_func=2):
        " Create a set of parameters for LAYLA"

        self.opti_problem = 'Lamination parameters matching'
        self.type_obj_func = type_obj_func

        # coefficients for the for the in-plane lampams
        self.n_outer_step = np.around(n_outer_step)
        if not(isinstance(n_outer_step, int)):
            raise ParametersDefinitionError(
                    'Attention, n_outer_step must be an integer!')
        if 1 > n_outer_step:
            raise ParametersDefinitionError(
                    'Attention, n_outer_step must be strictly positive!')

        ### 10% rule
        if not constraints.rule_10_percent:
            # penalty for the 10% rule based on ply count restrictions
            self.penalty_10_pc_switch = False
            # penalty for the 10% rule based on lamination parameter restrictions
            self.penalty_10_lampam_switch = False
            # Coefficient for the 10% rule penalty
            self.coeff_10 = 0
        else:
            # penalty for the 10% rule based on ply count restrictions
            self.penalty_10_pc_switch = penalty_10_pc_switch
            if not isinstance(penalty_10_pc_switch, bool):
                raise ParametersDefinitionError("""
Attention, penalty_10_pc_switch should be a boolean value!""")
            # penalty for the 10% rule based on lamination parameter restrictions
            self.penalty_10_lampam_switch = penalty_10_lampam_switch
            if not isinstance(penalty_10_lampam_switch, bool):
                raise ParametersDefinitionError("""
Attention, penalty_10_lampam_switch should be a boolean value!""")
            # Coefficient for the 10% rule penalty
            self.coeff_10 = coeff_10
            if not(isinstance(coeff_10, int) or isinstance(coeff_10, float)):
                raise ParametersDefinitionError("""
Attention, coeff_10 must be a number (float or integer)!""")
            if coeff_10 < 0:
                raise ParametersDefinitionError(
'The weight of penalty for the 10% rule must be a positive!')
            self.penalty_10_pc_switch = penalty_10_pc_switch
        constraints.penalty_10_pc_switch = penalty_10_pc_switch

        if penalty_10_pc_switch and penalty_10_lampam_switch:
            raise ParametersDefinitionError(
'You cannot use both penalties on ply counts and on lampam for the 10% rule!')

        if not constraints.bal:
            # penalty based on ply counts
            self.penalty_bal_switch = False
        else:
            # penalty based on ply counts
            self.penalty_bal_switch = penalty_bal_switch
            if not isinstance(penalty_bal_switch, bool):
                raise ParametersDefinitionError("""
Attention, penalty_bal_switch should be a boolean value!""")

        if not constraints.ipo:
            # penalty based on lamination parameters
            self.penalty_ipo_switch = False
            # coefficient for the penalty
            self.coeff_bal_ipo = 0
        else:
            # penalty based on lamination parameters
            self.penalty_ipo_switch = penalty_ipo_switch
            if not isinstance(penalty_ipo_switch, bool):
                raise ParametersDefinitionError("""
Attention, penalty_ipo_switch should be a boolean value!""")
            # coefficient for the penalty
            self.coeff_bal_ipo = coeff_bal_ipo
            if not(isinstance(coeff_bal_ipo, int) \
                   or isinstance(coeff_bal_ipo, float)):
                raise ParametersDefinitionError("""
Attention, coeff_bal_ipo must be a number (float or integer)!""")
            if coeff_bal_ipo < 0:
                raise ParametersDefinitionError(
'The weight of penalty for in-plane orthotropy must be a positive!')

        if penalty_ipo_switch and penalty_bal_switch:
            raise ParametersDefinitionError(
'You cannot use both the penalties for balance and in-plane orthotropy!')

        ### out-of-plane orthotropy
        if not constraints.oopo:
            # coefficient for the penalty based on lamination parameters
            self.coeff_oopo = 0
        else:
            # coefficient for the penalty based on lamination parameters
            self.coeff_oopo = coeff_oopo
            if not(isinstance(coeff_oopo, int) \
                   or isinstance(coeff_oopo, float)):
                raise ParametersDefinitionError("""
Attention, coeff_oopo must be a number (float or integer)!""")
            if coeff_oopo < 0:
                raise ParametersDefinitionError(
'The weight of penalty for out-of-plane orthotropy must be a positive!')

        ### repair to improve the convergence of in-plane lamination parameters
        # and of out-of-plane lamination parameters
        self.repair_membrane_switch = repair_membrane_switch
        if not isinstance(repair_membrane_switch, bool):
            raise ParametersDefinitionError("""
Attention, repair_membrane_switch should be a boolean value!""")
        self.repair_flexural_switch = repair_flexural_switch
        if not isinstance(repair_flexural_switch, bool):
            raise ParametersDefinitionError("""
Attention, repair_flexural_switch should be a boolean value!""")

        # coefficient for the proportion of the laminate thickness that can be
        # modified during the refinement for membrane properties in the repair
        # process
        self.p_A \
        = p_A
        if not isinstance(p_A, (float, int)):
            raise ParametersDefinitionError("""
Attention, p_A should have a numeric value!""")
        if not (0 <= p_A <= 100):
            raise ParametersDefinitionError("""
Attention, p_A must be between 0 and 100!""")

        # n_D1: number of plies in the last permutation
        # during repair for disorientation and/or contiguity
        self.n_D1 = n_D1
        if not isinstance(n_D1, int):
            raise ParametersDefinitionError("""
Attention, n_D1 must be an integer!""")

        # n_D2: number of ply shifts tested at each step of the
        # re-designing process during refinement of flexural properties
        self.n_D2 = n_D2
        if not isinstance(n_D2, int):
            raise ParametersDefinitionError("""
Attention, n_D2 must be an integer!""")

        # n_D3: number of times the algorithms 1 and 2 are repeated during the
        # flexural property refinement
        self.n_D3 = n_D3
        if not isinstance(n_D3, int):
            raise ParametersDefinitionError("""
Attention, n_D2 must be an integer!""")

        ### size of the groups
        # desired group size for smaller groups
        self.group_size_min = np.around(group_size_min)
        if not(isinstance(group_size_min, int)):
            raise ParametersDefinitionError(
'Attention, group_size_min must be an integer!')
        if group_size_min < 1:
            raise ParametersDefinitionError(
'Attention, group_size_min must be strictly positive!')
        if group_size_min < constraints.n_contig:
            raise ParametersDefinitionError(
'Attention, group_size_min must be at least equal to n_contiguity!')
        # maximum number of plies per group for each outer loop of LAYLA
        if isinstance(group_size_max, int):
            group_size_max = group_size_max*np.ones((n_outer_step,))
        self.group_size_max = group_size_max.astype(int)
        for el in group_size_max:
            if group_size_min > el:
                raise ParametersDefinitionError(
'Attention, group_size_min must smaller than group_size_max!')
        if group_size_max.size < n_outer_step:
            raise ParametersDefinitionError('''
Attention, the vector group_size_max should have as many elements as the
number of outer loops in LAYLA!''')

        # Branching limit for global pruning during ply orientation
        # optimisation
        self.global_branching_limit = np.around(global_branching_limit)
        if not isinstance(global_branching_limit, int):
            raise ParametersDefinitionError("""
Attention, global_branching_limit must be an integer!""")
        if global_branching_limit < 1:
            raise ParametersDefinitionError("""
Attention, global_branching_limit must be strictly positive!""")

        # Branching limit for global pruning at the penultimate level during
        # ply orientation optimisation
        self.global_branching_limit_p = np.around(global_branching_limit_p)
        if not isinstance(global_branching_limit_p, int):
            raise ParametersDefinitionError("""
Attention, global_branching_limit_p must be an integer!""")
        if global_branching_limit_p < 1:
            raise ParametersDefinitionError("""
Attention, global_branching_limit_p must be strictly positive!""")

        # Branching limit for local pruning during ply orientation
        # optimisation
        self.local_branching_limit = np.around(local_branching_limit)
        if not isinstance(local_branching_limit, int):
            raise ParametersDefinitionError("""
Attention, local_branching_limit must be an integer!""")
        if local_branching_limit < 1:
            raise ParametersDefinitionError("""
Attention, local_branching_limit must be strictly positive!""")

        # Lamination parameters sensitivities from the first-lebel optimiser
        if not(isinstance(first_level_sensitivities, np.ndarray)) \
        and first_level_sensitivities.size == 12 \
        and first_level_sensitivities.dtype == float:
            raise ParametersDefinitionError("""
Attention, first_level_sensitivities must be a vector with 12 float components!
""")
        if [False for elem in first_level_sensitivities if elem < 0]:
                    raise ParametersDefinitionError("""
Attention, the elements of first_level_sensitivities must be positive!
""")
        self.first_level_sensitivities = first_level_sensitivities

        # Lamination parameters to be considered in the multi-objective
        # functions
        if not(isinstance(lampam_to_be_optimised, np.ndarray)) \
        and lampam_to_be_optimised.size == 12 \
        and lampam_to_be_optimised.dtype == float:
            raise ParametersDefinitionError("""
Attention, lampam_to_be_optimised must be a vector with 12 float components!
""")
        if [False for elem in lampam_to_be_optimised if elem not in (0,1)]:
            raise ParametersDefinitionError("""
Attention, the elements of lampam_to_be_optimised must be either 0 or 1!
""")
        self.lampam_to_be_optimised = lampam_to_be_optimised

        # calculation of the multi-objective function lamination parameter
        # weightings
        (self.lampam_weightings_final, self.lampam_weightings_ini) \
        = self.weights_calculation(
                self.first_level_sensitivities,
                self.lampam_to_be_optimised, constraints,
                self.coeff_bal_ipo, self.coeff_oopo)

        if np.isclose(self.lampam_weightings_final[0:4],
                      np.array([0, 0, 0, 0], float)).all():
            self.weighting_finalA = self.lampam_weightings_final[0:4]
        else:
            self.weighting_finalA = self.lampam_weightings_final[0:4] / sum(
                self.lampam_weightings_final[0:4])

        if np.isclose(self.lampam_weightings_final[8:12],
                      np.array([0, 0, 0, 0], float)).all():
            self.weighting_finalD = self.lampam_weightings_final[8:12]
        else:
            self.weighting_finalD = self.lampam_weightings_final[8:12] / sum(
                self.lampam_weightings_final[8:12])


    def weights_calculation(
            self, first_level_sensitivities, lampam_to_be_optimised,
            constraints, coeff_bal_ipo, coeff_oopo):
        '''
        Calculation of the objective function lamination parameter weightings
        '''
        sensitivities = first_level_sensitivities * lampam_to_be_optimised
        if sum(sensitivities) == 0:
            raise ParametersDefinitionError(
'Attention, objective function lamination parameter weightings are all 0s!')
        sensitivities = sensitivities / sum(sensitivities)
        lampam_weightings_ini = np.copy(sensitivities)
#        print('sensitivities', sensitivities)

        # Filtering zero lamination parameter weightings
        opti_sensitivities = np.ones(12)
        if constraints.sym:
            opti_sensitivities[4:8] = 0
        if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
            opti_sensitivities[3]= 0
            opti_sensitivities[7]= 0
            opti_sensitivities[11]= 0
        # wlevel1 * wlevel2
        sensitivities = sensitivities * opti_sensitivities
#        print('sensitivities', sensitivities)

        # in-plane and out-of-plane orthotropy requirements
        if constraints.ipo:
            sensitivities[2] = 0
            sensitivities[3] = 0
        if constraints.oopo :
            sensitivities[10] = 0
            sensitivities[11] = 0
        s = np.average(sensitivities[sensitivities != 0])

        if constraints.ipo and self.penalty_ipo_switch and constraints.oopo:
            if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
                sensitivities[2] = s*coeff_bal_ipo
                sensitivities[10] = s*coeff_oopo
            else:
                sensitivities[2] = s*coeff_bal_ipo
                sensitivities[3] = s*coeff_bal_ipo
                sensitivities[10] = s*coeff_oopo
                sensitivities[11] = s*coeff_oopo
        elif constraints.ipo and self.penalty_ipo_switch:
            if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
                sensitivities[2] = s*coeff_bal_ipo
                sensitivities[10] = s*coeff_oopo
            else:
                sensitivities[2] = s*coeff_bal_ipo
                sensitivities[3] = s*coeff_bal_ipo
        elif constraints.oopo:
            if set(constraints.set_of_angles) == set([0, 45, -45, 90]):
                sensitivities[10] = s*coeff_oopo
            else:
                sensitivities[10] = s*coeff_oopo
                sensitivities[11] = s*coeff_oopo
#        print('sensitivities', sensitivities)

        return sensitivities/np.sum(sensitivities), lampam_weightings_ini


    def __repr__(self):
        " Display object "

        return (f"""
Optimiser parameters:

    Minimum size for the groups of plies: {self.group_size_min}
    Maximum size for the groups of plies: {self.group_size_max}
    Number of outer steps: {self.n_outer_step}
    Repair for in-plane lamination parameter convergence: {self.repair_membrane_switch}
    Repair for out-of-plane lamination parameter convergence: {self.repair_flexural_switch}
    Penalty for the 10% rule (on ply counts): {self.penalty_10_pc_switch}
    Penalty for the 10% rule (on lampams): {self.penalty_10_lampam_switch}
    Penalty for balance: {self.penalty_bal_switch}
    Penalty for in-plane orthotropy: {self.penalty_ipo_switch}
    Coefficient for the penalty for the 10\% rule: {self.coeff_10}
    Coefficient for the penalty for in-plane orthotropy or balance: {self.coeff_bal_ipo}
    Coefficient for the penalty for out-of-plane orthotropy: {self.coeff_oopo}
    Percent_thickness_repair_membrane: {self.p_A}%
    N_plies_last_permut_diso_contig: {self.n_D1}
    N_shifts_tested_flexural_repair: {self.n_D2}
    N_repeat_flexural_repair: {self.n_D3}
    In-plane orthotropy threshold: {self.threshold_ipo}
    Out-of-plane orthotropy threshold: {self.threshold_oopo}
    Type of objective function: norm {self.type_obj_func}
    Branching limits:
        - for global pruning during ply orientation optimisation: {self.global_branching_limit}
        - for global pruning at the last level of ply orientation optimisation: {self.global_branching_limit_p}
        - for local pruning during ply orientation optimisation: {self.local_branching_limit}
    In-plane orthotropy threshold: {self.threshold_ipo}
    Lampam_weightings_final:
        A:   {self.lampam_weightings_final[0]:.2f}   {self.lampam_weightings_final[1]:.2f}   {self.lampam_weightings_final[2]:.2f}   {self.lampam_weightings_final[3]:.2f}
        B:   {self.lampam_weightings_final[4]:.2f}   {self.lampam_weightings_final[5]:.2f}   {self.lampam_weightings_final[6]:.2f}   {self.lampam_weightings_final[7]:.2f}
        D:   {self.lampam_weightings_final[8]:.2f}   {self.lampam_weightings_final[9]:.2f}   {self.lampam_weightings_final[10]:.2f}   {self.lampam_weightings_final[11]:.2f}
        """)

class ParametersDefinitionError(Exception):
    pass

if __name__ == "__main__":
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
        percent_45=10,
        percent_90=10,
        percent_135=10,
        set_of_angles=[0, 45, -45, 90])
    parameters = Parameters(
        constraints=constraints,
        first_level_sensitivities = np.array([
            8, 4, 2, 1, 0, 12, 0, 13, 0, 0, 0, 0]),
        lampam_to_be_optimised=np.array([
            1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]))
    print(parameters)