# -*- coding: utf-8 -*-
"""
Class for design and manufacturing constraints
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

class Constraints():
    """
    An object for storing a set of design & manufacturing constraints
    """
    def __init__(self,
                 sym=False,
                 bal=False,
                 ipo=False,
                 oopo=False,
                 dam_tol=False,
                 rule_10_percent=False,
                 percent_0=0,
                 percent_45=0,
                 percent_90=0,
                 percent_135=0,
                 percent_45_135=0,
                 diso=False,
                 contig=False,
                 n_contig=1,
                 delta_angle=45,
                 dam_tol_rule=2,
                 set_of_angles=np.array([0, 45, 90, -45])):
        " Create object with inactive constraints by default"
        # Symmetry
        self.sym = sym
        if not isinstance(sym, bool):
            raise ConstraintDefinitionError("""
Attention, sym must be a boolean!""")

        # Damage tolerance
        self.dam_tol = dam_tol
        if not isinstance(dam_tol, bool):
            raise ConstraintDefinitionError("""
Attention, dam_tol must be a boolean!""")
        # which damage tolerance rule?
        if dam_tol:
            if dam_tol_rule not in [1, 2, 3]:
                raise ConstraintDefinitionError("""
Choose either:
    - damage tolerance rule 1: one outer ply at + or -45 deg at the laminate
    surfaces (2 plies in total)
    - damage tolerance rule 2: two [+45, -45] or [-45, +45] at the laminate
    surfaces (4 plies in total)
    - damage tolerance rule 3: [+45,-45] [-45,+45] [+45,+45] or [-45,-45] at
    the laminate surfaces (4 plies in total)""")
            self.dam_tol_rule = dam_tol_rule
        else:
            self.dam_tol_rule = 0

        # 10% rule: requirement on minimum percentage of plies
        # in the 0/+-45/90  deg directions
        self.rule_10_percent = rule_10_percent
        if not isinstance(rule_10_percent, bool):
            raise ConstraintDefinitionError("""
Attention, rule_10_percent must be a boolean!""")
        # set the limit percentages for the 10% rule
        self.set_percentages(
            percent_0, percent_45, percent_90, percent_135, percent_45_135)

        # Disorientation rule: the change of angles between two consecutive
        # plies should not exceed delta_angle
        self.diso = diso
        if not isinstance(diso, bool):
            raise ConstraintDefinitionError("""
Attention, diso must be a boolean!""")
        if diso:
            self.delta_angle = delta_angle
        else:
            self.delta_angle = 180

        # Contiguity rule: no more that 'n_contig' plies with same fibre
        # orientation should be next to each other
        self.contig = contig
        if not isinstance(contig, bool):
            raise ConstraintDefinitionError("""
Attention, contig must be a boolean!""")
        if not contig:
            self.n_contig = 1
            self.n_contig_c = 1e10
        else:
            self.n_contig = n_contig
            self.n_contig_c = n_contig
            if n_contig == 1:
                raise ConstraintDefinitionError("""
Not allowing two adjacent plies to have the same fibre orientation is too
restricting!
""")

        # Set the fibre orientations
        self.set_fibre_orientations(set_of_angles, rule_10_percent)

        # balance and orthrotropy requirements
        if not isinstance(ipo, bool):
            raise ConstraintDefinitionError(f"""
Attention, the in-plane orthotropy requirement {ipo} must be a boolean!""")
        if not isinstance(oopo, bool):
            raise ConstraintDefinitionError(f"""
Attention, out-of-plane orthotropy requirements {oopo} must be a boolean!""")
        if not isinstance(bal, bool):
            raise ConstraintDefinitionError(f"""
Attention, the balance requirement {bal} must be a boolean!""")
        if not isinstance(bal, bool):
            raise ConstraintDefinitionError(f"""
Attention, the balance requirement {bal} must be a boolean!""")
        self.ipo = ipo # In-plane orthotropy requirement
        self.oopo = oopo # Out-of-plane orthotropy requirement
        self.bal = bal # Balance requirements
        if self.bal:
            self.ipo = True

        if self.ipo:
            if self.sym:
                self.laminate_scheme = 'S'
                self.laminate_scheme_test = 'SB'
            else:
                self.laminate_scheme = 'U'
                self.laminate_scheme_test = 'B'
        else:
            if self.sym:
                self.laminate_scheme = 'S'
                self.laminate_scheme_test = 'S'
            else:
                self.laminate_scheme = 'U'
                self.laminate_scheme_test = 'U'

    def set_percentages(self, percent_0, percent_45, percent_90, percent_135,
                        percent_45_135):
        'sets the percentages for the 10% rule'
        # Minimum percentage of 0deg plies
        self.percent_0 = percent_0
        # Minimum percentage of 45deg plies
        self.percent_45 = percent_45
        # Minimum percentage of 90deg plies
        self.percent_90 = percent_90
        # Minimum percentage of -45deg plies
        self.percent_135 = percent_135
        # Minimum percentage of +-45deg plies
        self.percent_45_135 = max(percent_45_135, percent_45 + percent_135)
        if not isinstance(percent_0, (float, int)):
            raise ConstraintDefinitionError("""
Attention, percent_0 must be a number (float or integer)!""")
        if percent_0 < 0 or percent_0 > 100:
            raise ConstraintDefinitionError("""
Attention, percent_0 is a pecentage and must between 0 and 100!""")
        if not isinstance(percent_45, (float, int)):
            raise ConstraintDefinitionError("""
Attention, percent_45 must be a number (float or integer)!""")
        if percent_45 < 0 or percent_45 > 100:
            raise ConstraintDefinitionError("""
Attention, percent_45 is a pecentage and must between 0 and 100!""")
        if not isinstance(percent_90, (float, int)):
            raise ConstraintDefinitionError("""
Attention, percent_90 must be a number (float or integer)!""")
        if percent_90 < 0 or percent_90 > 100:
            raise ConstraintDefinitionError("""
Attention, percent_90 is a pecentage and must between 0 and 100!""")
        if not isinstance(percent_135, (float, int)):
            raise ConstraintDefinitionError("""
Attention, percent_135 must be a number (float or integer)!""")
        if percent_135 < 0 or percent_135 > 100:
            raise ConstraintDefinitionError("""
Attention, percent_135 is a pecentage and must between 0 and 100!""")
        if not isinstance(percent_45_135, (float, int)):
            raise ConstraintDefinitionError("""
Attention, percent_45_135 must be a number (float or integer)!""")
        if percent_45_135 < 0 or percent_45_135 > 100:
            raise ConstraintDefinitionError("""
Attention, percent_45_135 is a pecentage and must between 0 and 100!""")
        if percent_0 + percent_45 + percent_90 + percent_135 > 100 \
        or percent_0 + percent_90 + percent_45_135 > 100:
            print("""
Total percentage for the plies in the directions 0/+-45/90 greater than 100!
""")
        if self.rule_10_percent:
            self.percent_0 = self.percent_0/100
            self.percent_45 = self.percent_45/100
            self.percent_90 = self.percent_90/100
            self.percent_135 = self.percent_135/100
            self.percent_45_135 = self.percent_45_135/100
            self.percent_tot = max(
                self.percent_0 + self.percent_45 \
                + self.percent_90 + self.percent_135,
                self.percent_0 + self.percent_45_135 + self.percent_90)
        else:
            self.percent_0 = 0
            self.percent_45 = 0
            self.percent_90 = 0
            self.percent_135 = 0
            self.percent_45_135 = 0
            self.percent_tot = 0

    def set_fibre_orientations(self, set_of_angles, rule_10_percent):
        'sets the aloowed fibre orientations'
        # Allowed fibre orientations
        self.set_of_angles = np.unique(set_of_angles)
        if (self.set_of_angles > 90).any() \
        or (self.set_of_angles <= -90).any():
            raise Exception(r"""
The allowed fibre angles must be between -90 (excluded) and 90 degrees
(included).""")
        sett = set(self.set_of_angles)
        if (0 not in sett or 45 not in sett or 90 not in sett \
            or -45 not in sett) and rule_10_percent:
            raise Exception(r"""
The 10% rule is only applicable if the fibre orientations 0, +45, 90, -90
are allowed.""")
        # usefull data for the 10% rule application
        if 0 in sett:
            # index of the 0deg fibre direction in set_of_angles
            self.index0 = np.where(self.set_of_angles == 0)[0][0]
        if 45 in sett:
            # index of the 45deg fibre direction in set_of_angles
            self.index45 = np.where(self.set_of_angles == 45)[0][0]
        if 90 in sett:
            # index of the 90deg fibre direction in set_of_angles
            self.index90 = np.where(self.set_of_angles == 90)[0][0]
        if -45 in sett:
            # index of the -45deg fibre direction in set_of_angles
            self.index135 = np.where(self.set_of_angles == -45)[0][0]
        # usefull data for counting plies in each fibre direction
        angles_dict = dict() # Dictionary to retrieve angles
        ind_angles_dict = dict() # Dictionary to retrieve indices
        for index, angle in enumerate(self.set_of_angles):
            angles_dict[index] = angle
            ind_angles_dict[angle] = index
        if 90 in sett:
            ind_angles_dict[-90] = ind_angles_dict[90]
        self.angles_dict = angles_dict
        self.ind_angles_dict = ind_angles_dict

        # usefull data for positive angles only
        pos_angles = np.unique(np.abs(self.set_of_angles))
#        pos_angles_dict = dict() # Dictionary to retrieve angles
        indices_pos_angles_dict = dict() # Dictionary to retrieve indices
        for index, angle in enumerate(pos_angles):
#            pos_angles_dict[index] = angle
            indices_pos_angles_dict[angle] = index
        self.pos_angles = pos_angles
        self.indices_pos_angles_dict = indices_pos_angles_dict

        # check that angled plies all have their balanced counterpart
        for angle in self.set_of_angles:
            if angle != 90 and -angle not in self.set_of_angles:
                raise Exception(f"""
Missing input fibre orientation {-angle} to have both angle plies +-{angle}.
""")
        # Identification of the panels which are not balanced by counting the
        # difference of ply counts for the angled plies.
        angles_bal = self.set_of_angles[
            [index for index in range(self.set_of_angles.size) \
             if self.set_of_angles[index] > 0 \
             and self.set_of_angles[index] < 90]]
        # angles_bal: each row has three values:
        # - fibre orientation of angle ply theta
        # - index of +theta in constraints.set_of_angles
        # - index of -theta in constraints.set_of_angles
        self.angles_bal = np.array([[
            elem,
            self.ind_angles_dict[elem],
            self.ind_angles_dict[-elem]] for elem in angles_bal])

        # data used for repair for balance
        self.indices_bal = dict() # Dictionary to retrieve indices
        for index, angle in enumerate(self.angles_bal[:, 0]):
            self.indices_bal[angle] = index

        # number of allowed fibre orientations
        self.n_set_of_angles = len(set(set_of_angles))
        if self.n_set_of_angles != len((set_of_angles)):
            raise Exception("""
Repeated angles in the set of allowed fibre orientations set_of_angles""")

        # usefull data to implement the 10% rule
        self.angles_10 = [0, 90, 45, -45]
        # dictionary to retrieve indices related to the 10% rule
        self.indices_10 = dict()
        for index, angle in enumerate(self.angles_10):
            self.indices_10[angle] = index

        # usefull data to avoid repetitive caluclations of cosines and sines
        self.cos_sin = np.empty((self.n_set_of_angles, 4), float)
        for ind_angle, angle in enumerate(self.set_of_angles):
            self.cos_sin[ind_angle, :] = np.hstack((
                np.cos(np.deg2rad(2*float(angle))),
                np.cos(np.deg2rad(4*float(angle))),
                np.sin(np.deg2rad(2*float(angle))),
                np.sin(np.deg2rad(4*float(angle)))))

    def __repr__(self):
        " Display object "
        return f'''
Constraints:

    Symmetry: {self.sym}
    Balance: {self.bal}
    In-plane orthotropy requirement: {self.ipo}
    Out-of-plane orthotropy requirement: {self.oopo}
    Damage tolerance constraint: {self.dam_tol}
        Damage tolerance rule: {self.dam_tol_rule}
    10 percent rule: {self.rule_10_percent}
        percentage 0: {self.percent_0*100:.2f} %
        percentage 45: {self.percent_45*100:.2f} %
        percentage 90: {self.percent_90*100:.2f} %
        percentage -45: {self.percent_135*100:.2f} %
        percentage +-45: {self.percent_45_135*100:.2f} %
        percentage of plies concerned by the 10 % rule: {self.percent_tot*100:.2f} %
    Disorientation rule: {self.diso}
        delta_angle: {self.delta_angle}
    Contiguity rule: {self.contig}
        n_contig: {self.n_contig}
    Set of angles: {self.set_of_angles}
        Number of allowed fibre orientations: {self.n_set_of_angles}
'''

class ConstraintDefinitionError(Exception):
    " Errors during the constraints definition"

if __name__ == "__main__":
    constraints = Constraints()
    print(constraints)
