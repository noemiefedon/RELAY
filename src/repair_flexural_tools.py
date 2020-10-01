# -*- coding: utf-8 -*-
"""
repair for flexural properties

-  Block
     class for blocks of plies oriented at the same fibre orientation

- calc_delta_lampamD_swap & calc_delta_lampamD_swap_1
    returns the out-of-plane lamination parameters variation due to the swap
    of ply groups, taking into account the two symmetric parts for symmetric
    laminates - only one panel accounted for

- calc_delta_lampamD_swap_2
    returns the out-of-plane lamination parameters variation due to the swap
    of ply groups, taking into account the two symmetric parts for symmetric
    laminates - account for several panels

- find_list_blocks
    divides a stacking sequence into blocks of plies at the same
    fibre orientation
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import operator
import numpy as np
import numpy.matlib

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_lampam, print_ss, print_list_ss
from src.disorientation import is_diso
from src.lampam_functions import calc_lampam

def calc_delta_lampamD_swap(
        angle_first, angle_second, pos_first, pos_second, n_plies,
        constraints):
    '''
    returns the out-of-plane lamination parameters variation due to the
    modifications of some ply fibre orientations, taking into account the two
    symmetric parts for symmetric laminates - only one panel accounted for

    OUTPUTS

    - delta_lampam_D: out-of-plane partial lamination parameters

    INPUTS

    - angle_first: fibre orientation of the first group of plies
    - angle_second: fibre orientation of the second group of plies
    - pos_first: position of the plies in the first group
    - pos_second: position of the plies in the second group
    - n_plies: ply count of the laminate
    - constraints: set of constraints
    '''
    cos_sin_first = constraints.cos_sin[
        constraints.ind_angles_dict[angle_first]].reshape((4, 1))
    cos_sin_second = constraints.cos_sin[
        constraints.ind_angles_dict[angle_second]].reshape((4, 1))

    n_plies_first = pos_first.size
    n_plies_second = pos_second.size

    # vector of moments of area of the order 0 for each ply
    z_0_first = np.ones(n_plies_first)
    z_0_second = np.ones(n_plies_second)
    # vector of second moments of area for each ply
    z_2_first = np.array((
        (- n_plies / 2) * z_0_first + pos_first)**3 \
        - ((-n_plies / 2) * z_0_first + (pos_first - 1))**3)
    z_2_second = np.array((
        (- n_plies / 2) * z_0_second + pos_second)**3 \
        - ((-n_plies / 2) * z_0_second + (pos_second - 1))**3)

    if constraints.sym:
        ## REMOVE CONTRIBUTION OF FIRST BLOCK
        delta_lampam_D = -np.array([
            (8/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_first, 1, n_plies_first),
                z_2_first
            )]).reshape((4,))

        ## ADD CONTRIBUTION OF FIRST BLOCK MOVED TO POSITION OF SECOND BLOCK
        delta_lampam_D += np.array([
            (8/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_first, 1, n_plies_second),
                z_2_second
            )]).reshape((4,))

        ## REMOVE CONTRIBUTION OF SECOND BLOCK
        delta_lampam_D -= np.array([
            (8/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_second, 1, n_plies_second),
                z_2_second
            )]).reshape((4,))

        ## ADD CONTRIBUTION OF SECOND BLOCK MOVED TO POSITION OF FIRST BLOCK
        delta_lampam_D += np.array([
            (8/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_second, 1, n_plies_first),
                z_2_first
            )]).reshape((4,))
    else:
        ## REMOVE CONTRIBUTION OF FIRST BLOCK
        delta_lampam_D = -np.array([
            (4/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_first, 1, n_plies_first),
                z_2_first
            )]).reshape((4,))

        ## ADD CONTRIBUTION OF FIRST BLOCK MOVED TO POSITION OF SECOND BLOCK
        delta_lampam_D += np.array([
            (4/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_first, 1, n_plies_second),
                z_2_second
            )]).reshape((4,))

        ## REMOVE CONTRIBUTION OF SECOND BLOCK
        delta_lampam_D -= np.array([
            (4/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_second, 1, n_plies_second),
                z_2_second
            )]).reshape((4,))

        ## ADD CONTRIBUTION OF SECOND BLOCK MOVED TO POSITION OF FIRST BLOCK
        delta_lampam_D += np.array([
            (4/n_plies**3)*np.matmul(
                np.matlib.repmat(cos_sin_second, 1, n_plies_first),
                z_2_first
            )]).reshape((4,))

#    # Filter for numerical approximations
#    sett = set([0, 45, -45, 90, -90, 135, -135])
#    if np.all([el in sett  for el in constraints.set_of_angles]):
#        delta_lampam_D[3] = 0

    return delta_lampam_D


def calc_delta_lampamD_swap_1(
        angle_first, angle_second, n_first_ply, n_second_ply, n_plies_group,
        n_plies, constraints):
    '''
    returns the out-of-plane lamination parameters variation due to the swap
    of ply groups, taking into account the two symmetric parts for symmetric
    laminates - only one panel accounted for

    INPUTS

    - angle_first: fibre orientation of the first block of plies
    - angle_second: fibre orientation of the second block of plies
    - n_first_ply is the position of the first ply of the first group
    - n_second_ply is the position of the first ply of the second group
    - n_plies_group: ply count of the block of plies
    - n_plies: ply count of the laminate
    - constraints: set of constraints
    '''
    return calc_delta_lampamD_swap(
        angle_first,
        angle_second,
        pos_first=n_first_ply + np.arange(n_plies_group),
        pos_second=n_second_ply + np.arange(n_plies_group),
        n_plies=n_plies,
        constraints=constraints)


class Block():
    " An object for a block of plies oriented at the same fibre orientation"
    def __init__(self, ID, angle, n_block_of_plies, first_ply_pos,
                 n_plies, angle_before, angle_after, constraints):
        self.ID = ID
        self.angle = angle
        self.n_block_of_plies = n_block_of_plies
        self.first_ply_pos = first_ply_pos

        if first_ply_pos < n_plies / 2:
            position_1 = first_ply_pos
        else:
            position_1 = n_plies - first_ply_pos - 1
        distance_1 = 2 * abs(position_1 - (n_plies / 2)) / n_plies
        if first_ply_pos + n_block_of_plies - 1 < n_plies / 2:
            position_2 = first_ply_pos + n_block_of_plies - 1
        else:
            position_2 = n_plies - first_ply_pos - n_block_of_plies
        distance_2 = 2 * abs(position_2 - (n_plies / 2)) / n_plies
        self.distance_middle = max(distance_1, distance_2)

        self.neighbour_angles = []
        if angle_before is not None:
            self.neighbour_angles.append(angle_before)
        if angle_after is not None and angle_before != angle_after:
            self.neighbour_angles.append(angle_after)
        self.neighbour_angles = self.neighbour_angles

        self.calc_possible_angles(constraints)

    def calc_possible_angles(self, constraints):
        """
        finds the ply angles that the ply block can be changed to, whilst still
        satisfying disorientation and contiguity
        """
        possible_angles = []
        for ang in constraints.set_of_angles:
            if len(self.neighbour_angles) == 1 \
            and ang != self.angle \
            and ang not in self.neighbour_angles \
            and (not constraints.diso \
                 or is_diso(ang, self.neighbour_angles[0],
                            constraints.delta_angle)):
                possible_angles.append(ang)
            if len(self.neighbour_angles) == 2 \
            and ang != self.angle \
            and ang not in self.neighbour_angles \
            and (not constraints.diso \
                 or is_diso(ang, self.neighbour_angles[0],
                            constraints.delta_angle)) \
            and (not constraints.diso \
                 or is_diso(ang, self.neighbour_angles[1],
                            constraints.delta_angle)):
                possible_angles.append(ang)
        self.possible_angles = possible_angles

    def update_possible_angles(self, constraints, list_blocks, midply=None):
        " update the block IDs"
        self.neighbour_angles = []
        for block in list_blocks:
            if block.ID == self.ID + 1 or block.ID == self.ID - 1:
                self.neighbour_angles.append(block.angle)
        if self.ID == 0 and midply is not None:
            self.neighbour_angles.append(midply)
        self.calc_possible_angles(constraints)

        for block1 in list_blocks:
            # update the block / ID + 1
            if block1.ID == self.ID + 1:
                block1.neighbour_angles = []
                for block2 in list_blocks:
                    if block2.ID == block1.ID + 1 \
                    or block2.ID == block1.ID - 1:
                        block1.neighbour_angles.append(block2.angle)
                block1.calc_possible_angles(constraints)
            # update the block / ID - 1
            if block1.ID == self.ID - 1:
                block1.neighbour_angles = []
                for block2 in list_blocks:
                    if block2.ID == block1.ID + 1 \
                    or block2.ID == block1.ID - 1:
                        block1.neighbour_angles.append(block2.angle)
                block1.calc_possible_angles(constraints)


    def __repr__(self):
        " Display object "
        return f"""
Block of {self.n_block_of_plies} plies oriented at {self.angle} deg
    ID: {self.ID}
    First ply position: {self.first_ply_pos}
    Neighbour ply orientations: {self.neighbour_angles}
    Possible angles for a swap: {self.possible_angles}
    Normalised distance from the middle surface: {self.distance_middle}
    """

def find_list_blocks(ss_ref, n_plies, constraints):
    """
    divides a stacking sequence into blocks of plies at the same
    fibre orientation
    """
    if constraints.sym:
        if n_plies % 2:
            ind_start = n_plies // 2 - 1
            while ss_ref[ind_start] == ss_ref[n_plies // 2]:
                ind_start -= 1
        else:
            ind_start = n_plies // 2 - 1
    else:
        ind_start = n_plies - 1
    list_blocks = []
    ID = 0
    n_block_of_plies = 1
    while ind_start != 0:
        if ss_ref[ind_start] == ss_ref[ind_start - 1]:
            ind_start -= 1
            n_block_of_plies += 1
        else:
            if ID == 0 \
            and (not constraints.sym or (constraints.sym and n_plies % 2 == 0)):
                new_block = Block(
                    ID=ID,
                    angle=ss_ref[ind_start],
                    n_block_of_plies=n_block_of_plies,
                    first_ply_pos=ind_start,
                    n_plies=n_plies,
                    angle_before=ss_ref[ind_start - 1],
                    angle_after=None, # first group has only one neighbour
                    constraints=constraints)
            else:
                new_block = Block(
                    ID=ID,
                    angle=ss_ref[ind_start],
                    n_block_of_plies=n_block_of_plies,
                    first_ply_pos=ind_start,
                    n_plies=n_plies,
                    angle_before=ss_ref[ind_start - 1],
                    angle_after=ss_ref[ind_start + n_block_of_plies],
                    constraints=constraints)
            list_blocks.append(new_block)
            n_block_of_plies = 1
            ind_start -= 1
            ID += 1

    if ID == 0 and not (constraints.sym and n_plies % 2):
        new_block = Block(
            ID=ID,
            angle=ss_ref[ind_start],
            n_block_of_plies=n_block_of_plies,
            first_ply_pos=ind_start,
            n_plies=n_plies,
            angle_before=ss_ref[ind_start - 1],
            angle_after=None, # first group has only one neighbour
            constraints=constraints)
    else:
        new_block = Block(
            ID=ID,
            angle=ss_ref[ind_start],
            n_block_of_plies=n_block_of_plies,
            first_ply_pos=ind_start,
            n_plies=n_plies,
            angle_before=None, # last group has only one neighbour
            angle_after=ss_ref[ind_start + n_block_of_plies],
            constraints=constraints)
    list_blocks.append(new_block)

    list_blocks.sort(key=operator.attrgetter('distance_middle'), reverse=True)
    return list_blocks


if __name__ == "__main__":
    print('\n*** Test for the function find_list_blocks ***')
    constraints = Constraints(
        sym=False, contig=True, diso=True,
        n_contig=4, delta_angle=45,
        set_of_angles=[0, 45, -45, 90, 30, -30, 60, -60, 15, -15, 75, -75])
    print('Initial stacking sequence')
    ss = np.array([-75, -45, -60, -30, 0], int)
    if constraints.sym:
        ss = np.hstack((ss, 90, np.flip(ss)))
    print_ss(ss, 40)
    n_plies = ss.size
    list_blocks = find_list_blocks(ss, n_plies, constraints)
    for elem in list_blocks:
        print(elem)

    print('\n*** Test for the function calc_delta_lampamD_swap ***\n')
    constraints = Constraints(sym=False)
    print('Result:', calc_delta_lampamD_swap(
        angle_first=45,
        angle_second=-45,
        pos_first=np.array([26], int),
        pos_second=np.array([22], int),
        n_plies=30,
        constraints=constraints))

    print('\n*** Test for the function calc_delta_lampamD_swap_1 ***\n')
    print('Inputs:')
    ss_ini = np.array([0, 0, 0, 0, 0, 90, 45, 45, 45, 90, 90, 90], int)
    print('Initial stacking sequence')
    print_ss(ss_ini, 40)
    ss = np.array([45, 45, 45, 0, 0, 90, 0, 0, 0, 90, 90, 90], int)
    print('Final stacking sequence')
    print_ss(ss, 40)
    print('Expected result:',
          (calc_lampam(ss, constraints) - calc_lampam(ss_ini))[8:12])
    angle_first = 0
    angle_second = 45
    n_first_ply = 1
    n_second_ply = 7
    n_plies = 12
    n_plies_group = 3
    constraints = Constraints(sym=False)
    print('Result:', calc_delta_lampamD_swap_1(
        angle_first, angle_second, n_first_ply, n_second_ply,
        n_plies_group, n_plies, constraints))

