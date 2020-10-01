# -*- coding: utf-8 -*-
"""
Functions to calculate lamination parameters

- filter_lampam
    returns lamination parameters where a few numerical approximations have
    been filtered regarding the constraints for laminate symmetry and balance,
    and the fibre orientations used

- calc_lampam, calc_lampam_2
    calculates the lamination parameters of one or more laminates from their
    stacking sequence

- calc_lampam_sym
    returns the lamination parameters of one or more symmetric laminates
    with even ply counts
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\RELAY')
from src.constraints import Constraints
from src.pretty_print import print_lampam, print_ss, print_list_ss
from src.parameters import Parameters

def filter_lampam(lampam, constraints):
    """
    returns lamination parameters where a few numerical approximations have
    been filtered regarding the constraints for laminate symmetry and balance,
    and the fibre orientations used

    OUTPUTS

    - lampam: array storing the filtered lamination parameters

    INPUTS

    - lampam: array storing the laminate lamination parameters
    - constraints: design and manufacturing guidelines
    """
    if lampam.ndim == 1:
        if constraints.sym:
            lampam[4:8] = 0
#        if constraints.bal:
#            lampam[2:4] = 0
        if constraints.n_set_of_angles:
            sett = set([0, 45, -45, 90, -90, 135, -135])
            if np.all([el in sett  for el in constraints.set_of_angles]):
                lampam[3] = 0
                lampam[7] = 0
                lampam[11] = 0
    elif lampam.ndim == 2:
        if constraints.sym:
            lampam[:, 4:8] = 0
#        if constraints.bal:
#            lampam[:, 2:4] = 0
        if constraints.n_set_of_angles:
            sett = set([0, 45, -45, 90, -90, 135, -135])
            if np.all([el in sett  for el in constraints.set_of_angles]):
                lampam[:, 3] = 0
                lampam[:, 7] = 0
                lampam[:, 11] = 0
    else:
        raise Exception('This should not happen.')
    return lampam




def calc_lampam_2(ss):
    """
    returns the lamination parameters of one or more laminates

    OUTPUTS

    - lampam: laminate lamination parameters

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if isinstance(ss, list):
        lampam = np.zeros((len(ss), 12), float)
        for index in range(len(ss)):
            lampam[index] = calc_lampam_2(ss[index])
        return lampam
    if ss.ndim == 2 and ss.shape[0] > 1:
        lampam = np.zeros((ss.shape[0], 12), float)
        for index in range(ss.shape[0]):
            lampam[index] = calc_lampam_2(ss[index])
        return lampam

    n_plies_in_panels = np.size(ss) # laminate ply count

    theta2 = np.deg2rad(2*ss.astype(float))
    theta4 = 2*theta2
    cos_sin = np.concatenate((
        np.cos(theta2),
        np.cos(theta4),
        np.sin(theta2),
        np.sin(theta4))).reshape((4, n_plies_in_panels))

    for_the_top = np.arange(n_plies_in_panels)
    z_0 = np.ones(n_plies_in_panels)
    z_2 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**3 \
        - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
    z_1 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**2 \
        - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**2

    return np.array([
        (1/n_plies_in_panels)*np.matmul(cos_sin, z_0),
        (2/n_plies_in_panels**2)*np.matmul(cos_sin, z_1),
        (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

def calc_lampam(ss, constraints=None):
    """
    returns the lamination parameters of one or more laminates

    OUTPUTS

    - lampam: laminate lamination parameters

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if constraints is None:
        return calc_lampam_2(ss)

    if isinstance(ss, list):
        lampam = np.zeros((len(ss), 12), float)
        for index in range(len(ss)):
            lampam[index] = calc_lampam(ss[index], constraints)
        return lampam
    if ss.ndim == 2 and ss.shape[0] > 1:
        lampam = np.zeros((ss.shape[0], 12), float)
        for index in range(ss.shape[0]):
            lampam[index] = calc_lampam(ss[index], constraints)
        return lampam
    n_plies_in_panels = np.size(ss) # laminate ply count

    if not constraints.sym:
        cos_sin = np.empty((4, n_plies_in_panels), float)
        for ind in range(n_plies_in_panels):
            cos_sin[:, ind] = np.copy(constraints.cos_sin[
                constraints.ind_angles_dict[ss[ind]]].reshape((4, )))

        for_the_top = np.arange(n_plies_in_panels)
        z_0 = np.ones(n_plies_in_panels)
        z_2 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**3 \
            - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
        z_1 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**2 \
            - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**2
        return np.array([
            (1/n_plies_in_panels)*np.matmul(cos_sin, z_0),
            (2/n_plies_in_panels**2)*np.matmul(cos_sin, z_1),
            (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

    cos_sin = np.empty((4, np.size(ss) // 2), float)
    for ind in range(np.size(ss) // 2):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4,))

    for_the_top = np.arange(np.size(ss) // 2)
    z_0 = np.ones(np.size(ss) // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3
    lampam = np.array([
        (2/n_plies_in_panels)*np.matmul(cos_sin, z_0),
        np.array([0, 0, 0, 0]),
        (8/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

    if np.size(ss) % 2:
        cos_sin_mid = constraints.cos_sin[
            constraints.ind_angles_dict[ss[n_plies_in_panels // 2]]]
        lampam += np.array([
            (1/n_plies_in_panels)*cos_sin_mid,
            np.zeros((4,), dtype=float),
            (1/n_plies_in_panels**3)*cos_sin_mid]).reshape(12)
    return lampam

def calc_lampam_sym(ss, constraints):
    """
    returns the lamination parameters of one or more symmetric laminates
    with even ply counts from thehalf stacking sequences

    OUTPUTS

    - lampam: laminate lamination parameters

    INPUTS

    - ss: half laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if isinstance(ss, list):
        lampam = np.zeros((len(ss), 12), float)
        for index in range(len(ss)):
            lampam[index] = calc_lampam_sym(ss[index], constraints)
        return lampam
    if ss.ndim == 2 and ss.shape[0] > 1:
        lampam = np.zeros((ss.shape[0], 12), float)
        for index in range(ss.shape[0]):
            lampam[index] = calc_lampam_sym(ss[index], constraints)
        return lampam

    n_plies_in_panels = 2 * np.size(ss) # laminate ply count

    cos_sin = np.empty((4, n_plies_in_panels // 2), float)
    for ind in range(n_plies_in_panels // 2):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4, ))

    for_the_top = np.arange(n_plies_in_panels // 2)
    z_0 = np.ones(n_plies_in_panels // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3
    lampam = np.array([
        (2 / n_plies_in_panels)*np.matmul(cos_sin, z_0),
        np.array([0, 0, 0, 0]),
        (8 / n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)
    return lampam



if __name__ == "__main__":
    constraints = Constraints(
        sym=False,
        set_of_angles=np.array([0, 45, 90, -45]))
    parameters = Parameters(
        constraints=constraints, group_size_min=10, group_size_max=20)

    print('\n*** Test for the function filter_lampam ***\n')
#    lampam = np.arange(1, 13)
#    lampam = np.arange(1, 25).reshape((2, 12))
#    print('Lamination parameters:\n')
#    lampam = filter_lampam(lampam, constraints)
#    print_lampam(lampam[1])

    print('\n*** Test for the function calc_lampam ***\n')
    print('Input stacking sequence:\n')
    ss = np.array([45, 90, 45, 45, 0, -45, -45, 0, 90, -45])
    print(f'{ss}\n')
    print('Lamination parameters:\n')
    lampam = calc_lampam(ss)
    print_lampam(lampam)
