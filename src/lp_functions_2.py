# -*- coding: utf-8 -*-
"""
Functions calculating lamination parameters
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def calc_lampamA(ss, constraints):
    """
    returns the in-plane lamination parameters of one or more laminates

    Beware: stacking sequences should be complete for symmetric laminates !

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if isinstance(ss, list):
        lampamA = np.empty((len(ss), 4), float)
        for ind_panel in range(len(ss)):
            lampamA[ind_panel] = calc_lampamA(ss[ind_panel], constraints)
        return lampamA

    n_plies_in_panels = np.size(ss) # laminate ply count
    cos_sin = np.zeros((4,), float)
    if not constraints.sym:
        for ind in range(n_plies_in_panels):
            cos_sin += constraints.cos_sin[
                constraints.ind_angles_dict[ss[ind]]].reshape((4, ))
        return (1 / n_plies_in_panels) * cos_sin
    for ind in range(np.size(ss) // 2):
        cos_sin += constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4,))
    if np.size(ss) % 2:
        cos_sin += 0.5 * constraints.cos_sin[
            constraints.ind_angles_dict[
                ss[n_plies_in_panels // 2]]].reshape((4,))
    return (2 / n_plies_in_panels) * cos_sin


def calc_lampamD(ss, constraints):
    """
    returns the out-of-plane lamination parameters of a laminate

    Beware: the stacking sequence should be complete for symmetric laminates !

    OUTPUTS

    - lampam: laminate out-of-plane lamination parameters

    INPUTS

    - ss: laminate stacking sequence
    - constraints: design and manufacturing guidelines
    """
    n_plies_in_panels = np.size(ss) # laminate ply count

    if not constraints.sym:
        cos_sin = np.empty((4, n_plies_in_panels), float)
        for ind in range(n_plies_in_panels):
            cos_sin[:, ind] = constraints.cos_sin[
                constraints.ind_angles_dict[ss[ind]]].reshape((4, ))

        for_the_top = np.arange(n_plies_in_panels)
        z_0 = np.ones(n_plies_in_panels)
        z_2 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**3 \
            - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
        return np.array([
            (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape((4,))

    cos_sin = np.empty((4, np.size(ss) // 2), float)
    for ind in range(np.size(ss) // 2):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4,))

    for_the_top = np.arange(np.size(ss) // 2)
    z_0 = np.ones(np.size(ss) // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3
    lampam = np.array([
        (8/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape((4,))
    if np.size(ss) % 2:
        cos_sin_mid = constraints.cos_sin[
            constraints.ind_angles_dict[ss[n_plies_in_panels // 2]]]
        lampam += np.array([
            (1/n_plies_in_panels**3)*cos_sin_mid]).reshape((4,))
    return lampam
