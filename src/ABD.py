# -*- coding: utf-8 -*-
"""
Functions for stiffness calculations based on the Classical Lamination Theory

- A_from_lampam
    calculates the in-plane stiffnesses from lamination parameters

- B_from_lampam
    calculates the coupling stiffnesses from lamination parameters

- D_from_lampam
    calculates the out-of-plane stiffnesses from lamination parameters

- a_from_lampam
    calculates the in-plane compliances from lamination parameters

- ad_from_lampam
    calculates the in-plane andout-of-plane compliances from lamination
    parameters

- d_from_lampam
    'calculates the out-of-plane compliances from lamination parameters

- filter_ABD
    filters the A, B, D stiffness matrices for 0 coefficients
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np

def filter_ABD(A=None, B=None, D=None, ply_t=None, sym=False):
    """
    filters the A, B, D stiffness matrices for 0 coefficients
    """
    B_is_filtered = False
    if sym and B is not None:
        B =  np.zeros((3, 3), float)
        B_is_filtered = True

    if A is not None:
        Amax = np.max(np.abs(A))
        A[abs(A) < 1e-10 * Amax] = 0

        if B is not None and not B_is_filtered:
            if ply_t is None:
                raise Exception('Ply thickness nedded for B filtering')
            Bmax = Amax * ply_t / 4
            B[abs(B) < 1e-10 * Bmax] = 0

    if D is not None:
        Dmax = np.max(np.abs(D))
        D[abs(D) < 1e-10 * Dmax] = 0

        if B is not None and not B_is_filtered:
            if ply_t is None:
                raise Exception('Ply thickness nedded for B filtering')
            Bmax = 3 * Dmax / ply_t
            B[abs(B) < 1e-10 * Bmax] = 0

    return 0




def A_from_lampam(lampam, mat):
    'calculates the in-plane stiffnesses from lamination parameters'
    A11 = mat.ply_t*(mat.U1 + mat.U2*lampam[0] + mat.U3*lampam[1])
    A12 = mat.ply_t*(- mat.U3*lampam[1] + mat.U4)
    A22 = mat.ply_t*(mat.U1 - mat.U2*lampam[0] + mat.U3*lampam[1])
    A66 = mat.ply_t*(- mat.U3*lampam[1] + mat.U5)
    A16 = mat.ply_t*(0.5*mat.U2*lampam[2] + mat.U3*lampam[3])
    A26 = mat.ply_t*(0.5*mat.U2*lampam[2] - mat.U3*lampam[3])
    return np.array([[A11, A12, A16],
                     [A12, A22, A26],
                     [A16, A26, A66]])


def B_from_lampam(lampam, mat, sym=False):
    'calculates the coupling stiffnesses from lamination parameters'
    if sym:
        return np.array([[0, 0, 0],
                         [0, 0, 0],
                         [0, 0, 0]])
    a = mat.ply_t**2/4
    B11 = a*(mat.U2*lampam[4] + mat.U3*lampam[5])
    B12 = a*(- mat.U3*lampam[5])
    B22 = a*(- mat.U2*lampam[4] + mat.U3*lampam[5])
    B66 = a*(- mat.U3*lampam[5])
    B16 = a*(0.5*mat.U2*lampam[6] + mat.U3*lampam[7])
    B26 = a*(0.5*mat.U2*lampam[6] - mat.U3*lampam[7])
    return np.array([[B11, B12, B16],
                     [B12, B22, B26],
                     [B16, B26, B66]])


def D_from_lampam(lampam, mat):
    'calculates the out-of-plane stiffnesses from lamination parameters'
    a = mat.ply_t**3/12
    D11 = a*(mat.U1 + mat.U2*lampam[8] + mat.U3*lampam[9])
    D12 = a*(- mat.U3*lampam[9] + mat.U4)
    D22 = a*(mat.U1 - mat.U2*lampam[8] + mat.U3*lampam[9])
    D66 = a*(- mat.U3*lampam[9] + mat.U5)
    D16 = a*(0.5*mat.U2*lampam[10] + mat.U3*lampam[11])
    D26 = a*(0.5*mat.U2*lampam[10] - mat.U3*lampam[11])
    return np.array([[D11, D12, D16],
                     [D12, D22, D26],
                     [D16, D26, D66]])


def a_from_lampam(lampam, mat, sym=False):
    'calculates the in-plane compliances from lamination parameters'
    if sym:
        A = A_from_lampam(lampam, mat)
        a = np.linalg.inv(A)
    else:
        A = A_from_lampam(lampam, mat)
        B = B_from_lampam(lampam, mat, sym)
        D = D_from_lampam(lampam, mat)
        a = np.linalg.inv(A - B@(np.linalg.inv(D))@B)
    return a


def ad_from_lampam(lampam, mat, sym=False):
    """
    calculates the in-plane andout-of-plane compliances from lamination
    parameters
    """
    if sym:
        A = A_from_lampam(lampam, mat)
        a = np.linalg.inv(A)
        D = D_from_lampam(lampam, mat)
        d = np.linalg.inv(D)
    else:
        A = A_from_lampam(lampam, mat)
        B = B_from_lampam(lampam, mat, sym)
        D = D_from_lampam(lampam, mat)
        a = np.linalg.inv(A - B@(np.linalg.inv(D))@B)
        d = np.linalg.inv(D - B@(np.linalg.inv(A))@B)
    return a, d


def d_from_lampam(lampam, mat, sym=False):
    'calculates the out-of-plane compliances from lamination parameters'
    if sym:
        D = D_from_lampam(lampam, mat)
        d = np.linalg.inv(D)
    else:
        A = A_from_lampam(lampam, mat)
        B = B_from_lampam(lampam, mat, sym)
        D = D_from_lampam(lampam, mat)
        d = np.linalg.inv(D - B@(np.linalg.inv(A))@B)
    return d


if __name__ == "__main__":
    print('*** Test for the functions of the module ABD ***\n')
    import sys
    sys.path.append(r'C:\BELLA_and_LAYLA')
    from CLA.lampam_functions import calc_lampam
    from BELLA.materials import Material
    ss = np.array([0, 45, 45, -45, -45, 90, 0, 90])
    ss = np.hstack((ss, np.flip(ss, axis=0)))
    lampam = calc_lampam(ss)
    E11 = 130e9
    E22 = 9e9
    nu12 = 0.3
    G12 = 4e9
    threshold = 0.01
    mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12)
    sym = False
    A = A_from_lampam(lampam, mat)
    B = B_from_lampam(lampam, mat, sym)
    D = D_from_lampam(lampam, mat)
    filter_ABD(A, B, D, sym=sym, ply_t=0.0002)
    a = a_from_lampam(lampam, mat, sym)
    d = d_from_lampam(lampam, mat, sym)
    a2, d2 = ad_from_lampam(lampam, mat, sym)
    print('A = ', A)
    print('B = ', B)
    print('D = ', D)
    print('a = ', a)
    print('a2 = ', a2)
    print('d = ', d)
    print('d2 = ', d2)
