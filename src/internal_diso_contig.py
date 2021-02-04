# -*- coding: utf-8 -*-
"""
Function performing internal sorting for disorientation and contiguity
'within a sublaminate'

Created on Mon Jan 29 12:00:18 2018

@author: Noemie Fedon
"""
import numpy as np

def internal_diso_contig(angle, constraints, angle2 = None):
    '''
    returns only the stacking sequences that satisfy the contiguity and
    disorientation rules

    OUTPUTS

    - angle: the selected sublaminate stacking sequences line by
    line
    - angle2: the selected sublaminate stacking sequences line by
    line if a second sublaminate is given as input for angle2

    INPUTS

    - angle: the first sublaminate stacking sequences
    - angle:2 matrix storing the second sublaminate stacking sequences
    '''
    if angle.ndim == 1:
        angle = angle.reshape((1, angle.size))
    n_plies_group = angle.shape[1]

    delta = constraints.delta_angle

    # TO ENSURE DISORIENTATION
    if constraints.diso:
        # To ensure the disorientation constraint within groups of plies
        if angle2 is None:
            a = angle.shape[0]
            for ii in range(a)[::-1]:
                jj = 0
                while jj<=n_plies_group-2:
                    if abs(angle[ii, jj + 1]-angle[ii, jj])>delta \
                    and abs(180+angle[ii, jj + 1]-angle[ii, jj])>delta  \
                    and abs(angle[ii, jj + 1]-angle[ii, jj]-180)>delta:
#                        print(angle[ii, jj + 1], angle[ii, jj])
                        angle = np.delete(angle, np.s_[ii], axis=0)
                        break
                    else:
                        jj=jj+1
        else:
            a = angle.shape[0]
            for ii in range(a)[::-1]:
                jj = 0
                while jj<=n_plies_group-2:
                    if abs(angle[ii, jj + 1]-angle[ii, jj])>delta \
                    and abs(180+angle[ii, jj + 1]-angle[ii, jj])>delta \
                    and abs(angle[ii, jj + 1]-angle[ii, jj]-180)>delta:
                        angle = np.delete(angle, np.s_[ii], axis=0)
                        angle2 = np.delete(angle2, np.s_[ii], axis=0)
                        break
                    else:
                        jj=jj + 1
    # TO ENSURE CONTIGUITY
    if constraints.contig:
        # To ensure the contiguity constraint within groups of plies
        if not angle2 is None:
            if constraints.n_contig < n_plies_group:
                diff = n_plies_group-constraints.n_contig
                if constraints.n_contig == 2:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 3:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 1] \
                            and angle[ii, jj + 3]==angle[ii, jj + 1]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 4:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 1] \
                            and angle[ii, jj + 2]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 5:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 6:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5] \
                            and angle[ii, jj]==angle[ii, jj + 6]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                angle2 = np.delete(angle2, np.s_[ii], axis=0)
                                break
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
        else:
            if constraints.n_contig < n_plies_group:
                diff = n_plies_group-constraints.n_contig
                if constraints.n_contig == 2:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 3:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 4:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 5:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break
                elif constraints.n_contig == 6:
                    a = angle.shape[0]
                    for ii in range(a)[::-1]:
                        for jj in np.arange(diff):
                            if angle[ii, jj]==angle[ii, jj + 1] \
                            and angle[ii, jj]==angle[ii, jj + 2] \
                            and angle[ii, jj]==angle[ii, jj + 3] \
                            and angle[ii, jj]==angle[ii, jj + 4] \
                            and angle[ii, jj]==angle[ii, jj + 5] \
                            and angle[ii, jj]==angle[ii, jj + 6]:
                                angle = np.delete(angle, np.s_[ii], axis=0)
                                break
                else:
                    raise Exception(
                            'constraints.n_contig must be 2, 3, 4 or 5')
    return angle, angle2


def internal_diso_contig_mod(ss, constraints, table_int, middle_ply = 0):
    '''
returns True for stacking sequences satisfying the extended disorientation
and contiguity rule, False otherwise
    '''

    if constraints.laminate_scheme == 'SB':
        if table_int.shape[0] == 1:
            ss, _ = internal_diso_contig(ss, constraints)
            return bool(ss.size)
        else:
            ss1, _ = internal_diso_contig(
                ss[0: table_int[-1, -1] -1], constraints)
            constraints.n_contig += 1
            ss2, _ = internal_diso_contig(ss, constraints)
            constraints.n_contig -= 1
            return bool(ss1.size) and bool(ss2.size)
    elif constraints.laminate_scheme == 'U':
        ss, _ = internal_diso_contig(ss, constraints)
        return bool(ss.size)
    elif constraints.laminate_scheme == 'S':
        ss, _ = internal_diso_contig(ss, constraints)
        return bool(ss.size)
    elif constraints.laminate_scheme == 'B':
        if table_int.shape[0] <= 2:
            ss, _ = internal_diso_contig(ss, constraints)
            return bool(ss.size)
        elif table_int.shape[0] == 3:
            ss1, _ = internal_diso_contig(
                ss[0: int(table_int[-1,-2])-1], constraints)
            ss2, _ = internal_diso_contig(
                ss[table_int[-1, -2] + table_int[-1, -1] -1:], constraints)
            constraints.n_contig += 1
            ss3, _ = internal_diso_contig(ss, constraints)
            constraints.n_contig -= 1
            return bool(ss1.size) and bool(ss2.size) and bool(ss3.size)
        else:
            a = np.argmin(table_int[-2:,2])
            if a == 0:
                beg = table_int[-2,2]
            else:
                beg = table_int[-1,2]
            ss1, _ = internal_diso_contig(ss[0: int(beg) -1], constraints)
            ss2, _ = internal_diso_contig(
                ss[int(beg) + table_int[-1, -1] + table_int[-2, -1] -1:],
                constraints)
            constraints.n_contig += 1
            ss3, _ = internal_diso_contig(ss, constraints)
            constraints.n_contig -= 1
#                print('ss1', 0, int(beg))
#                print_ss(np.reshape(ss1, (ss1.size,)))
#
#                print('ss2', int(beg) + table_int[-1, -1] \
#                           + table_int[-2, -1] -1)
#                print_ss(np.reshape(ss2, (ss2.size,)))
#
#                print('ss3')
#                print_ss(np.reshape(ss3, (ss3.size,)))
            return bool(ss1.size) and bool(ss2.size) and bool(ss3.size)



if __name__ == "__main__":
    'Test'

    import sys
    sys.path.append(r'C:\BELLA_and_LAYLA')
    from LAYLA_V01.constraints import Constraints
    from divers.pretty_print import print_ss, print_list_ss

    constraints = Constraints()
    constraints.diso = True
    delta = 45
    constraints.contig = True
    constraints.n_contig = 5

    print('*** Test for the function internal_diso_contig ***\n')
    print('Input stacking sequences:\n')
    ss = np.array([[0, 0, 45, 0, 0, 45, 90, 90, 45, 90],
                   [0, 0, 0, 0, 0, 45, 90, 90, 45, 90]])
    print_list_ss(ss)
    test = internal_diso_contig(ss, constraints)[0]
    if test.shape[0]:
        print('Stacking sequences satisfying the rule:\n')
        print_list_ss(test)
    else:
        print('No stacking sequence satisfy the rule\n')


    print('\n*** Test for the function internal_diso_contig_mod ***\n')

    constraints = Constraints()
    constraints.diso = True
    constraints.delta_angle = 45
    constraints.contig = True
    constraints.n_contig = 5
    constraints.sym = True
    constraints.bal = True

    table_int = np.array([  [8  ,1, 25 , 8],
                                 [8, 73 ,57 , 8],
                                 [8 , 9 ,33  ,8],
                                 [8 ,65, 49 , 8],
                                 [8 ,17 ,41 , 8]])

    print('Input stacking sequence:\n')
    ss = np.array([  0, 0, 0, 0, 0, 45, 90, 90, 90, 90,
                  45, 0, 0, 0, 0, 0, -45, 90, 90, 45,
                  45, 90, -45, -45, 0, 0, 0, 0, 0, -45,
                  90, 90, 90, 90, -45, 0, 0, 0, 0, 0,
                  45, 45, 90, -45, -45, 90, 90, 45, 45, 0,
                 -45, 90, 90, 90, 45, 0, 90, 90, 90, 45,
            0, 0, 0, 0, 0, -45, 90, 90, 90, 45,
            0, -45, 0, 0, 0, 0, -45, 90, 90, 90,])
    print_ss(ss, 8)
    test = internal_diso_contig_mod(ss, constraints, table_int)

    print('Modified internal disorientation and contiguity rule verified?')
    print(test)




