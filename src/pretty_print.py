# -*- coding: utf-8 -*-
"""
Pretty prints

    - print_lampam: prints all the components of a lamination parameters array
    - print_ss: prints all the components of a stacking sequence array
    - print_list: prints all the components of a list of floats
    - print_list_ss: prints all the components of a list of stacking sequences
"""

import numpy as np

def print_lampam(lampam_1, lampam_2=None, diff=False):
    " Prints all the components of a lamination parameters array"
    if lampam_2 is None:
        lampam_1 = lampam_1.reshape(12)
        print(f'lampam 1 :{lampam_1[0]:.5f}')
        print(f'lampam 2 :{lampam_1[1]:.5f}')
        print(f'lampam 3 :{lampam_1[2]:.5f}')
        print(f'lampam 4 :{lampam_1[3]:.5f}')
        print('\033[37m', end='')
        print(f'lampam 5 :{lampam_1[4]:.5f}')
        print(f'lampam 6 :{lampam_1[5]:.5f}')
        print(f'lampam 7 :{lampam_1[6]:.5f}')
        print(f'lampam 8 :{lampam_1[7]:.5f}')
        print('\x1b[0m', end='')
        print(f'lampam 9 :{lampam_1[8]:.5f}')
        print(f'lampam 10:{lampam_1[9]:.5f}')
        print(f'lampam 11:{lampam_1[10]:.5f}')
        print(f'lampam 12:{lampam_1[11]:.5f}\n')
    else:
        lampam_1 = lampam_1.reshape(12)
        lampam_2 = lampam_2.reshape(12)
        if not diff:
            print('lampam 1  : ', end='')
            if lampam_1[0] > 0:
                print(' ', end='')
            print(f'{lampam_1[0] :.5f} : ', end='')
            if lampam_2[0] > 0:
                print(' ', end='')
            print(f'{lampam_2[0]:.5f}')

            print('lampam 2  : ', end='')
            if lampam_1[1] > 0:
                print(' ', end='')
            print(f'{lampam_1[1] :.5f} : ', end='')
            if lampam_2[1] > 0:
                print(' ', end='')
            print(f'{lampam_2[1]:.5f}')

            print('lampam 3  : ', end='')
            if lampam_1[2] > 0:
                print(' ', end='')
            print(f'{lampam_1[2] :.5f} : ', end='')
            if lampam_2[2] > 0:
                print(' ', end='')
            print(f'{lampam_2[2]:.5f}')

            print('lampam 4  : ', end='')
            if lampam_1[3] > 0:
                print(' ', end='')
            print(f'{lampam_1[3] :.5f} : ', end='')
            if lampam_2[3] > 0:
                print(' ', end='')
            print(f'{lampam_2[3]:.5f}')

            print('\033[37m', end='')

            print('lampam 5  : ', end='')
            if lampam_1[4] > 0:
                print(' ', end='')
            print(f'{lampam_1[4] :.5f} : ', end='')
            if lampam_2[4] > 0:
                print(' ', end='')
            print(f'{lampam_2[4]:.5f}')

            print('lampam 6  : ', end='')
            if lampam_1[5] > 0:
                print(' ', end='')
            print(f'{lampam_1[5] :.5f} : ', end='')
            if lampam_2[5] > 0:
                print(' ', end='')
            print(f'{lampam_2[5]:.5f}')

            print('lampam 7  : ', end='')
            if lampam_1[6] > 0:
                print(' ', end='')
            print(f'{lampam_1[6] :.5f} : ', end='')
            if lampam_2[6] > 0:
                print(' ', end='')
            print(f'{lampam_2[6]:.5f}')

            print('lampam 8  : ', end='')
            if lampam_1[7] > 0:
                print(' ', end='')
            print(f'{lampam_1[7] :.5f} : ', end='')
            if lampam_2[7] > 0:
                print(' ', end='')
            print(f'{lampam_2[7]:.5f}')

            print('\x1b[0m', end='')

            print('lampam 9  : ', end='')
            if lampam_1[8] > 0:
                print(' ', end='')
            print(f'{lampam_1[8] :.5f} : ', end='')
            if lampam_2[8] > 0:
                print(' ', end='')
            print(f'{lampam_2[8]:.5f}')

            print('lampam 10 : ', end='')
            if lampam_1[9] > 0:
                print(' ', end='')
            print(f'{lampam_1[9] :.5f} : ', end='')
            if lampam_2[9] > 0:
                print(' ', end='')
            print(f'{lampam_2[9]:.5f}')

            print('lampam 11 : ', end='')
            if lampam_1[10] > 0:
                print(' ', end='')
            print(f'{lampam_1[10] :.5f} : ', end='')
            if lampam_2[10] > 0:
                print(' ', end='')
            print(f'{lampam_2[10]:.5f}')

            print('lampam 12 : ', end='')
            if lampam_1[11] > 0:
                print(' ', end='')
            print(f'{lampam_1[11] :.5f} : ', end='')
            if lampam_2[11] > 0:
                print(' ', end='')
            print(f'{lampam_2[11]:.5f}')
        else:
            print(f'lampam 1 :{lampam_1[0]:.5f}:{lampam_2[0]:.5f}:{lampam_2[0] - lampam_1[0]:.5f}')
            print(f'lampam 2 :{lampam_1[1]:.5f}:{lampam_2[1]:.5f}:{lampam_2[1] - lampam_1[1]:.5f}')
            print(f'lampam 3 :{lampam_1[2]:.5f}:{lampam_2[2]:.5f}:{lampam_2[2] - lampam_1[2]:.5f}')
            print(f'lampam 4 :{lampam_1[3]:.5f}:{lampam_2[3]:.5f}:{lampam_2[3] - lampam_1[3]:.5f}')
            print('\033[37m', end='')
            print(f'lampam 5 :{lampam_1[4]:.5f}:{lampam_2[4]:.5f}:{lampam_2[4] - lampam_1[4]:.5f}')
            print(f'lampam 6 :{lampam_1[5]:.5f}:{lampam_2[5]:.5f}:{lampam_2[5] - lampam_1[5]:.5f}')
            print(f'lampam 7 :{lampam_1[6]:.5f}:{lampam_2[6]:.5f}:{lampam_2[6] - lampam_1[6]:.5f}')
            print(f'lampam 8 :{lampam_1[7]:.5f}:{lampam_2[7]:.5f}:{lampam_2[7] - lampam_1[7]:.5f}')
            print('\x1b[0m', end='')
            print(f'lampam 9 :{lampam_1[8]:.5f}:{lampam_2[8]:.5f}:{lampam_2[8] - lampam_1[8]:.5f}')
            print(f'lampam 10:{lampam_1[9]:.5f}:{lampam_2[9]:.5f}:{lampam_2[9] - lampam_1[9]:.5f}')
            print(f'lampam 11:{lampam_1[10]:.5f}:{lampam_2[10]:.5f}:{lampam_2[10] -lampam_1[10]:.5f}')
            print(f'lampam 12:{lampam_1[11]:.5f}:{lampam_2[11]:.5f}:{lampam_2[11] -lampam_1[11]:.5f}')
            print()


def print_ss(ss1, elem_per_line=200):
    " Prints all the components of a stacking sequence"
    ss1 = np.copy(ss1)
    ss1.reshape((ss1.size,))
    for ind in range(ss1.size):
        if ind % elem_per_line == 0:
            print(f'{ss1[ind]:3d}, ', end='')
        elif ind % elem_per_line == elem_per_line - 1:
            print(f'{ss1[ind]:3d}, ')
        else:
            print(f'{ss1[ind]:3d}, ', end='')
    print(' ')



def print_list(ss1, elem_per_line=200):
    "Prints all the components of a list of floats"
    if isinstance(ss1, list):
        length = len(ss1)
    else:
        length = ss1.size
    for ind in range(length):
        if ind % elem_per_line == 0:
            print(f'{ss1[ind]:0.5f}, ', end='')
        elif ind % elem_per_line == elem_per_line - 1:
            print(f'{ss1[ind]:0.5f}, ')
        else:
            print(f'{ss1[ind]:0.5f}, ', end='')
    print(' ')
    print(' ')


def print_list_ss(ss1, elem_per_line=200):
    #print(type(ss1), type(elem_per_line))
    " Prints all the components of a list of stacking sequences"
    if isinstance(ss1, list):
        length = len(ss1)
    else:
        length = ss1.shape[0]
    for iii in range(length):
        print_ss(ss1[iii], elem_per_line)
    print(' ')

if __name__ == "__main__":
    print('\n*** Test for the function print_lampam ***')
    print_lampam(np.arange(12))
    print('\n*** Test for the function print_ss ***')
    print_ss(np.arange(-20, 20), elem_per_line=40)
    print('\n*** Test for the function print_list ***')
    print_list(np.arange(20), elem_per_line=5)
    print('\n*** Test for the function print_list_ss ***')
    print_list_ss([np.arange(20), -np.arange(20)], elem_per_line=5)
