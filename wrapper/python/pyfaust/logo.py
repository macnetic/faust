import matplotlib.pyplot as plt
from random import random
from tempfile import gettempdir
from os.path import join
from pyfaust import Faust
import numpy as np

"""
This module is for generating the FAµST logo as .svg or simply as a
pyfaust.Faust object.
"""

plt.rcParams.update({
        "text.usetex": True,
#        "font.family": "Helvetica"
        "font.family": "stixsans"

})

def gen_line(array, x1, y1, x2, y2, value='rand', thickness=1,
              shift_links=True, shift=1):
    """
    Writes horizontal or vertical line in array.

    Args:
        array: output array.
        x1: start x-coord in array rows.
        x2: end x-coord in array rows.
        y1: start x-coord in array columns.
        y2: end x-coord in array columns.
        value: 'max' or 'rand' > 0 (might also be seen as the color).
        thickness: thickness of the line.
        shift_links: True for an alternate shift of the links that form the line.
        shift: the size of the shift set by shift_links.
    """
    x1, y1, x2, y2 = int(x1), int(y1), int(x2), int(y2)
    if x2 < x1:
        x1, x2 = x2, x1
    if y2 < y1:
        y1, y2 = y2, y1
    for x in [x1, x2]:
        if x < 0 or x >= array.shape[0]:
            raise ValueError('Overflow array row dimension')

    for y in [y1, y2]:
        if y < 0 or y >= array.shape[1]:
            raise ValueError('Overflow array column dimension')

    val = 1 if value == 'max' else (random() * .75 + .25)
    if x1 == x2 and y2 != y1:
        # draw vertical row
        n = y2 - y1 + 1
        if shift_links:
            array[x1+shift:thickness+x1+shift, y1:y1+n:2] = val
            array[x1:thickness+x1, y1+1:y1+n:2] = val
        else:
            array[x1:x1 + n, y1] = val
    elif y1 == y2 and x1 != x2:
        # draw horizontal row
        n = x2 - x1 + 1
        if shift_links:
            array[x1:x1+n:2, y1+shift:thickness+y1+shift] = val
            array[x1+1:x1+n:2, y1:thickness+y1] = val
        else:
            array[x1, y1:y1 + n] = val
    else:
        raise ValueError('can\'t draw anything but horizontal/vertical rows')

def gen_F(base_size, margin, fontsize):
    """
    Generates array of letter F.
    """
    F = np.zeros((base_size, base_size))
    abs_margin = int(base_size * margin)
    x1 = int(abs_margin)
    x2 = base_size - int(abs_margin) - 1
    y1 = y2 = abs_margin
    gen_line(F, x1, y1, x2, y2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    gen_line(F, y1, x1, y2, x2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    mid_offset = base_size // 2 - abs_margin
    gen_line(F, y1 + mid_offset, x1, y2 + mid_offset, x2 - mid_offset // 2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    return F

def gen_A(base_size, margin, fontsize):
    """
    Generates array of letter A.
    """
    A = np. zeros((base_size, base_size))
    abs_margin = int(base_size * margin)
    x1 = int(abs_margin)
    x2 = base_size - int(abs_margin) - 1
    y1 = y2 = abs_margin
    # left vert line
    gen_line(A, x1, y1, x2, y2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    right_offset = base_size - y1 - abs_margin
    # right vertical line
    gen_line(A, x1, right_offset, x2, right_offset, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    # up horizontal line
    gen_line(A, y1, x1, y2, x2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    mid_offset = base_size // 2 - abs_margin
    # middle horizontal line
    gen_line(A, y1 + mid_offset, x1, y2 + mid_offset, x2 - mid_offset // 2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    return A

def gen_MU(base_size, margin, fontsize):
    """
    Generates array of letter \MU.
    """
    MU = np. zeros((base_size, base_size))
    abs_margin = int(base_size * margin)
    x1 = int(abs_margin)
    x2 = base_size - int(abs_margin) - 1
    y1 = y2 = abs_margin
    # left vert line
    gen_line(MU, x1 + abs_margin, y1, x2 + abs_margin, y2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    right_offset = base_size - y1 - abs_margin
    # right vertical line
    gen_line(MU, x1+ abs_margin, right_offset, x2 - abs_margin+ abs_margin, right_offset, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    mid_offset = base_size // 2 - abs_margin
    # middle horizontal line
    gen_line(MU, y1 + mid_offset+ abs_margin, x1, y2 + mid_offset+ abs_margin, x2 - abs_margin, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    return MU

def gen_S(base_size, margin, fontsize):
    """
    Generates array of letter S.
    """
    S = np. zeros((base_size, base_size))
    abs_margin = int(base_size * margin)
    x1 = int(abs_margin)
    x2 = base_size - int(abs_margin) - 1
    y1 = y2 = abs_margin
    mid_offset = base_size // 2 - abs_margin
    # left vertical line
    gen_line(S, x1, y1, x1 + mid_offset, y2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    # up horizontal line
    gen_line(S, y1, x1, y2, x2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    # middle horizontal line
    gen_line(S, y1 + mid_offset - abs_margin // 2, x1, y2 + mid_offset -
              abs_margin // 2, x2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    # right vertical line
    gen_line(S, x1 + mid_offset - abs_margin // 2, base_size - y1 - int(abs_margin * 5 / 4),
              x1 + 2 * mid_offset - abs_margin // 2, base_size - y2 - int(abs_margin * 5 / 4),
              shift_links=True, thickness=abs_margin, shift=abs_margin//2)
    # bottom horizontal line
    gen_line(S, base_size - y1 - abs_margin, x1, base_size - y2 -
              abs_margin, x2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    return S

def gen_T(base_size, margin, fontsize):
    T = np. zeros((base_size, base_size))
    abs_margin = int(base_size * margin)
    x1 = int(abs_margin)
    x2 = base_size - int(abs_margin) - 1
    y1 = y2 = abs_margin
    mid_offset = int(base_size / 2 - abs_margin / 2)
    # middle vert line
    gen_line(T, x1, mid_offset, x2, mid_offset, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    right_offset = base_size - y1 - abs_margin
    # up horizontal line
    gen_line(T, y1, x1, y2, x2, shift_links=True, thickness=abs_margin,
              shift=abs_margin//2)
    return T

def draw_F(base_size, margin, fontsize):
    F = gen_F(base_size, margin, fontsize)
    plt.xlabel(r'\textbf{\textit{Flexible}}', fontweight='bold', fontsize=fontsize)
    plt.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)
    plt.tick_params(
        axis='y',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelleft=False)

    plt.imshow(F, aspect='equal')

def draw_A(base_size, margin, fontsize):
    A = gen_A(base_size, margin, fontsize)
    plt.xlabel(r'\textbf{\textit{Approximate}}', fontweight='bold', fontsize=fontsize,
           fontname='stixsans')
    plt.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)
    plt.tick_params(
        axis='y',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelleft=False)

    plt.imshow(A, aspect='equal')

def draw_MU(base_size, margin, fontsize):
    MU = gen_MU(base_size, margin, fontsize)
    plt.xlabel(r'\textbf{\textit{MUltilayer}}', fontweight='bold', fontsize=fontsize,
           fontname='stixsans')
    plt.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)
    plt.tick_params(
        axis='y',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelleft=False)
    plt.imshow(MU, aspect='equal')


def draw_S(base_size, margin, fontsize):
    S = gen_S(base_size, margin, fontsize)
    plt.xlabel(r'\textbf{\textit{Sparse}}', fontweight='bold', fontsize=fontsize,
           fontname='stixsans')
    plt.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)
    plt.tick_params(
        axis='y',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelleft=False)

    plt.imshow(S, aspect='equal')

def draw_T(base_size, margin, fontsize):
    """
    Generates array of letter T.
    """
    T = gen_T(base_size, margin, fontsize)
    plt.xlabel(r'\textbf{\textit{Transforms}}', fontweight='bold', fontsize=fontsize,
           fontname='stixsans')
    plt.tick_params(
            axis='x',
            which='both',
            bottom=False,
            top=False,
            labelbottom=False)
    plt.tick_params(
        axis='y',
        which='both',
        bottom=False,
        top=False,
        left=False,
        labelleft=False)

    plt.imshow(T, aspect='equal')

def gen_faust_logo(base_size=50, margin=.15, fontsize=12):
    a = [base_size, margin, fontsize]
    F, A, MU, S, T = gen_F(*a), gen_A(*a), gen_MU(*a), gen_S(*a), gen_T(*a)
    logoF = Faust([F, A, MU, S, T])
    return logoF

def draw_faust(base_size=50, margin=.15, fontsize=12):
    """
    Outputs faust_logo.svg (with margins) and faust_logo-tight.svg (without
    margins) in /tmp directory.

    Args:
        base_size: width and height of the letter array to draw.
        fontsize: the font size to be used to write acronym.
        margin: in [0, 1].
    """
    tmpdir = gettempdir()
    plt.subplot(161)
    draw_F(base_size, margin, fontsize)
    plt.subplot(162)
    draw_A(base_size, margin, fontsize)
    plt.subplot(163)
    draw_MU(base_size, margin, fontsize)
    plt.subplot(164)
    draw_S(base_size, margin, fontsize)
    plt.subplot(165)
    draw_T(base_size, margin, fontsize)
    plt.tight_layout()
    plt.savefig(join(tmpdir, 'faust_logo.svg'), transparent=True)
    plt.savefig(join(tmpdir, 'faust_logo-tight.svg'), bbox_inches='tight', transparent=True)


if __name__ == '__main__':
    draw_faust()
    logo_F = gen_faust_logo()
    logo_F.imshow()
    print("logo_F:", logo_F)
    plt.show()
