import sys
import os

from mpmath import *
mp.prec = 500

def full_sum_length(t):
     return int(sqrt(t/(2*pi)))

command_template = 'zetacalc --t={t} --start={start} --length={length} --N={N} --delta={delta} --number_of_threads=0 --terse'

def commands_for_split_computation(t):
    if t < 1e20:
        num_splits = 1
    elif t < 1e22:
        num_splits = 10
    elif t < 1e26:
        num_splits = 100
    elif t < 1e28:
        num_splits = 1000
    else:
        num_splits = 10000

    full_length = full_sum_length(t)
    split_length = full_length/num_splits
    start = 1
    for k in range(num_splits - 1):
        command = command_template.format(t=t, start=start, length=split_length, N=1000, delta=.04)
        print command
        start = start + split_length

    final_length = full_length - start + 1
    command = command_template.format(t=t, start=start, length=final_length, N=1000, delta=.04)
    print command



if __name__ == '__main__':
    t = mpf(sys.argv[1])
    commands_for_split_computation(t)
