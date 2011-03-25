import sys
import os
import subprocess

location = "/home/bober/math/experiments/theta_sums2/data/rs_sums/"

def rearrange_files():
    files = os.listdir(location)
    files.sort()
    for filename in files:
        f = open(location + filename)
        t = float(f.readline().strip())
        delta = float(f.readline().strip())

        f.close()
        if delta < .03:
            os.rename(location + filename, "short_ranges/" + filename)


def print_info():
    files = os.listdir(location)
    files.sort()
    for filename in files:
        f = open(location + filename)
        t = float(f.readline().strip())
        delta = float(f.readline().strip())

        f.close()
        print filename, t, delta

def check_RH():
    files = os.listdir(location)
    files.sort()
    for filename in files:
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--check_RH',
                                          '--terse' ])

        print filename, "\t", return_value

def print_min_midpoint_values():
    files = os.listdir(location)
    files.sort()
    for filename in files:
        print filename, "...\t",
        sys.stdout.flush()
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--min_midpoint_value',
                                          '--terse' ])
        sys.stdout.flush()

def print_largest_S_values():
    files = os.listdir(location)
    files.sort()
    for filename in files:
        print filename, "...\t",
        sys.stdout.flush()
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--largest_S_value',
                                          '--terse' ])
        sys.stdout.flush()



if __name__ == "__main__":
    print_min_midpoint_values()
