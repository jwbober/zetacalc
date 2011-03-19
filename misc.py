import sys
import os
import subprocess

def rearrange_files():
    for filename in os.listdir("rs_sums"):
        f = open("rs_sums/" + filename)
        t = float(f.readline().strip())
        delta = float(f.readline().strip())

        f.close()
        if delta < .03:
            os.rename("rs_sums/" + filename, "short_ranges/" + filename)


def print_info():
    for filename in os.listdir("rs_sums"):
        f = open("rs_sums/" + filename)
        t = float(f.readline().strip())
        delta = float(f.readline().strip())

        f.close()
        print filename, t, delta

def check_RH():
    for filename in os.listdir("rs_sums"):
        return_value = subprocess.call([  'blfi',
                                          '--filename', 'rs_sums/' + filename,
                                          '--check_RH',
                                          '--terse' ])

        print filename, "\t", return_value

def print_min_midpoint_values():
    for filename in os.listdir("rs_sums"):
        print filename, "...\t",
        sys.stdout.flush()
        return_value = subprocess.call([  'blfi',
                                          '--filename', 'rs_sums/' + filename,
                                          '--min_midpoint_value',
                                          '--terse' ])
        sys.stdout.flush()

def print_largest_S_values():
    for filename in os.listdir("rs_sums"):
        print filename, "...\t",
        sys.stdout.flush()
        return_value = subprocess.call([  'blfi',
                                          '--filename', 'rs_sums/' + filename,
                                          '--largest_S_value',
                                          '--terse' ])
        sys.stdout.flush()



if __name__ == "__main__":
    check_RH()
