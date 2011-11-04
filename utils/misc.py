import sys
import os
import subprocess

location = "/home/bober/math/zeta_computations/rs_sums/"

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

def all_zeros():
    output_location = "/home/bober/zeta_computations/"
    files = os.listdir(location)
    files.sort()
    for filename in files:
        print "processing", filename
        sys.stdout.flush()
        output_file = open(output_location + "zeros/" + filename[:-5] + "zeros", 'w')
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--zeros',
                                          '--terse' ],
                                        stdout=output_file)
        output_file.close()

def write_S_values():
    output_location = "/home/bober/math/zeta_computations/"
    files = os.listdir(location)
    files.sort()
    for filename in files:
        print "processing", filename
        sys.stdout.flush()
        output_file = open(output_location + "S/" + filename[:-5] + "Svalues", 'w')
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--S',
                                          '--terse' ],
                                        stdout=output_file)
        output_file.close()

def write_Z_values():
    output_location = "/home/bober/math/zeta_computations/"
    files = os.listdir(location)
    files.sort()
    for filename in files:
        print "processing", filename
        sys.stdout.flush()
        output_file = open(output_location + "values/" + filename[:-5] + "graph", 'w')
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--values',
                                          '--spacing', '.0004' ],
                                        stdout=output_file)
        output_file.close()

def write_max_info():
    output_location = "/home/bober/math/zeta_computations/"
    files = os.listdir(location)
    files.sort()
    output_file = open(output_location + "max_values_temp", 'a')
    for filename in files:
        print "processing", filename
        sys.stdout.flush()
        return_value = subprocess.call([  'blfi',
                                          '--filename', location + filename,
                                          '--maxmin',
                                          '--terse'],
                                        stdout=output_file)
        sys.stdout.flush()

    outfile_file.close()
if __name__ == "__main__":
    #write_max_info()
    #write_Z_values()
    print_largest_S_values()
