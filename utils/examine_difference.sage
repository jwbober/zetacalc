import sys

def examine_difference(filename1, filename2):
    file1 = open(filename1)
    file2 = open(filename2)

    contents1 = file1.read().split('\n')
    contents2 = file2.read().split('\n')

    for n in range(1, len(contents1)):
        z1 = CC(contents1[n])
        z2 = CC(contents2[n])

        print z1, z1 - z2

if __name__ == "__main__":
    examine_difference(sys.argv[1], sys.argv[2])
