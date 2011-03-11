import sys

file1 = open(sys.argv[1], 'r')
file2 = open(sys.argv[2], 'r')

contents1 = file1.read().split('\n')
contents2 = file2.read().split('\n')

for n in range(1, len(contents1)):
    a1 = contents1[n].strip()[1:-1]
    a2 = contents2[n].strip()[1:-1]
 
    if(a1 == ""):
        continue

    x1, y1 = a1.split(',')
    x1 = float(x1)
    y1 = float(y1)


    x2, y2 = a2.split(',')
    x2 = float(x2)
    y2 = float(y2)

    print x1, y1, x1 - x2, y1 - y2
