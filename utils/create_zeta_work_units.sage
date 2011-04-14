import sys
import os
import fcntl

smallest_work_unit_size = 2500000000
default_number_of_work_units = 2000

def main():
    t = RealField(200)(sys.argv[1])
    N = sys.argv[2]
    delta = sys.argv[3]
    output_location = sys.argv[4]
    output_name = sys.argv[5] + "_" + sys.argv[1]

    todo_location = os.path.abspath(os.path.join(output_location, output_name, 'todo'))
    inprogress_location = os.path.abspath(os.path.join(output_location, output_name, "inprogress"))
    finished_location = os.path.abspath(os.path.join(output_location, output_name, "finished"))

    lockfile = open(os.path.join(output_location, "lockfile"), "w")
    fcntl.lockf(lockfile, fcntl.LOCK_EX)

    if not os.path.exists( os.path.join(output_location, output_name) ):
        os.mkdir( os.path.join(output_location, output_name))
    if not os.path.exists(todo_location):
        os.mkdir(todo_location)
    if not os.path.exists(finished_location):
        os.mkdir(finished_location)
    if not os.path.exists(inprogress_location):
        os.mkdir(inprogress_location)


    endpoint = floor(sqrt(t/(2 * pi)))

    work_unit_size = floor(endpoint/default_number_of_work_units)

    work_unit_size = max(work_unit_size, smallest_work_unit_size)

    if endpoint % work_unit_size == 0:
        number_of_work_units = endpoint/work_unit_size
    else:
        number_of_work_units = floor(endpoint/work_unit_size) + 1

    count = 1
    while count < number_of_work_units:
        start = (count - 1) * work_unit_size + 1
        filename = os.path.join(todo_location, "work_unit%04d" % count)
        print "writing work unit to", filename
        outfile = open(filename, 'w')
        outfile.write(t.str(truncate=False, no_sci = 2, skip_zeroes = True))
        outfile.write(" ")
        outfile.write(start.str())
        outfile.write(" ")
        outfile.write(work_unit_size.str())
        outfile.write(" ")
        outfile.write(N)
        outfile.write(" ")
        outfile.write(delta)
        outfile.write(" ")
        outfile.write(os.path.join(finished_location, "work_unit%04d" % count))
        outfile.close()
        count = count + 1

    # there is still one more work unit to write
    # its length may be different than the others

    final_work_unit_size = endpoint - (count - 1) * work_unit_size

    start = (count - 1) * work_unit_size + 1
    filename = os.path.join(todo_location, "work_unit%04d" % count)
    outfile = open(filename, 'w')
    print "writing work unit to", filename
    outfile.write(t.str(truncate=False, no_sci = 2, skip_zeroes = True))
    outfile.write(" ")
    outfile.write(start.str())
    outfile.write(" ")
    outfile.write(final_work_unit_size.str())
    outfile.write(" ")
    outfile.write(N)
    outfile.write(" ")
    outfile.write(delta)
    outfile.write(" ")
    outfile.write(os.path.join(finished_location, "work_unit%04d" % count))
    outfile.close()

    lockfile.close()


if __name__=="__main__":
    main()
