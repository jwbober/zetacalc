import sys
import os
import os.path
import fcntl
import time
import subprocess

def main():
    work_location = sys.argv[1]
    print work_location

    while 1:
        lockfile = open(os.path.join(work_location, 'lockfile'), "w")
        fcntl.lockf(lockfile, fcntl.LOCK_EX)

        possible_work_locations = os.listdir(work_location)
        possible_work_locations.sort()

        found_work = False

        for filename in possible_work_locations:
            if not os.path.isdir(os.path.join(work_location, filename)):
                continue
            if not os.path.exists( os.path.join(work_location, filename, "todo") ):
                continue
            possible_work_units = os.listdir( os.path.join(work_location, filename, "todo") )
            if len(possible_work_units) == 0:
                continue
            found_work = True
            possible_work_units.sort()
            work_unit_original_location = os.path.join(work_location, filename, "todo", possible_work_units[0])
            work_unit_new_location = os.path.join(work_location, filename, "inprogress", possible_work_units[0])
            print "beginning work on", work_unit_new_location
            os.rename(work_unit_original_location, work_unit_new_location)

            break

        lockfile.close()

        if found_work:
            return_code = subprocess.call(['./zeta', work_unit_new_location])
            if(return_code != 0):
                print "problem calling ./zeta, returning work unit to 'todo'."
                os.rename(work_unit_new_location, work_unit_original_location)
            else:
                print "successful, removing 'inprogress' file"
                os.remove(work_unit_new_location)
        else:
            print "Found no work to do, sleeping for a minute"
            time.sleep(60)


if __name__ == '__main__':
    main()
