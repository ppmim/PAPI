from __future__ import print_function
import os
import time

major = 1
minor = 2

rlfile = 'version.py'
backup = 'version.py.bak'
path = os.path.dirname(os.path.abspath(__file__))

def make_release():

    try:
        release = time.strftime("%Y%m%d%H%M%S", time.gmtime(time.time()))
        release = int(release)

        cwd = os.getcwd()
        os.chdir(path)

        if os.path.exists(backup):
            os.remove(backup)

        if os.path.exists(rlfile):
            os.rename(rlfile, backup)

        with open(rlfile, 'w') as out_f:
            out_f.write("# this file was automatically generated\n")
            out_f.write("major = %d\n" % major)
            out_f.write("minor = %d\n" % minor)
            out_f.write("release = %d\n" % release)
            out_f.write("\n")
            out_f.write("__version__ = '%d.%d.%d' % (major, minor, release)\n")
            out_f.write("\n")
    finally:
        os.chdir(cwd)
        return "%d.%d.%d" % (major, minor, release)
    
if __name__ == "__main__":
    print(make_release())
    
