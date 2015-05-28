#! /usr/bin/env python

# Author: Victor Terron (c) 2014
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

import os.path
import pipes
import subprocess

def check_open(path, process):
    """ 
    Check whether the file is open by a process.

    https://stackoverflow.com/q/17321930/184363.

    """

    proc  = pipes.quote(process)
    file_ = pipes.quote(os.path.basename(path))

    # Find the files currently opened by the process (or processes) with the
    # specified name. Then, for each process, list its open files and check
    # whether one of them matches 'path'.

    args = ['pgrep', proc]
    try:
        output = subprocess.check_output(args)
    except subprocess.CalledProcessError:
        # There is no process
        return False

    PIDs = (int(line) for line in output.splitlines())

    for pid_ in PIDs:
        cmd = "ls -l /proc/{0}/fd | grep {1}".format(pid_, file_)
        try:
            # grep returns zero if selected lines are found, so check_output()
            # will raise an exception except when there is a match.
            subprocess.check_output(cmd, shell=True)
            return True
        except subprocess.CalledProcessError:
            pass
    # else suite is executed after the for, but only if the for terminates 
    # normally (not by a break).
    else:
        return False


if __name__ == "__main__":

    # Example:
    import tempfile
    import os
    fd, path = tempfile.mkstemp()
    os.close(fd)

    print check_open(path, 'python')
    with open(path, "wt") as fd:
        print check_open(path, 'python')
    print check_open(path, 'python')
