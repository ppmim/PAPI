#! /usr/bin/env python

# Copyright (c) 2009 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC

# This file is part of IAAraf.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


def prefix():
    """ Prefix that will be added to each line printed to the standard output.

    This method is used as a global variable, so that the string to be prefixed
    to each line printed to standard output can be easily modified. If, for
    example, this method is set to return ">> ", the following line...
    print prefix() + "Dayvan Cowboy, by Boards of Canada."
    ... would display ">> Dayvan Cowboy, by Boards of Canada."

    """

    return ">> "


def error_exit_message():
    """ Return the error message to be printed if the execution is aborted. """

    return prefix() + "Execution aborted."


def red_string(string):
    """ Return the string in red color, using the ANSI escape sequences.
    
    Keyword arguments:
    string - string to be printed in red color

    """

    # These codes are used as seen here:
    # (1) http://www.justlinux.com/forum/archive/index.php/t-107493.html
    # (2) http://stackoverflow.com/questions/287871/print-in-terminal-with-\
    # colors-using-python

    reset = "\033[0;0m"    
    CSI = "\x1B["
    return CSI + "31;40m" + string + reset


def bold_string(string):
    """ Return the string in bold, using the ANSI escape sequences.
    
    Keyword arguments:
    string - string to be printed in bold

    """
   
    bold = "\033[1m"
    reset = "\033[0;0m"   
    return bold + string + reset
