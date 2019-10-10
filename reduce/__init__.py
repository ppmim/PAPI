# Copyright (c) 2009-2012 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''PANIC data processing system (PAPI)'''

import logging

__version__ = '1.2.0'

# Top level NullHandler
logging.getLogger("papi").addHandler(logging.NullHandler())

# modules
from .calDark import *
from .calDarkModel import *
from .threadsmod import *
from .calTwFlat import *
from .applyDarkFlat import *
from .calGainMap import *
from .calSuperFlat import *
from .calDomeFlat import *
from .dxtalk import *
from .eval_focus_serie import *
