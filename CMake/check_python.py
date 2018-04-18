##############################################################################
##                              Description:                                ##
##       script to check if all the Python needed module are installed      ## 
##                i.e numpy  and cython                                     ##
##                                                                          ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.gforge.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##############################################################################


import sys
#getting the version of python
output_value=0

print("python version : "+str(sys.version_info.major)+'.'+str(sys.version_info.minor))



###### NECESSARY MODULE ######
print("*** PYTHON NECESSARY MODULES ***")
try:
    print('NUMPY PYTHON module :')
    import numpy 
    print('module installed')
    print("version : "+numpy.__version__)
    if (numpy.__version__ < '1.9.2'):
        print("WARNING your numpy version may be too old, only tested with 1.9.2 version")
except ImportError:
    print('module is missing !')
    output_value=-1
    
try:
    print('CYTHON PYTHON module : ')
    import Cython  
    print('module installed')
    print("version : "+Cython.__version__)
except ImportError:
    print('module is missing !')
    output_value=-1


###### OPTIONAL MODULE ######
print("*** PYTHON OPTIONAL MODULES ***")
try:
    print('SCIPY PYTHON module :')
    import scipy as sp
    print('module installed')
    print("version : "+sp.__version__)

except ImportError:
    print('module Scipy not present, no time comparison with scipy will be made')
    output_value=1
    
try:
    print('MATPLOTLIB PYTHON module :')
    import matplotlib as plt
    print('module installed')
    print("version : "+plt.__version__)

except ImportError:
    print('module matplotlib not present, no display will be made')
    
    
exit(output_value)
