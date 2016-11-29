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
print "python version : "+str(sys.version_info.major)+'.'+str(sys.version_info.minor)
output_value=0


###### NECESSARY MODULE ######
try:
    print('NUMPY PYTHON module :') 
    import numpy 
    print('module installed')
    print "version : "+numpy.__version__
except ImportError:
    print('module is missing !')
    output_value=-1
    
try:
    print('CYTHON PYTHON module : ')
    import Cython  
    print('module installed')
    print "version : "+Cython.__version__
except ImportError:
    print('module is missing !')
    output_value=-1


###### OPTIONAL MODULE ######
try:
    print('SCIPY PYTHON module :') 
    import scipy as sp
    print('module installed')
    print "version : "+sp.__version__  

except ImportError:
    print('module Scipy not present, no time comparison with scipy will be made')
    output_value=1
    
exit(output_value)
