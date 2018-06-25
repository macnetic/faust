::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::                              Description:                                ::::
::::  batch script to run matlab test. Take 2 arguments :                     ::::
::::   -1st arg is a char representing the list of matlab command that will   :::: 
::::      be executed														  ::::
::::   -2nd arg is the name of the text file 								  ::::
::::	where all the display in the matlab command window will be written in ::::
::::    in order to display this file into the current terminal				  ::::	
::::																		  ::::	 
::::  For more information on the FAuST Project, please visit the website     ::::
::::  of the project : <http://faust.inria.fr>                         ::::
::::                                                                          ::::
::::                              License:                                    ::::
::::  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ::::
::::                      Luc Le Magoarou, Remi Gribonval                     ::::
::::                      INRIA Rennes, FRANCE                                ::::
::::                      http://www.inria.fr/                                ::::
::::                                                                          ::::
::::  The FAuST Toolbox is distributed under the terms of the GNU Affero      ::::
::::  General Public License.                                                 ::::
::::  This program is free software: you can redistribute it and/or modify    ::::
::::  it under the terms of the GNU Affero General Public License as          ::::
::::  published by the Free Software Foundation.                              ::::
::::                                                                          ::::
::::  This program is distributed in the hope that it will be useful, but     ::::
::::  WITHOUT ANY WARRANTY; without even the implied warranty of              ::::
::::  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ::::
::::  See the GNU Affero General Public License for more details.             ::::
::::                                                                          ::::
::::  You should have received a copy of the GNU Affero General Public        ::::
::::  License along with this program.                                        ::::
::::  If not, see <http://www.gnu.org/licenses/>.                             ::::
::::                                                                          ::::
::::                             Contacts:                                    ::::
::::      Nicolas Bellot  : nicolas.bellot@inria.fr                           ::::
::::      Adrien Leman    : adrien.leman@inria.fr                             ::::
::::      Thomas Gautrais : thomas.gautrais@inria.fr                          ::::
::::      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ::::
::::      Remi Gribonval  : remi.gribonval@inria.fr                           ::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
@ECHO OFF

:: par defaut le test est bon
set ERRORLEVEL=0
 
set matlab_command_core=%1
set outputfile=%2%.txt
::transform the "/" character into "\" to have a valid path
set outputfile=%outputfile:/=\%

echo "OUTPUTFILE %outputfile%"

set "matlab_command_begin= delete %outputfile%;diary %outputfile%;try;testpass=0;"
set "matlab_command_end= catch ME ;testpass=2;disp(getReport(ME));end; diary OFF;exit(testpass);"
set matlab_command=%matlab_command_begin%%matlab_command_core%%matlab_command_end%
::echo DEBUT : %debut%
::echo MIDLLE : %middle%
::echo FIN : %fin%
::echo TEST : %test%

::transform the "@" character into "," 




set matlab_command=%matlab_command:@=,%

::transform the "$" character into ";"
set matlab_command=%matlab_command:$=;%


echo MATLAB COMMAND : %matlab_command%

:: launch the test in matlab an wait for matlab to finish
matlab -nodesktop -wait -r """%matlab_command%"""

:: si le code renvoyer par la derniere cmd est sup a 2 => erreur
IF ERRORLEVEL 2 set ERRORLEVEL=2




:: display the content of the MATLAB Command Window stored in a txtfile
echo ****** MATLAB WINDOW ******* 
type %outputfile%
IF ERRORLEVEL 1 set ERRORLEVEL=1

ECHO ERROR_LEVEL %ERRORLEVEL%
::renvoie le code erreur
exit %ERRORLEVEL%
