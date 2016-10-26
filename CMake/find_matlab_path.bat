::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::                              Description:                                ::::
::::  batch script that will find where matlab is on your computer			  ::::
::::                                                                          ::::
::::  For more information on the FAuST Project, please visit the website     ::::
::::  of the project : <http://faust.gforge.inria.fr>                         ::::
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

@echo off
set "find_exe=matlab.exe"

IF NOT exist "tmp\*.*"  mkdir tmp

(where matlab.exe) > tmp/tmpPathMatlab.txt 
(where /R "C:\\Program Files\\MATLAB" matlab.exe) >> tmp\\tmpPathMatlab.txt 
(where /R "C:\\Program Files (x86)\\MATLAB" matlab.exe) >> tmp\\tmpPathMatlab.txt


REM set MATLAB_EXE_DIR_TMP=
REM set /p MATLAB_EXE_DIR_TMP=<tmpPathMatlab.txt

REM set /a cpt=1
REM for /f "tokens=*" %%a In (tmpPathMatlab.txt) do (
REM echo --
REM echo %%a
REM set MATLAB_EXE_DIR_TMP%cpt%=%%a
REM set /a cpt=%cpt%+1
REM echo -- 
REM )
REM echo environment variable is defined for matlab Path : %MATLAB_EXE_DIR_TMP%
