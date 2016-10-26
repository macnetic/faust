::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::::                              Description:                                ::::
::::  batch file is useful to find the extension of the mexfile used by Matlab::::
::::  This script takes 2 arguments in input :								  :::: 	
::::  	-1st :  the number of the method that will be used					  ::::
::::				0 : use the Matlab command compute('arch')				  ::::
::::					 and write the result in a file specified by arg2     ::::
::::                1 : use the Matlab variable mexext                        ::::
::::					and write the result in a file specified by arg2	  ::::
:::: 				2 : use the name of a directory specified in arg2         :::: 
::::						 which is equal to MATLABROOT                     ::::
::::   -2nd : depending on the method used									  ::::
::::																		  ::::	 
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

if "%~1"=="" goto end
if "%~2"=="" goto end

if "%1"=="0" goto matlab_lib_subdir
if "%1"=="1" goto matlab_ext
if "%1"=="2" goto matlab_lib_subdir2

:matlab_lib_subdir
matlab -wait -noFigureWindows -nosplash -r """fid=fopen('%2','w');fprintf(fid,'%%s',computer('arch'));fclose(fid);exit"""
more %2
goto end 

:matlab_ext
matlab -wait -noFigureWindows -nosplash -r """fid=fopen('%2','w');fprintf(fid,'%%s',mexext);fclose(fid);exit"""
more %2
goto end 

:matlab_lib_subdir2
dir /B "%~2\extern\lib" | findstr /R /C:"^win[16|32|64]"
goto end 


:end

@echo on
