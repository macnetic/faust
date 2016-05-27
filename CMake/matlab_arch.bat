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
