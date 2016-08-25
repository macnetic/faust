@echo off
set "find_exe=matlab.exe"

(where matlab.exe) > .\tmp\tmpPathMatlab.txt 
(where /R "C:\\Program Files\\MATLAB" matlab.exe) >> .\tmp\tmpPathMatlab.txt 
(where /R "C:\\Program Files (x86)\\MATLAB" matlab.exe) >> .\tmp\tmpPathMatlab.txt

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
