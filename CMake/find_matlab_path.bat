@echo off
set "find_exe=matlab.exe"

(where matlab.exe) > logPath.txt 
(where /R "C:\\Program Files\\MATLAB" matlab.exe) >> logPath.txt 
(where /R "C:\\Program Files (x86)\\MATLAB" matlab.exe) >> logPath.txt


set MATLAB_EXE_DIR_TMP=
set /p MATLAB_EXE_DIR_TMP=<logPath.txt

REM for /f "delims=" %%i in ('type logPath.txt') do (set MATLAB_DIR_TMP=%%i && echo %%i)

echo environment variable is defined for matlab Path : %MATLAB_EXE_DIR_TMP%
