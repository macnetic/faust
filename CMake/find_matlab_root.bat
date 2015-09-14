@echo off
set "find_exe=matlab.exe"
set MATLAB_DIR_TMP=
(for /f "delims=" %%i in ('where %find_exe%') do set MATLAB_DIR_TMP=%%i)>NUL 2>&1 
rem echo %MATLAB_DIR_TMP%
setlocal enableDelayedExpansion
set ^"str=%MATLAB_DIR_TMP%"
for /F "delims=" %%a in (^"!str:\bin^=^

!^") do (
set "MATLAB_ROOT_DIR=%%a"
goto :break
)
:break

rem echo "MATLAB" %MATLAB_ROOT_DIR%
echo|set /p=%MATLAB_ROOT_DIR%
@echo on