@echo off

for /f "delims=" %%i in ('where matlab.exe') do set MATLAB_DIR_TMP=%%i

setlocal enableDelayedExpansion
set ^"str=%MATLAB_DIR_TMP%"
for /F "delims=" %%a in (^"!str:\bin^=^

!^") do (
set "MATLAB_ROOT_DIR=%%a"
goto :break
)
:break

echo|set /p=%MATLAB_ROOT_DIR%
