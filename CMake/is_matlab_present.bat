@echo off
rem set "find_exe=jgtjio.exe"
set "find_exe=matlab.exe"
setlocal enableDelayedExpansion

(for /f "delims=" %%i in ('where %find_exe%') do set MATLAB_DIR_TMP=%%i)>NUL 2>&1
rem echo %MATLAB_DIR_TMP%

rem set MATLAB_DIR_TMP="fbfbezbf"
rem echo MATLAB : %MATLAB_DIR_TMP%


if "%MATLAB_DIR_TMP%" == "" (
	rem echo VIDE
	echo|set /p=0
) else (
	
	rem echo PAS VIDE
	rem enleve les espaces
	set  "MATLAB_DIR_TMP=!MATLAB_DIR_TMP: =!"
	rem enleve bin
	set "MATLAB_DIR_WITHOUT_BIN=!MATLAB_DIR_TMP:\bin\=!"
	
	rem echo MATLAB = !MATLAB_DIR_TMP!
	rem echo MATLAB_without BIN = !MATLAB_DIR_WITHOUT_BIN!
	if "!MATLAB_DIR_TMP!" == "!MATLAB_DIR_WITHOUT_BIN!" (
	echo|set /p=0
	) else (
	echo|set /p=1
	rem echo bin present
	)

)
@echo on
