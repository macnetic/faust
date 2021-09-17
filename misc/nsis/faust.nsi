Name "Faust-@CPACK_PACKAGE_VERSION@"


OutFile "faust-@CPACK_PACKAGE_VERSION@-amd64.exe"

; default install. dir.
InstallDir "$PROGRAMFILES64\Faust"

; admin rights
RequestExecutionLevel admin

; license content
LicenseText "Faust License (other secondary licenses will be installed here: $INSTDIR\doc)"
LicenseData "@PROJECT_BINARY_DIR@\..\license.txt"

;------------ Installer Pages
Page directory
Page License
Page instfiles

;------------------------------------- useful macro
; macro from http://nsis.sourceforge.net/StrRep
!define StrRep "!insertmacro StrRep"
!macro StrRep output string old new
    Push `${string}`
    Push `${old}`
    Push `${new}`
    !ifdef __UNINSTALL__
        Call un.StrRep
    !else
        Call StrRep
    !endif
    Pop ${output}
!macroend
 
!macro Func_StrRep un
    Function ${un}StrRep
        Exch $R2 ;new
        Exch 1
        Exch $R1 ;old
        Exch 2
        Exch $R0 ;string
        Push $R3
        Push $R4
        Push $R5
        Push $R6
        Push $R7
        Push $R8
        Push $R9
 
        StrCpy $R3 0
        StrLen $R4 $R1
        StrLen $R6 $R0
        StrLen $R9 $R2
        loop:
            StrCpy $R5 $R0 $R4 $R3
            StrCmp $R5 $R1 found
            StrCmp $R3 $R6 done
            IntOp $R3 $R3 + 1 ;move offset by 1 to check the next character
            Goto loop
        found:
            StrCpy $R5 $R0 $R3
            IntOp $R8 $R3 + $R4
            StrCpy $R7 $R0 "" $R8
            StrCpy $R0 $R5$R2$R7
            StrLen $R6 $R0
            IntOp $R3 $R3 + $R9 ;move offset by length of the replacement string
            Goto loop
        done:
 
        Pop $R9
        Pop $R8
        Pop $R7
        Pop $R6
        Pop $R5
        Pop $R4
        Pop $R3
        Push $R0
        Push $R1
        Pop $R0
        Pop $R1
        Pop $R0
        Pop $R2
        Exch $R1
    FunctionEnd
!macroend
!insertmacro Func_StrRep ""
;----------- Files to install

Section "" ; no component so name not needed

  CreateDirectory $INSTDIR\matlab
  CreateDirectory $INSTDIR\python
  CreateDirectory $INSTDIR\doc\html

  WriteUninstaller "$INSTDIR\uninstall.exe"

  ; install gpu mod lib
  CreateDirectory $INSTDIR\lib
  SetOutPath $INSTDIR\lib
  File @PROJECT_BINARY_DIR@\..\gpu_mod\build\gm.dll ; error if USE_GPU_MOD to OFF

  ; install python wrapper
  SetOutPath $INSTDIR\python\pyfaust
  File /r @PROJECT_BINARY_DIR@\wrapper\python\pyfaust\*py
  SetOutPath $INSTDIR\python
  File @PROJECT_BINARY_DIR@\wrapper\python\*pyd
  File @PROJECT_BINARY_DIR@\wrapper\python\*pxd

  ; install matlab wrapper
  SetOutPath $INSTDIR\matlab
  File /nonfatal /r /x old_matlab @PROJECT_BINARY_DIR@\wrapper\matlab\*.m @PROJECT_BINARY_DIR@\wrapper\matlab\*.@MEX_EXT@ @PROJECT_BINARY_DIR@\wrapper\matlab\*.mat
  ; nonfatal useful in case of data *.mat not used/present (because they are downloaded at installation)


  ; check the python version matches python major.minor version used to build the the wrapper shared library
  Exec "python --version | python -c $\"import re; ver = input(); exit(0) if re.match('Python @WIN_PY_VER@', ver) else exit(1)$\""
  IfErrors 0 +2
  MessageBox MB_OK "Error: this version of FAuST is pre-compiled for Python @WIN_PY_VER@ which must be installed and configured as the default python on your system (i.e. it must be available as $\"python$\" command in the PATH environment variable)." /SD IDOK IDOK data_dl
  MessageBox MB_OK "The pyfaust wrapper will be installed for Python @WIN_PY_VER@." /SD IDOK

  ; post install pyfaust auto-setup in environment (works only if python is installed in path)
  ${StrRep} '$0' $TEMP '\' '\\'
  Exec "python -c $\"import site;dir=site.getsitepackages()[1];f=open('$0\\tmp_site_pkg', 'w');f.write(dir);f.close()$\""
  IfErrors 0 +2
  MessageBox MB_OK "Error: no python found into your PATH environment variable. You'll have to do the Faust setup manually (you'll see how in the documentation)." /SD IDOK IDOK data_dl
  MessageBox MB_OK "Faust installed in your python environment (the version found into your PATH environment variable)." /SD IDOK

  FileOpen $1 "$TEMP\tmp_site_pkg" r
  FileRead $1 $2
  FileClose $1
  ;MessageBox MB_OK "$2"

  SetOutPath $2\pyfaust
  File /r @PROJECT_BINARY_DIR@\wrapper\python\pyfaust\*py
  SetOutPath $2\pyfaust\lib
  File @PROJECT_BINARY_DIR@\..\gpu_mod\build\gm.dll
  ; CreateShortCut pyfaust @PROJECT_BINARY_DIR@\wrapper\python\pyfaust ; tested and it can't be like a linux symlink
  SetOutPath $2
  File @PROJECT_BINARY_DIR@\wrapper\python\*pyd
  ;File @PROJECT_BINARY_DIR@\wrapper\python\*pxd

  ; add data path in __init__.py (both into site-packages and $INSTDIR)
  FileOpen $1 "$INSTDIR\pyfaust\__init__.py" a
  ; go at the end of file (but with append mode is that necessary ?)
  FileSeek $1 0 END
  ; add the install path
  FileWrite $1 "$\r$\n_NSI_INSTALL_PATH='$INSTDIR'"
  FileClose $1

  ; do the same thing on the copy
  FileOpen $1 "$2\pyfaust\__init__.py" a
  FileSeek $1 0 END
  FileWrite $1 "$\r$\n_NSI_INSTALL_PATH='$INSTDIR'"
  FileClose $1

  ExecWait "python -m pip install @PYFAUST_PYTHON_REQUIREMENTS@"
  ; ExecWait doesn't work with this command, if eventually pip install command fails, the user will be noticed when importing pyfaust
  ;IfErrors 0 +2
  ;MessageBox MB_OK "Error: failed partly or totally to install the pyfaust python packages through pip, please install them manually to get a workable pyfaust, list of packages: @PYFAUST_PYTHON_REQUIREMENTS@." IDOK data_dl

  ; =================================================

  data_dl:
  ; download data into matlab wrapper data folder
  ; create data folder
  ; try with python first
  ${StrRep} '$3' $INSTDIR '\' '\\'
  Exec "python -c $\"from os import mkdir; mkdir('$3\\matlab\\data')$\""
  ; download data
  ClearErrors ; in case the data dir was already existing
  Exec "python $2\pyfaust\datadl.py $\"$INSTDIR\matlab\data$\""
  IfErrors 0 after_data_dl
  MessageBox MB_OK "Downloading FAuST data with python seems to have failed (or maybe it's already done), now trying with powershell." /SD IDOK
  ClearErrors
  ExecWait "powershell -WindowStyle Hidden Invoke-WebRequest $\"@REMOTE_DATA_URL@/@REMOTE_DATA_FILE@$\" -O $\"$TEMP\@REMOTE_DATA_FILE@$\"" ; ExecWait because unzipping needs download finished
  Exec "powershell -WindowStyle Hidden Expand-Archive -Force $\"$TEMP\@REMOTE_DATA_FILE@$\" '$INSTDIR\matlab\data'" ; output folder data is created auto. ; simple quote used to avoid powershell to think there is two arguments when we meant one argument for the dest. path (double quote doesn't allow that).
  IfErrors 0 after_data_dl
  MessageBox MB_OK "Error downloading FAuST data (or maybe it's already done). You'll need to download manually (please check the documentation)." /SD IDOK IDOK after_data_dl

  after_data_dl:

  ; post install matfaust auto-setup
  ; find matlab in PATH, to get the value of matlabroot and edit the startup.m
  ExecWait "matlab -nojvm -r $\"fd = fopen('$TEMP\matlabroot','w'); fwrite(fd, matlabroot); fclose(fd); exit $\""
  ifErrors loc_matlab_man 0
  ; unfortunately ExecWait is not working, maybe because matlab spawns subprocesses
  ; so wait manually until the output file exists or a timeout is reached
  StrCpy $R0 0
  wait_matlab:
      ifFileExists "$TEMP\matlabroot" +4
    sleep 1000
    IntOp $R0 $R0 + 1
    StrCmp $R0 60 0 wait_matlab ; sleep 60 secs max
  ifFileExists "$TEMP\matlabroot" +1 loc_matlab_man
  FileOpen $1 "$TEMP\matlabroot" r
  ifErrors loc_matlab_man +1
  FileRead $1 "$R9"
  FileClose $1
  StrCmp $R9 "" loc_matlab_man
  ifFileExists "$R9" +1 loc_matlab_man
  Call matlabFoundCb
  MessageBox MB_OK 'Faust installed in $R9 (matlab was found in PATH environment variable).' /SD IDOK
  goto continue

  loc_matlab_man: ; manually locate matlab
  ; (in case matlab is not found in PATH) find matlab install parent dir
  !include "FileFunc.nsh" ; for Locate
  !include "WordFunc.nsh" ; for WordFind
  ${locate} "$PROGRAMFILES64"  "/L=D /M=MATLAB /G=0 /S=" "matlabInstallDirFoundCb" ; assuming user installed matlab in the default path
  StrCmp $R3 "" 0 init_loop
  ${Locate} "E:\" "/L=D /M=MATLAB /G=0 /S=" "matlabInstallDirFoundCb" ; if not default path, searching in E:\ (VM conf.)
  StrCmp $R3 "" 0 init_loop
  ${Locate} "F:\" "/L=D /M=MATLAB /G=0 /S=" "matlabInstallDirFoundCb" ; searching in F:\ just to test an non-existent disk on VM (TODO: delete)
  StrCmp $R3 "" fatal_error init_loop
  init_loop:
  ; loop into matlab directories (possible to have multiple versions)
  StrCpy $R0 1
  ; https://en.wikipedia.org/wiki/MATLAB#Release_history
  StrCpy $R1 "R2016a R2016b R2017a R2017b R2018a R2018b"
  loop:
      ${WordFind} $R1 " " +$R0 $R2
      IntOp $R0 $R0 + 1
      StrCmp $R0 8 continue ; it must be the last index + 2 !! e.g. if there are 6 versions to test, it must be 8
      ;MessageBox MB_OK "Version looked for: $R2 in $R3 (R0=$R0)"
      StrCpy $R4 ""
      ${Locate} $R3 "/L=D /M=$R2 /G=0 /S=" "matlabFoundCb"
      ;MessageBox MB_OK "R4=$R4"
      StrCmp $R4 "" loop +1
      MessageBox MB_OK 'Faust installed in $R4.' /SD IDOK
      ; MessageBox MB_YESNO 'Do you want to continue searching another version of Matlab to install Faust for ?' IDYES +2
      ;goto continue
      goto loop

  goto continue
  ;IfErrors 0 +2
  ;MessageBox MB_OK "Error" IDOK +2
  ;MessageBox MB_OK "$$R0=$R0"

  fatal_error:
      MessageBox MB_OK "Matlab installation path is not in default $PROGRAMFILES64. You will have to do the Faust setup on your own (you'll see how in the documentation)."  /SD IDOK

  continue:   

  ; install the doc
  SetOutPath $INSTDIR\doc\html
  File @PROJECT_BINARY_DIR@\doc\html\*

  ;ExecShell open "file:///$INSTDIR/doc/html/index.html" ; enable this link again when the bug issue #99 is fixed (doxypypy failing on windows)
  ExecShell open "https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/"
SectionEnd

Function matlabFoundCb
    StrCpy $R4 $R9

    ; as said above ExecWait won't work with matlab, wait manually
    ExecWait "matlab -nojvm -r $\"F=matfaust.rand(5,5);save(F, '$TEMP\nsisFFF.mat');exit$\""
    wait_faust_save:
	    ; do not edit matlabrc.m and startup.m if matfaust.rand is already usable (because it implies matlabrc.m and startup.m were already properly edited)
	    ifFileExists "$TEMP\nsisFFF.mat" matlabFoundCbEnd
	    sleep 1000
	    IntOp $R0 $R0 + 1
	    StrCmp $R0 15 0 wait_faust_save ; sleeps 15 secs max
    FileOpen $1 "$R4\toolbox\local\startup.m" a
    FileSeek $1 0 END ; do not erase start of file (but risk to add Faust path multiple times)
    ;FileWrite $1 "$\r$\naddpath(genpath('$INSTDIR\matlab'));$\r$\nmatfaust.enable_gpu_mod('silent', true)"
    FileWrite $1 "$\r$\nmatfaust_path = which('matfaust.version');if(isempty(matfaust_path));addpath('$INSTDIR\matlab');addpath('$INSTDIR\matlab\mex');addpath('$INSTDIR\matlab\mex\Release');addpath('$INSTDIR\matlab\tools');addpath('$INSTDIR\matlab\data');matfaust.enable_gpu_mod('silent', true);end"
    FileClose $1

    FileOpen $1 "$R4\toolbox\local\matlabrc.m" a
    FileSeek $1 0 END ; do not erase start of file (but risk to add Faust path multiple times)
    FileWrite $1 "$\r$\nif(exist('$INSTDIR'))oldpwd=pwd;startup_path=fullfile(matlabroot,'toolbox','local');cd(startup_path);startup_found=which('startup');if(isempty(startup_found))rehash('toolbox');startup_found=which('startup');if(isempty(startup_found)) disp('Faust startup code can''t be executed -- script not found. Please consult the online documentation to resolve this issue.');else;startup;end;else;startup;end;cd(oldpwd);end"
    FileClose $1

    # disable because it's redundant
;    FileOpen $1 "$DOCUMENTS\MATLAB\startup.m" a
;    IfErrors done
;    FileSeek $1 0 END
;    ;FileWrite $1 "$\r$\naddpath(genpath('$INSTDIR\matlab'));$\r$\nmatfaust.enable_gpu_mod('silent', true)"
;    FileWrite $1 "matfaust_path = which('matfaust.version');$\r$\nif(isempty(matfaust_path))"
;    FileWrite $1 "$\r$\naddpath('$INSTDIR\matlab');$\r$\naddpath('$INSTDIR\matlab\mex');$\r$\naddpath('$INSTDIR\matlab\tools');$\r$\naddpath('$INSTDIR\matlab\data')$\r$\nmatfaust.enable_gpu_mod('silent', true);"
;    FileWrite $1 "$\r$\nend"
;    FileClose $1
;    done:

    ;MessageBox MB_OK '$R0$\n$\nFaust bound into $R4.'
    ;MessageBox MB_YESNO 'Faust installed for $R4. Do you want to continue searching another version of Matlab to install Faust for ?' IDYES +2
    matlabFoundCbEnd:
    Delete "$TEMP\nsisFFF.mat"
    StrCpy $0 StopLocate

    Push $0
FunctionEnd

Function matlabInstallDirFoundCb
      StrCpy $R3 $R9
    ;MessageBox MB_OK 'Matlab install. dir: $R3'
    Push StopLocate
FunctionEnd


section "uninstall"
	; first, delete the uninstaller
	Delete "$INSTDIR\uninstall.exe"
	; do not RMDir /r $INSTDIR (danger: https://nsis.sourceforge.io/Reference/RMDir)
	RMDir /r "$INSTDIR\doc"
	RMDir /r "$INSTDIR\lib"
	RMDir /r "$INSTDIR\python"
	RMDir /r "$INSTDIR\matlab"
	; delete default install dir (doesn't work if user chose a custom dir)
	RMDir /r "$PROGRAMFILES64\Faust"
sectionEnd

