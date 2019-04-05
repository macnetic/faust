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

  ; install python wrapper
  SetOutPath $INSTDIR\python\pyfaust
  File @PROJECT_BINARY_DIR@\wrapper\python\*py
  SetOutPath $INSTDIR\python
  File @PROJECT_BINARY_DIR@\wrapper\python\*pyd
  File @PROJECT_BINARY_DIR@\wrapper\python\*pxd

  ; install matlab wrapper
  SetOutPath $INSTDIR\matlab
  File /nonfatal /r /x old_matlab @PROJECT_BINARY_DIR@\wrapper\matlab\*.m @PROJECT_BINARY_DIR@\wrapper\matlab\*.@MEX_EXT@ @PROJECT_BINARY_DIR@\wrapper\matlab\*.mat
  ; nonfatal useful in case of data *.mat not used/present (because they are downloaded at installation)



  ; post install pyfaust auto-setup in environment (only works if python is installed in path) 
  ${StrRep} '$0' $TEMP '\' '\\'
  Exec "python -c $\"import site;dir=site.getsitepackages()[1];f=open('$0\\tmp_site_pkg', 'w');f.write(dir);f.close()$\""
  IfErrors 0 +2
  MessageBox MB_OK "Error: no python found into your PATH environment variable. You'll have to do the Faust setup manually (you'll see how in the documentation)." IDOK +2
  MessageBox MB_OK "Faust installed in your python environment (the version found into your PATH environment variable)."

  FileOpen $1 "$TEMP\tmp_site_pkg" r
  FileRead $1 $2
  FileClose $1
  ;MessageBox MB_OK "$2"

  SetOutPath $2\pyfaust 
  File /r @PROJECT_BINARY_DIR@\wrapper\python\pyfaust\*py
  SetOutPath $2
  File @PROJECT_BINARY_DIR@\wrapper\python\*pyd
  File @PROJECT_BINARY_DIR@\wrapper\python\*pxd

  ;add data path in __init__.py
  FileOpen $1 "$2\pyfaust\__init__.py" a
  ; go at the end of file (but with append mode is that necessary ?)
  FileSeek $1 0 END
  ; add the install path
  FileWrite $1 "$\r$\n_NSI_INSTALL_PATH='$INSTDIR'"
  FileClose $1

  ; download data into matlab wrapper data folder
  ; create data folder
  ; try with python first
  ${StrRep} '$3' $INSTDIR '\' '\\'
  Exec "python -c $\"from os import mkdir; mkdir('$3\\matlab\\data')$\""
  ; download data
  Exec "python $2\pyfaust\datadl.py $\"$INSTDIR\matlab\data$\""
  IfErrors 0 after_data_dl
  MessageBox MB_OK "Downloading FAuST data with python has failed, now trying with powershell."
  ClearErrors
  ExecWait "powershell Invoke-WebRequest $\"@REMOTE_DATA_URL@/@REMOTE_DATA_FILE@$\" -O $\"$TEMP\@REMOTE_DATA_FILE@$\"" ; ExecWait because unzipping needs download finished
  Exec "powershell Expand-Archive -Force $\"$TEMP\@REMOTE_DATA_FILE@$\" '$INSTDIR\matlab\data'" ; output folder data is created auto. ; simple quote used to avoid powershell to think there is two arguments when we meant one argument for the dest. path (double quote doesn't allow that).
  IfErrors 0 after_data_dl
  MessageBox MB_OK "Error downloading FAuST data. You'll need to download manually (please check the documentation)." IDOK after_data_dl

  after_data_dl:

  ; post install matfaust auto-setup
  !include "FileFunc.nsh" ; for Locate
  !include "WordFunc.nsh" ; for WordFind
  ; find matlab install parent dir ; TODO: change E:\ to default $PROGRAMFILES64
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
  StrCpy $R1 "R2016a R2016b R2017a R2017b R2018a"
  loop:
	  ${WordFind} $R1 " " +$R0 $R2
	  IntOp $R0 $R0 + 1
	  ; MessageBox MB_OK "Version looked for: $R2"
	  StrCmp $R0 6 continue
	  StrCpy $R4 ""
	  ${Locate} $R3 "/L=D /M=$R2 /G=0 /S=" "matlabFoundCb"
	  StrCmp $R4 "" loop +1
	  MessageBox MB_OK 'Faust installed for $R4.'
	  ; MessageBox MB_YESNO 'Do you want to continue searching another version of Matlab to install Faust for ?' IDYES +2
	  ;goto continue
	  goto loop

  goto continue
  ;IfErrors 0 +2
  ;MessageBox MB_OK "Error" IDOK +2
  ;MessageBox MB_OK "$$R0=$R0"

  fatal_error:
  	MessageBox MB_OK "Matlab installation path is not in default $PROGRAMFILES64. You will have to do the Faust setup on your own (you'll see how in the documentation)."

  continue:	

  ; install the doc
  SetOutPath $INSTDIR\doc\html
  File @PROJECT_BINARY_DIR@\doc\html\*

  ExecShell open "file:///$INSTDIR/doc/html/index.html"

SectionEnd

Function matlabFoundCb
	StrCpy $R4 $R9

	FileOpen $1 "$R4\toolbox\local\startup.m" a
	FileSeek $1 0 END ; do not erase start of file (but risk to add Faust path multiple times)
	FileWrite $1 "$\r$\naddpath(genpath('$INSTDIR\matlab'))"
	FileClose $1

	FileOpen $1 "$DOCUMENTS\MATLAB\startup.m" w
	IfErrors done
	FileSeek $1 0 END
	FileWrite $1 "$\r$\naddpath(genpath('$INSTDIR\matlab'))"
	FileClose $1
	done:

	;MessageBox MB_OK '$R0$\n$\nFaust bound into $R4.' 
	;MessageBox MB_YESNO 'Faust installed for $R4. Do you want to continue searching another version of Matlab to install Faust for ?' IDYES +2
	StrCpy $0 StopLocate

	Push $0
FunctionEnd

Function matlabInstallDirFoundCb
  	StrCpy $R3 $R9
	;MessageBox MB_OK 'Matlab install. dir: $R3'
	Push StopLocate
FunctionEnd
