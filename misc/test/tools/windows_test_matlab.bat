@ECHO OFF

:: par defaut le test est bon
set ERRORLEVEL=0
 
set matlab_command_core=%1
set outputfile=%2%.txt
::transform the "/" character into "\" to have a valid path
set outputfile=%outputfile:/=\%

echo "OUTPUTFILE %outputfile%"

set "matlab_command_begin= delete %outputfile%;diary %outputfile%;try;testpass=0;"
set "matlab_command_end= catch ME ;testpass=2;disp(getReport(ME));end; diary OFF;exit(testpass);"
set matlab_command=%matlab_command_begin%%matlab_command_core%%matlab_command_end%
::echo DEBUT : %debut%
::echo MIDLLE : %middle%
::echo FIN : %fin%
::echo TEST : %test%

::transform the "@" character into "," 
set matlab_command=%matlab_command:@=,%

::transform the "$" character into ";"
set matlab_command=%matlab_command:$=;%


echo MATLAB COMMAND : %matlab_command%

:: launch the test in matlab an wait for matlab to finish
matlab -nodesktop -wait -r """%matlab_command%"""

:: si le code renvoyer par la derniere cmd est sup a 2 => erreur
IF ERRORLEVEL 2 set ERRORLEVEL=2




:: display the content of the MATLAB Command Window stored in a txtfile
echo ****** MATLAB WINDOW ******* 
type %outputfile%
IF ERRORLEVEL 1 set ERRORLEVEL=1

ECHO ERROR_LEVEL %ERRORLEVEL%
::renvoie le code erreur
exit %ERRORLEVEL%
