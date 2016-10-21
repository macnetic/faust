@ECHO OFF

:: par defaut le test est bon
set TESTPASS=0 
set middle=%1
set "debut= delete MATLAB_TEST.txt;diary MATLAB_TEST.txt;try;testpass=0;"
set "fin= catch ME ;testpass=2;disp(getReport(ME)); end ; disp (testpass); diary OFF;exit(testpass);"
set test=%debut%%middle%%fin%
::echo DEBUT : %debut%
::echo MIDLLE : %middle%
::echo FIN : %fin%
::echo TEST : %test%

::transform the "@" character into "," 
set test=%test:@=,%

::transform the "$" character into ";"
set test=%test:$=;%


echo MATLAB COMMAND : %test%

:: launch the test in matlab an wait for matlab to finish
matlab -nojvm -wait -r """%test%"""

:: si le code renvoyer par la derniere cmd est sup a 2 on met TESPASS a 2 et erreur
IF ERRORLEVEL 2 set TESTPASS=2




:: display the content of the MATLAB Command Window stored in a txtfile
echo ****** MATLAB WINDOW ******* 
type MATLAB_TEST.txt

ECHO ERROR_LEVEL %TESTPASS%
::renvoie le code erreur
exit %TESTPASS%
