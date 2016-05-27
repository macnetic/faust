##################################################################
############## README : run nightly faust project ################
##################################################################

This tool aims to launch automatically a shell-bash SCRIPT : 
 
1- 	First at all, please chek out only the files presents in 
	following directory : devcpp/test/CTest_nightly 
	svn checkout --username XXX https://scm.gforge.inria.fr/authscm/XXX/svn/faust/trunk/devcpp/misc/CTest_nightly/ ./

2-	Configure your personnal run shell "run_nightly_PLATFORM.sh" with corresponding export of library and PATH used.

3- 	configure the "crontab" tool for UNIX or shulder for windows
	crontab -e : 
	example :
# Run all day at 21h01 the nightly test and put results on CDASH
05 14  * * * /usr/bin/svn up /home/aleman/WORK/FAUST/faust_nightly/
35 14  * * * /home/aleman/WORK/FAUST/faust_nightly/run_nightly_XXX.sh
##

4- 	You have the possibility to set the time of the precedent svn in the file 
CTestConfig.cmake


