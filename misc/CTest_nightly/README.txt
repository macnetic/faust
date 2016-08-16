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
# Run all day at 14h35 the nightly test and put results on CDASH
05 14  * * * /usr/bin/svn up /home/aleman/WORK/FAUST/faust_nightly/
35 14  * * * /home/aleman/WORK/FAUST/faust_nightly/run_nightly_XXX.sh

# run all day and send a mail 
40 11  * * * /Users/ci/CTest_nightly/run_nightly_OS-X.sh 2>&1| mail -s "Cron job execution" leman.adrien@gmail.com
##

4- 	You have the possibility to set the time of the precedent svn in the file 
CTestConfig.cmake


5- CTest OPTION: There are three types of dashboard submissions: 
    -Experimental means the current state of the project. An experimental submission can be performed at any time, usually interactively from the current working copy of a developer.

    -Nightly is similar to experimental, except that the source tree will be set to the state it was in at a specific nightly time. This ensures that all "nightly" submissions correspond to the state of the project at the same point in time. "Nightly" builds are usually done automatically at a preset time of day.

    -Continuous means that the source tree is updated to the latest revision, and a build / test cycle is performed only if any files were actually updated. Like "Nightly" builds, "Continuous" ones are usually done automatically and repeatedly in intervals.

To select the ctest mode, change the files faustTest.cmake to 
	CTEST_START("Experimental")
	CTEST_START("Nightly")
	CTEST_START("Continuous")


