##################################################################
############## README : Continuous Integration 		##############
##############          of the FAuST project 		##############
##################################################################


##################################################################
############## CDASH 	 							##############
##################################################################
The Continuous Integration for the project FAUST is based on the CDASH tool (see. http://www.cdash.org/)
The building and test are available on the public link: 
http://cdash.inria.fr/CDash/index.php?project=faust

CDASH is configured in the Cmakefile.txt files. 

 
CDASH is deployed on the CI inria platform. For that you must have an login on ci inria. 
Then go to the  https://ci.inria.fr/project/faust/show to manage the CI: 
In slave, you can add/modify the virtual Machine (linux mac, windows). 
If you want to get more actions / information on slaves (add a disk, create a template, etc.), you can access CloudStack using the same credentials as the CI portal. The Domain must be ci/faust. CloudStack 
https://sesi-cloud-ctl1.inria.fr/client/


https://ci.inria.fr/

connect to CI of inria: 
https://ci.inria.fr/project/faust/ to manage your virtual machine. 



This tool aims to launch automatically a shell-bash SCRIPT on a machine (local, or virtual ...) : 
 
1- 	First at all, please chek out only the files presents in 
	following directory : devcpp/misc/CTest_nightly 
	svn checkout --username XXX https://scm.gforge.inria.fr/authscm/XXX/svn/faust/trunk/devcpp/misc/CTest_nightly/ ./

2-	Create and Configure your personnal run script (".sh" for Unix and ".bat" for Windows) "run_nightly_PLATFORM.sh" with corresponding configuration. (see existing run_nighlty_X.sh files for more informations)

3- 	configure the "crontab" tool for UNIX or "shulder" tool for windows
	crontab -e : 
	example :
# Run all day at 14h35 the nightly test and put results on CDASH
05 14  * * * /usr/bin/svn up /home/aleman/WORK/FAUST/faust_nightly/
35 14  * * * /home/aleman/WORK/FAUST/faust_nightly/run_nightly_XXX.sh

# run all day and send a mail 
40 11  * * * /Users/ci/CTest_nightly/run_nightly_OS-X.sh 2>&1| mail -s "Cron job execution" youremail@inria.fr
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




##################################################################
##################### Acces to the virtual machine
##################################################################
cf. https://ci.inria.fr/project/faust/slaves  --> onget connect

Windows : 
- Windows MinGW Openblas:
$ ssh <login>@ci-ssh.inria.fr -L 3380:Windows-MinGW:3389

- VM windowsTWO 7: MinGW without Openblas
Password of new VM windowsTWO is  bF8brkgqs
1 - Create an SSH tunnel
$ ssh <login>@ci-ssh.inria.fr -L 3380:windowsTWO:3389 
2- Open graphical interface : 
$ rdesktop -k fr -g 90% -u ci -p ci 127.0.0.1:3380 

- Password of new VM faust-windows7VisualStudio is  yX4arrrub
Password has been reset to xH8cteihd
1- 
$ ssh <login>@ci-ssh.inria.fr -L 3380:faust-windows7VisualStudio:3389
2- 
$ rdesktop -k fr -g 90% -u ci -p ci 127.0.0.1:3380 

MAC: 
1- 
	$ ssh ci@mac.ci
2- 
	Password of new VM faust-mac-copy is  mD4pszrtb
	$ ssh ci@faust-mac-copy.ci

3- 
	$ ssh ci@MAC-Xcode.ci


LINUX: 
pour machine virtuelle unbuntu 
$ ssh -X ci@faust-ubuntu14.ci
ci@faust-ubuntu14.ci's password: ci

$ ssh ci@faust-ubuntu-Copy.ci
password: ci


######################################################################################""
Pour ajouter un data disk sur les machine virtuelles :
							1- https://sesi-cloud-ctl1.inria.fr pour ajouter un volume à l'instance
							2- stop and restart
							3- mount the volume
								
In case of UNIX platform : 
	lsblk --> list les disk et les partitions correspondantes
	admin@mariadb:~$ sudo mkfs.ext4 /dev/vda
	Allocating group tables: done
	Writing inode tables: done
	Creating journal (32768 blocks): done
	Writing superblocks and filesystem accounting information: done

	admin@mariadb:~$ sudo mkdir /mnt/database

	admin@mariadb:~$ sudo mount /dev/vda /mnt/database/

	admin@mariadb:~$ df -h
	Filesystem      Size  Used Avail Use% Mounted on
	/dev/vda1       158G  1.1G  150G   1% /
	/dev/vdb        9.7T   38M  9.3T   1% /mnt/database

In case of Windows platform : disk managment option.
	https://wiki.inria.fr/ciportal/FAQ#
	I_added_an_extra_hard_drive_to_the_VM_but_it_is_not_detected_by_the_OS
	il faut installer les drivers scsi via I added an extra hard drive to the VM but it is not detected by the OS
	You may have a problem with the SCSI driver.
    go to the device manager (Right click on My Computer -> Manage -> Device Manager)
    select the SCSI controller
    install the latest viostor driver from Redhad (https://fedoraproject.org/wiki/Windows_Virtio_Drivers#Direct_download). Be careful to select the right driver (windows version & architecture) or windows will crash very badly (do not even use automatic detection, select the driver manually from the right folder).
##################################################################











