#! /bin/bash


LOG_DIR="$DEBUG_DIR/log"

export COMPILE_SPMAT=1
export COMPILE_GPU=1
#export COMPILE_TIMERS=1
export OPENBLAS_NUM_THREADS=`cat /proc/cpuinfo |grep processor |wc -l`

cd $DEBUG_DIR

rm -f *.tmp

declare -a real_type=("float" "double")


for precision in "${real_type[@]}"
do
	sed -i "s/typedef[[:space:]][[:alpha:]]\+[[:space:]]\([[:alpha:]]\+\)/typedef $precision \1/g"  $BUILD_DIR_GPU/hierarchical_fact_test_cu.cpp $BUILD_DIR_GPU/MEG_fact_cu.cpp $BUILD_DIR_CPU/testing/src/hierarchical_fact_test.cpp $BUILD_DIR_CPU/testing/src/MEG_fact.cpp 

	cd $BUILD_DIR_GPU
	make cleanall 
	make  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log" 2>&1
	cd $DEBUG_DIR
	$BUILD_DIR_GPU/hierarchical_fact_test_cu.out | tee -a "$LOG_DIR/faust_hier_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"
	#find $FAUST_DATA_TEST -name "config_*rows*cols*.mat" -exec $BUILD_DIR_GPU/MEG_fact_cu.out {} \;  | tee -a "$LOG_DIR/meg_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"
	$BUILD_DIR_GPU/MEG_fact_cu.out | tee -a "$LOG_DIR/meg_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"

	cd $BUILD_DIR_CPU 
	make clean
	make  meg faust_hier >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log" 2>&1
	cd $DEBUG_DIR
	$BUILD_DIR_CPU/testing/bin/faust_hier | tee -a "$LOG_DIR/faust_hier_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
	#find $FAUST_DATA_TEST -name "config_*rows*cols*.mat" -exec $BUILD_DIR_CPU/testing/bin/meg {} \;   | tee -a "$LOG_DIR/meg_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
	$BUILD_DIR_CPU/testing/bin/meg | tee -a "$LOG_DIR/meg_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
done
