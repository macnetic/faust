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
	make -j$OPENBLAS_NUM_THREADS  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log" 2>&1
	cd $DEBUG_DIR
#$BUILD_DIR_GPU/hierarchical_fact_test_cu.out | tee -a "$LOG_DIR/faust_hier_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"
	$BUILD_DIR_GPU/MEG_fact_cu.out /home/tgautrai/faust2/test/data/config_4rows_4cols_4facts.mat /home/tgautrai/faust2/test/data/config_4rows_4cols_10facts.mat /home/tgautrai/faust2/test/data/config_8rows_8cols_4facts.mat /home/tgautrai/faust2/test/data/config_8rows_8cols_10facts.mat /home/tgautrai/faust2/test/data/config_16rows_16cols_4facts.mat /home/tgautrai/faust2/test/data/config_16rows_16cols_10facts.mat /home/tgautrai/faust2/test/data/config_32rows_32cols_4facts.mat /home/tgautrai/faust2/test/data/config_32rows_32cols_10facts.mat /home/tgautrai/faust2/test/data/config_64rows_64cols_4facts.mat /home/tgautrai/faust2/test/data/config_64rows_64cols_10facts.mat /home/tgautrai/faust2/test/data/config_128rows_128cols_4facts.mat /home/tgautrai/faust2/test/data/config_128rows_128cols_10facts.mat /home/tgautrai/faust2/test/data/config_256rows_256cols_4facts.mat /home/tgautrai/faust2/test/data/config_256rows_256cols_10facts.mat /home/tgautrai/faust2/test/data/config_512rows_512cols_4facts.mat /home/tgautrai/faust2/test/data/config_512rows_512cols_10facts.mat /home/tgautrai/faust2/test/data/config_1024rows_1024cols_4facts.mat /home/tgautrai/faust2/test/data/config_1024rows_1024cols_10facts.mat  | tee -a "$LOG_DIR/meg_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"

	cd $BUILD_DIR_CPU 
	make clean
	make -j$OPENBLAS_NUM_THREADS  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log" 2>&1
	cd $DEBUG_DIR
#$BUILD_DIR_CPU/testing/bin/faust_hier | tee -a "$LOG_DIR/faust_hier_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
	$BUILD_DIR_CPU/testing/bin/meg /home/tgautrai/faust2/test/data/config_4rows_4cols_4facts.mat /home/tgautrai/faust2/test/data/config_4rows_4cols_10facts.mat /home/tgautrai/faust2/test/data/config_8rows_8cols_4facts.mat /home/tgautrai/faust2/test/data/config_8rows_8cols_10facts.mat /home/tgautrai/faust2/test/data/config_16rows_16cols_4facts.mat /home/tgautrai/faust2/test/data/config_16rows_16cols_10facts.mat /home/tgautrai/faust2/test/data/config_32rows_32cols_4facts.mat /home/tgautrai/faust2/test/data/config_32rows_32cols_10facts.mat /home/tgautrai/faust2/test/data/config_64rows_64cols_4facts.mat /home/tgautrai/faust2/test/data/config_64rows_64cols_10facts.mat /home/tgautrai/faust2/test/data/config_128rows_128cols_4facts.mat /home/tgautrai/faust2/test/data/config_128rows_128cols_10facts.mat /home/tgautrai/faust2/test/data/config_256rows_256cols_4facts.mat /home/tgautrai/faust2/test/data/config_256rows_256cols_10facts.mat /home/tgautrai/faust2/test/data/config_512rows_512cols_4facts.mat /home/tgautrai/faust2/test/data/config_512rows_512cols_10facts.mat /home/tgautrai/faust2/test/data/config_1024rows_1024cols_4facts.mat /home/tgautrai/faust2/test/data/config_1024rows_1024cols_10facts.mat | tee -a "$LOG_DIR/meg_${precision}_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
done
