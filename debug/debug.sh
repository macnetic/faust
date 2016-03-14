#! /bin/bash


LOG_DIR="$DEBUG_DIR/log"

export COMPILE_SPMAT=1
export COMPILE_GPU=1
export COMPILE_TIMERS=1
export OPENBLAS_NUM_THREADS=`cat /proc/cpuinfo |grep processor |wc -l`

cd $DEBUG_DIR

rm -f *.tmp

cd $BUILD_DIR_GPU
make cleanall 
make -j$OPENBLAS_NUM_THREADS  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log" 2>&1
cd $DEBUG_DIR
#$BUILD_DIR_GPU/hierarchical_fact_test_cu.out | tee -a "$LOG_DIR/faust_hier_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"
$BUILD_DIR_GPU/MEG_fact_cu.out | tee -a "$LOG_DIR/meg_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_GPU.log"

cd $BUILD_DIR_CPU 
make clean
make -j$OPENBLAS_NUM_THREADS  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log" 2>&1
cd $DEBUG_DIR
#$BUILD_DIR_CPU/testing/bin/faust_hier | tee -a "$LOG_DIR/faust_hier_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
$BUILD_DIR_CPU/testing/bin/meg | tee -a "$LOG_DIR/meg_`date '+%Y-%m-%d_%H-%M-%S'`_`hostname`_`whoami`_CPU.log"
