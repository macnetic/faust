#! /bin/bash

DEBUG_DIR='/home/tgautrai/faust2/debug'
BUILD_DIR_GPU='/home/tgautrai/faust2/test/src'
BUILD_DIR_CPU='/home/tgautrai/faust2/build'

LOG_DIR="$DEBUG_DIR/log"
DATE_DEBUT_SCRIPT=`date '+%Y-%m-%d_%H-%M-%S'`
LOG_GPU_FILE="$DEBUG_DIR/log/compilation_${DATE_DEBUT_SCRIPT}_GPU.log"
LOG_CPU_FILE="$DEBUG_DIR/log/compilation_${DATE_DEBUT_SCRIPT}_CPU.log"
LOG_GPU_OUTPUT="$DEBUG_DIR/log/execution_${DATE_DEBUT_SCRIPT}_GPU.log"
LOG_CPU_OUTPUT="$DEBUG_DIR/log/execution_${DATE_DEBUT_SCRIPT}_CPU.log"

NOM_FICHIER="$VAR/`date '+%Y-%m-%d_%H-%M-%S_CPU.log'`"

export COMPILE_SPMAT=1
export COMPILE_GPU=1
export COMPILE_TIMERS=1
export OPENBLAS_NUM_THREADS=`cat /proc/cpuinfo |grep processor |wc -l`

cd $DEBUG_DIR

rm -f *.tmp

cd $BUILD_DIR_GPU
make cleanall 
make -j8  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_GPU.log" 2>&1
cd $DEBUG_DIR
#$BUILD_DIR_GPU/hierarchical_fact_test_cu.out | tee -a "$LOG_DIR/faust_hier_`date '+%Y-%m-%d_%H-%M-%S'`_GPU.log"
$BUILD_DIR_GPU/MEG_fact_cu.out | tee -a "$LOG_DIR/meg_`date '+%Y-%m-%d_%H-%M-%S'`_GPU.log"

cd $BUILD_DIR_CPU 
make clean
make -j8  >> "$LOG_DIR/compilation_`date '+%Y-%m-%d_%H-%M-%S'`_CPU.log" 2>&1
cd $DEBUG_DIR
#$BUILD_DIR_CPU/testing/bin/faust_hier | tee -a "$LOG_DIR/faust_hier_`date '+%Y-%m-%d_%H-%M-%S'`_CPU.log"
$BUILD_DIR_CPU/testing/bin/meg | tee -a "$LOG_DIR/meg_`date '+%Y-%m-%d_%H-%M-%S'`_CPU.log"
