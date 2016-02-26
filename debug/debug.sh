#! /bin/bash

DEBUG_DIR='/home/tgautrai/faust2/debug'
BUILD_DIR_GPU='/home/tgautrai/faust2/test/src'
BUILD_DIR_CPU='/home/tgautrai/faust2/build'

LOG_GPU_FILE="$DEBUG_DIR/log_compil_gpu"
LOG_CPU_FILE="$DEBUG_DIR/log_compil_cpu"
LOG_GPU_OUTPUT="$DEBUG_DIR/log_output_gpu"
LOG_CPU_OUTPUT="$DEBUG_DIR/log_output_cpu"


cd $DEBUG_DIR

rm -f *.tmp

cd $BUILD_DIR_GPU
make cleanall 
date '+---------------------------------- %Y-%m-%d %H:%M:%S ----------------------------------' >> $LOG_GPU_FILE
make -j8  >>$LOG_GPU_FILE 2>&1
cd $DEBUG_DIR
date '+---------------------------------- %Y-%m-%d %H:%M:%S ----------------------------------' >> $LOG_GPU_OUTPUT
$BUILD_DIR_GPU/hierarchical_fact_test_cu.out | tee -a $LOG_GPU_OUTPUT

cd $BUILD_DIR_CPU 
make clean
date '+---------------------------------- %Y-%m-%d %H:%M:%S ----------------------------------' >> $LOG_CPU_FILE
make -j8  >>$LOG_CPU_FILE 2>&1
cd $DEBUG_DIR
date '+---------------------------------- %Y-%m-%d %H:%M:%S ----------------------------------' >> $LOG_CPU_OUTPUT
$BUILD_DIR_CPU/testing/bin/faust_hier | tee -a $LOG_CPU_OUTPUT
