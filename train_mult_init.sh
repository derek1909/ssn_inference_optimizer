#!/bin/bash

N_INSTANCES=10
MAX_ITER_1=250
T_MAX=0.5
LAMBDA_VARIANCE=2.0
N_TARGETS=5
TGTS_DIR='targets'
MAIN_DIR='results'
for ((i=1; i<=N_INSTANCES; i++))
do
    mkdir $MAIN_DIR/instance_$i
    cp $TGTS_DIR/* $MAIN_DIR/instance_$i
    echo -e "\nTraining iter: $i"
    ./train.native -d $MAIN_DIR/instance_$i -m 50 -n_targ $N_TARGETS -lambda_noise 0.0 -t_max $T_MAX -max_iter $MAX_ITER_1 -lambda_var $LAMBDA_VARIANCE -rand_init true -n_trials 50 -save_every 50 -eta 0.005 -tau_eta 0.02
    ./train.native -d $MAIN_DIR/instance_$i -m 50 -n_targ $N_TARGETS -lambda_noise 0.0 -t_max $T_MAX -max_iter $MAX_ITER_1 -lambda_var $LAMBDA_VARIANCE -rand_init false -save_every 50 -reuse $MAIN_DIR/instance_$i/state_vector.bin -no_progression -tau_eta 0.02 -lambda_slow 0.1
done

