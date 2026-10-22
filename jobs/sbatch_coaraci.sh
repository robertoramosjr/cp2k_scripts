#!/bin/bash
# sbatch with the coaraci whole-node defaults: par48-x, 3 days, 1 task x 48 threads
# (local CP2K 2026.2 ssmp), all node memory, CP2K_CLUSTER=coaraci exported (no guessing). Anything you pass comes after the
# defaults and overrides them (e.g. -t 1-00:00:00, or --ntasks-per-node=24 --cpus-per-task=2;
# tasks x cpus must stay 48, common_coaraci.sh refuses anything else).
#   ~/work_cp2k/cp2k_scripts/jobs/sbatch_coaraci.sh -J zro2-seq ~/work_cp2k/cp2k_scripts/jobs/job_sequential.sh LIST
exec sbatch --export=ALL,CP2K_CLUSTER=coaraci -p par48-x -t 3-00:00:00 -N 1 --ntasks-per-node=1 --cpus-per-task=48 --mem-per-cpu=2400M "$@"
