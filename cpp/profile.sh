#!/bin/bash
#SBATCH --exclusive

function filter_average_median() {
  "$@" | grep -i -E "average|median"
}

function for_core_count() {
  w=$1
  shift 1
  echo "Workers: $w"
  taskset -c 0-$((w - 1)) "$@"
}

function for_all_cores() {
  max_w=$1
  shift 1
  for ((w = 1; w <= $max_w; w++)); do
    for_core_count $w "$@"
  done
}

function for_all_grainsizes() {
  max_grain=$1
  shift 1
  for ((g = 1; g <= $max_grain; g *= 2)); do
    echo "Grainsize: $g"
    CILK_GRAINSIZE=$g "$@"
  done
}

function for_all_power2_params() {
  start=$1
  end=$2
  shift 2
  for ((p = $start; p <= $end; p *= 2)); do
    echo "Param: $p"
    "$@ $p"
  done
}


function for_all_params() {
  start=$1
  end=$2
  shift 2
  for ((p = $start; p <= $end; p++)); do
    echo "Param: $p"
    "$@ $p"
  done
}

echo "command: '$@'"
"$@"
