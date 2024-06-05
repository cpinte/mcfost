#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4g
#SBATCH --time=160:00:00
#SBATCH --job-name=mcfost
#SBATCH --mail-type=ALL


# results to be uploaded to ipag-nfs.u-ga.fr:webpage/mcfost/
#export mcfost=$MCFOST_INSTALL/../src/mcfost
export OMP_NUM_THREADS=1 # 1 cpu to compute reference
export mcfost=$(pwd)/../src/mcfost
export MCFOST_UTILS=$(pwd)/../utils

for dir in test_data/*; do
    param=`basename "$dir".para`
    pushd "$dir"

    echo "---------------------------------------------"
    pwd
    echo "---------------------------------------------"

    rm -rf data_* *.tmp
    if [ "$param" == "discF_00500.para" ]; then
        opt="-phantom discF_00500 -not_random_Voronoi"
    else
        opt=""
    fi
    $mcfost $param $opt -mol
    $mcfost $param $opt -img 1.0
    $mcfost $param $opt -img 10
    $mcfost $param $opt -img 100
    $mcfost $param $opt -img 1000
    rm -f *.tmp
    popd
done
