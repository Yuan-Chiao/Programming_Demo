#!/bin/bash
#SBATCH --job-name=${job_key}
#SBATCH -p gpuq
#SBATCH --gres=gpu:2
#SBATCH --output /home/yuan/Desktop/DeepLearning/Placenta_seg/run.log

work_dir="/home/yuan/Desktop/DeepLearning/Placenta_seg/"
job_key="Pla_Seg"

module load tensorflow python
echo "----------------------------------------------------"
date
echo $work_dir
cd $work_dir
echo NAME: $job_key
echo "---------------------start--------------------------"
python ${work_dir}/Main_dp.py
echo "----------------------done--------------------------"
date
