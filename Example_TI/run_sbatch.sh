#!/bin/bash

temps=(300 650 1000)
for tmp_T in "${temps[@]}"
do
cp sbatch.cmd TI_${tmp_T}
cd TI_${tmp_T}
for file in ../configs/*
do
Name=$(basename $file)
cp sbatch.cmd tmp.cmd
sed -i "s/config/${Name}/g" tmp.cmd
sbatch tmp.cmd
done
cd -
done