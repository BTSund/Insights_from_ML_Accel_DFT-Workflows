#!/bin/bash

# Define the source directory
CONFIG_DIR="../../configs"

temps=(300 650 1000)
for tmp_T in "${temps[@]}"
do
mkdir TI_${tmp_T}
cp -r starterpack TI_${tmp_T}
sed -i "s/set_temp/${tmp_T}/g" TI_${tmp_T}/starterpack/*
cd TI_${tmp_T}

for rep in {0..4}
do
mkdir "rep_${rep}"
cd "rep_${rep}"
for File in "$CONFIG_DIR"/*
do
    # Extract phase name (e.g., Diamond, SH)
    Name=$(basename "$File")
    echo "Processing Phase: $Name"

    # --- PHASE-SPECIFIC LOGIC ---
    # Define which pressures to run based on the folder name
    case "$Name" in
        "CON_Dia")
            Press=(0 10 15)
            ;;
        "CON_B-Sn"|"CON_IMMA")
            Press=(10 12 14 16 18 20 22)
            ;;
        "CON_sh")
            Press=(18 20 22 30 34 38 40)
            ;;
        "*") 
            Press=(30 34 38 40 50 54 58)
            ;;
        "CON_opt_HCP") 
            Press=(83 88)
            ;;
        "CON_dhcp") 
            Press=(62 66 70 74 78 82 86 90 95 100)
            ;;
        "CON_FCC") 
            Press=( 70 74 78 82 86 90 95 100)
            ;;
        *)
            # Default pressures if the phase name isn't specifically listed
            echo $Name
            ;;
    esac
    # ----------------------------

    mkdir -p "$Name"
    cp "$File" "$Name/data.lammps"
    cd "$Name"

    for tmp_P in "${Press[@]}"; do
        # Convert GPa to atm (1 GPa ≈ 9869.23 atm)
        P_atm=$(( tmp_P * 9869 ))
        
        # Zero-fill pressure for clean directory naming (e.g., Press_010)
        Name2=$(printf "Press_%03d" "$tmp_P")
        
        mkdir -p "$Name2"
        cp ../../starterpack/* "$Name2/"
        cp data.lammps "$Name2/"
        
        # Update the LAMMPS input file with the calculated pressure
        sed -i "s/set_press/$P_atm/g" "$Name2/lmp.in"
        sed -i "s/rand_1/10${rep}/g" "$Name2/lmp.in"
        sed -i "s/rand_2/1${rep}/g" "$Name2/TI.in"
    done
    
    cd ..
    done
    cd ..
done
cd .. 
done