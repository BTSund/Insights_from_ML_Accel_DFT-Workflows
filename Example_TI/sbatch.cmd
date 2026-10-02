#!/bin/bash
#SBATCH -J TI
#SBATCH -N 5
#SBATCH --ntasks-per-node 48
#SBATCH -t 24:00:00
#SBATCH -p [queue]
#SBATCH -o stdoutmsg

module load cmake intel impi

# Capture the base directory to avoid navigation errors
BASE_DIR=$(pwd)

for n in {0..4}
do
    # Corrected variable name to 'n'
    OFFSET=$(( n * 48 ))
    Dir="rep_${n}"

    # Run each replica in the background
    (
        if [ ! -d "$Dir" ]; then echo "Directory $Dir not found"; exit; fi
        cd "$Dir" || exit

        for Folder in config/Press*/
        do
            cd "$BASE_DIR/$Dir/$Folder" || continue

            # Check if job is already done
            if [[ -f backward_fl.dat ]]; then
                last_line=$(tail -n 1 backward_fl.dat)
                if [[ "$last_line" == *"30000"* ]]; then
                    echo "Skipping $Folder (completed)"
                    continue
                fi
            fi

            # --- Calculation 1 ---
            ibrun -n 48 -o $OFFSET /work2/09981/btsund/stampede3/codes/New_Lammps-test/etc/lmp/exe/lmp_mpi_chimes -i lmp.in > out.NVT.lammps

            # Extract averages
            lx=$(awk 'NR==1 { print $2 }' averages.txt)
            ly=$(awk 'NR==1 { print $3 }' averages.txt)
            lz=$(awk 'NR==1 { print $4 }' averages.txt)
            xy=$(awk 'NR==1 { print $5 }' averages.txt)
            xz=$(awk 'NR==1 { print $6 }' averages.txt)
            yz=$(awk 'NR==1 { print $7 }' averages.txt)

            # Update TI.in
            sed -i "s/set_lx/${lx}/g" TI.in
            sed -i "s/set_ly/${ly}/g" TI.in
            sed -i "s/set_lz/${lz}/g" TI.in
            sed -i "s/set_xy/${xy}/g" TI.in
            sed -i "s/set_xz/${xz}/g" TI.in
            sed -i "s/set_yz/${yz}/g" TI.in

            # --- Calculation 2 ---
            ibrun -n 48 -o $OFFSET /work2/09981/btsund/stampede3/codes/New_Lammps-test/etc/lmp/exe/lmp_mpi_chimes -i TI.in > out.TI.lammps
        done
    ) &
done

# CRITICAL: Wait for all background subshells to finish
wait
echo "All tasks complete."