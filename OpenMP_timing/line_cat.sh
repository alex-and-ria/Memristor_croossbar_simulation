
































#!/bin/bash
set -Eeuo pipefail

set -x



shopt -s extglob
files_to_merge=(./OpenMP_tst_"$1"x"$2"_+([0-9]).csv)

exec 3>OpenMP_tst_"$1"x"$2".csv


while true; do
     IFS=,
     for file in "${files_to_merge[@]}"; do
          if read -r line < "$file"; then
               printf "%s," "$line" >&3
          
          else
               break 2
          fi
          
          
          
     
     done
     printf "\n" >&3





done





exec 3>$-
printf '\n%s\n' "${files_to_merge[@]}"
