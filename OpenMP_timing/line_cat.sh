
































#!/bin/bash
set -Eeuo pipefail

#set -x



shopt -s extglob
files_to_merge=(./OpenMP_tst_"$1"x"$2"_+([0-9]).csv)

#exec 3>OpenMP_tst_"$1"x"$2".csv
#Open file descriptors
fds=()
for file in "${files_to_merge[@]}"; do
    exec {fd}< "$file"
    fds+=("$fd")
done

#Open output file descriptor
#out_fl=OpenMP_tst_"$1"x"$2".csv
exec {res_fl}> OpenMP_tst_"$1"x"$2".csv


IFS=
while true; do
     curr_line=();
     for fd in "${fds[@]}"; do
     
          if read -r line <&"$fd"; then
               curr_line+=("$line,")
          
          else
               break 2
          fi
          
     
     
     
     
     
     done
     printf "%s" "${curr_line[@]}" >&"$res_fl"
     printf "\n">&"$res_fl"




done


# Close file descriptors
for fd in "${fds[@]}"; do
    exec {fd}<&-
done

exec {res_fl}>&-


printf '\n%s\n' "${files_to_merge[@]}"
