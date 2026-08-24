
































#!/bin/bash
set -Eeuo pipefail

set -x




for ((batch_size=1;batch_size<=1;batch_size++));do
for ((m=256;m<=256;m*=2)); do
     for((n=m;n==m;n*=2)); do
          exec {out_file}> "${m}x${n}_${batch_size}".csv
          printf "\ndecomp, L_sol, U_sol,type,online,offline,nds_n_req,sol,sol_full,sol_agg\n" >&${out_file}
          for((nds_n_req=1;nds_n_req<=2*m*n-1;nds_n_req+=100)); do
               rm -f tmp.csv
               for ((q=0;q<5;q++)); do
                    /share/reconfig/matlab/bin/matlab -batch "m=$m; n=$n; batch_size=$batch_size; nds_n_req=$nds_n_req; speed_tst" >>tmp.csv
               
               done
               sed -z 's/ans = \n\n    "\n     //g; s/"\n\n,/,/g; s/\nans = \n\n    ",/,/g; s/"\n\n//gi; s/\[[^]]*\]..\n//gi' tmp.csv | awk -v nds_n_req=$nds_n_req '
               BEGIN{
                    FS=","
                    OFS=","
               }
               {
                    for (i = 1; i <= NF; i++){
                         sums[i-1]+=$i
                         counts[i-1]++
                    }

               }
               END{
                    for(i=0;i<NF;i++){
                         printf "%s%g", (i==0)?"":"," ,sums[i]/counts[i]
                    }
                    printf("\n")
                    
               }
               ' >&${out_file}
          done
          exec {out_file}>&-
     
     
     
     
     done




done
done

:<<''

sed -z 's/ans = \n\n    "\n     /John/g; s/"\n\n,/qq,/g;  s/\nans = \n\n    ",/ww,/g; s/"\n/qw,/g;'
sed -z 's/ans = \n\n    "\n     / /g; s/"\n\n,/,/g;  s/\nans = \n\n    ",/,/g;   s/"\n\n//g;'
-z 's/ans = \n\n    "\n     //g; s/"\n\n,/,/g; s/\nans = \n\n    ",/,/g; s/"\n\n//gi'
sed -z 's/ans = \n\n    "\n     /John /g; s/"\n\n,/qq,/g; s/\nans = \n\n    ",/ww,/g; s/"\n\n/ee/gi; s/\[[^]]*\]..\n/QQ/gi' tmp.csv


awk '

BEGIN{
     printf("\ndecomp, L_sol, U_sol,type,online,ofline,sol,sol_agg\n")
     printf("\n%u,",nds_n_req)
     FS=","
     OFS=","
}
{
     for (i = 1; i <= NF; i++){
          sums[i-1]+=$i
          counts[i-1]++
     }

}
END{
     for(i=0;i<NF;i++){
          printf "\n%s%g", (i==0)?"":"," ,sums[i]/counts[i]
     
     }
     printf "\n"
     

}






' q.csv


