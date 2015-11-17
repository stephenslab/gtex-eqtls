source ~/.bashrc
function GetTSSCoords() {
  zcat $1 | awk 'BEGIN{OFS="\t"}{print $2,$3,$3+1,$1,1000,"+"}' | tail -n+2 | sort -k1,1V -k2,2g | bgzip > $2
  tabix -p bed $2
}
function GetSNPCoords() {
  if [[ $3 -eq 1 ]]; then
     mkdir -p SNPCoordsTmp
     for i in `ls $1`; do
         echo -e '#!/bin/bash\nzcat '"$1/$i"' | cut -f1 | tail -n+2 > SNPCoordsTmp/'"$i"'.txt' # | sbatch -J ID."$i" -o Extract."$i".o%j
     done
  fi
  if [[ $3 -eq 2 ]]; then
     cat SNPCoordsTmp/* | awk -F"_" '{print $1"\t"$2"\t"$2+1"\t"$0}' | sort -k1,1g -k2,2g | uniq | bgzip > $2
     tabix -p bed $2
  fi
}