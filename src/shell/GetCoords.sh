function GetTSSCoords() {
  zcat $1 | awk 'BEGIN{OFS="\t"}{print $2,$3,$3+1,$1,1000,"+"}' | tail -n+2 | sort -k1,1V -k2,2g | bgzip > $2
  tabix -p bed $2
}
function GetSNPCoords() {
  if [[ $1 -eq 1 ]]; then
      zcat $2 | cut -f1 | tail -n+2 > $3
  fi
  if [[ $1 -eq 2 ]]; then
     cat ${@:3} | awk -F"_" '{print $1"\t"$2"\t"$2+1"\t"$0}' | sort -k1,1g -k2,2g | uniq | bgzip > $2
     tabix -p bed $2
  fi
}
