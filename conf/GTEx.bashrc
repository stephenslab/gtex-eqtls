module load libtool/2.4
module load zlib/1.2
module load gsl/1.16 
source ~/.bashrc
export DataDir=/project/mstephens/gtex/analysis/april2015/eqtl_data
export DataPrefix=$DataDir/GTEx_Analysis_2015-01-12_
export InputDir=/project/mstephens/gtex/analysis/april2015/A.Input
export ArchiveDir=/project/mstephens/gtex/analysis/april2015/1.Archive
export TmpDir=/scratch/midway/gaow
export ConfDir=$HOME/GIT/gtex-eqtls/conf
export SrcDir=$HOME/GIT/gtex-eqtls/src
export LogDir=/project/mstephens/gtex/analysis/april2015/0.Log
export DBDir=/project/mstephens/database
export BFDir=/project/mstephens/gtex/analysis/april2015/C.eqtlbma_hm/RawBFs
function GetTSSCoords() {
  zcat $1 | awk 'BEGIN{OFS="\t"}{print $2,$3,$3+1,$1,1000,"+"}' | tail -n+2 | sort -k1,1V -k2,2g | bgzip > $2
  tabix -p bed $2
}
function GetSNPCoords() {
  if [[ $3 -eq 1 ]]; then
     mkdir -p SNPCoordsTmp
     for i in `ls $1`; do
         echo -e '#!/bin/bash\nzcat '"$1/$i"' | cut -f1 | tail -n+2 > SNPCoordsTmp/'"$i"'.txt' | sbatch -J ID."$i" -o $LogDir/Extract."$i".o%j
     done
  fi
  if [[ $3 -eq 2 ]]; then
     cat SNPCoordsTmp/* | awk -F"_" '{print $1"\t"$2"\t"$2+1"\t"$0}' | sort -k1,1g -k2,2g | uniq | bgzip > $2 
     tabix -p bed $2
  fi
}