#!/bin/bash
if [ $# -lt 3 ]; then
    echo "usage: zcat snp_coords.bed.gz | bash $0 /project/mstephens/database/dbSNP/All_20150605.vcf.gz tss_coords.bed.gz 100000 | gzip > snp-gene-pairs.gz"
    exit 0
fi
while IFS='' read -r line || [[ -n $line ]]; do
	coord=`echo $line | awk '{print $1":"$2"-"$2}'`
	varid=`echo $line | awk '{print $4}'`
	rsid=`tabix $1 $coord | awk '{print $1"_"$2"_"$4"_"$5"_b37""\t"$3}' | grep -w $varid | cut -f2`
	cisg=`bedtools window -w $3 -a <(echo $line | sed 's/ /\t/g') -b $2 | cut -f8`
    if [[ -z $cisg ]]; then
        cisg='-'
    fi
	if [[ -n $rsid ]]; then
		rsid=`echo $rsid | xargs -n1 | sort -u | xargs | sed 's/ /,/g'`
		cisg=`echo $cisg | xargs -n1 | sort -u | xargs | sed 's/ /,/g'`
		echo -e "$varid"'\t'"$rsid"'\t'"$cisg"
	fi
done
