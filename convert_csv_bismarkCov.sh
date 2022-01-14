for f in *.csv.gz; do fname=${f%".csv.gz"}; zcat $f | awk -vOFS='\t' '{ if (NR>1 && $4>0) print($1,$2,$2,$3/$4,$3,$4-$3) }' > ${fname}.bed; done
