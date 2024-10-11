#!/bin/bash 

GENES=$1
TABFILE=$2
OUTFILE=$3

echo "Gene" | paste - <( head -1 $TABFILE) > $OUTFILE
while read gene; do
	grep $gene $TABFILE | awk -F"\t" '{print gn,$0}' gn=$gene >> $OUTFILE
done < $GENES
