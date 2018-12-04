#!/bin/bash

PREFIX=WD40
ASIZE=42
NCLUS=175

mkdir -p output

sed -n '/	LOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER$/,/^AVSFLD/p' ${PREFIX}.dlg > tmp  #get the best conformation of each cluster from dlg file

for ii in $(seq 1 $NCLUS)
do
	sed -n '/^USER    Cluster Rank = '$ii'/,/^ENDMDL/p' tmp > ${PREFIX}_cluster_${ii}
	echo -e "REMARKS Cluster Rank = $ii" > Ligand_${ii}.pdb
	grep '^ATOM' ${PREFIX}_cluster_${ii} >> Ligand_${ii}.pdb
	echo -e "TER" >> Ligand_${ii}.pdb
done

mv ${PREFIX}_cluster_* Ligand_*.pdb output