#!/bin/bash

read idseq

#must have no [] in list!
#$idseq2=${idseq:1:${#idseq}-1}

for i in $idseq; do
	echo rsync://rsync.cmbi.ru.nl/ddsp/${i}
	rsync -avz rsync://rsync.cmbi.ru.nl/ddsp/${i+"dssp"}  DSSP/
done
