#!/bin/bash

read idseq
#input the sequence (missing first brackets….

for i in ${idseq}; do
	title=“rsync://rsync.cmbi.ru.nl/dssp/${i:1:4}*”
	echo $title
	rsync -avz ${title:1:35}  DSSP/
done
