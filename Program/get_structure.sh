#!/bin/bash

read idseq

for i in $idseq; do
	rsync -avz rsync://rsync.cmbi.ru.nl/dssp/${i}* /DSSP/
done
