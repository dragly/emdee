#!/bin/bash
FILES=data*.xyz
for f in $FILES
do
    (head -n 2 $f && (tail -n +3 $f | sort)) >> $f.xyz
done
