#!/bin/csh
set list = (`ls *.csh | grep tourne `)
set j = 0
foreach i ($list)
echo "Lance $i"
bsub -q 1nh $i
end
