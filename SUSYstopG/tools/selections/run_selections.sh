#!/bin/bash

dir=/afs/cern.ch/work/m/mewu/public
folder=76X_new
DIR=$dir/$folder

for i in $(ls $DIR/);do
    a=\"$DIR/$i/\"
    #a=\"76X/WJetsToLNu\"
    b=\"$i\"
    test=$DIR/$i
    if [ -d "$test" ];then
	root -l <<EOF
.L ./selections/selections.cxx+
Selections s($a, $b);
s.skimming();
.q
EOF
    fi
done

# a=\"76X/BulkGravToZZToZlepZinv_narrow_1000\"
# b=\"BulkGravToZZToZlepZinv_narrow_1000\"
# root -l <<EOF
#     .L ./selections/selections.cxx+
#     Selections s($a, $b);
#     s.skimming();
#     .q
#     EOF
# EOF
