#!/bin/sh


for file in list*.txt;
do
 grep "*        "  $file | awk {'print $4","$6","$8'}  > new_${file}
done
