#!/bin/bash -e
### Check sector subdirectory structure created by "create_parallel.py"
### Usage: $0 <run_dir> <output_file>
###    or  $0 <run_dir> <output_file> <ifcheck_condorerr> <ifcheck_outputexist> <ifcheck_pdfoutputexist>
### Note: <outputfile> should be the expected full name of the merged file, i.e. in the form of <order>.<filename_extension>
###       <ifcheck_condorerr> <ifcheck_outputexist> <ifcheck_pdfoutputexist> are optional
###       and take value 1 or 0 marking whether checking for them"

### <>

checkusage(){
echo "Usage: `basename $0` <run_dir> <output_file>"
echo "   or  `basename $0` <run_dir> <output_file> <ifcheck_condorerr> <ifcheck_outputexist> <ifcheck_pdfoutputexist>"
echo " Note: <ifcheck_condorerr> <ifcheck_outputexist> <ifcheck_pdfoutputexist> are optional"
echo "       and take value 1 or 0 marking whether checking for them"
exit 1
}
[ $# -lt 2 ] && checkusage

### default to check them all
CondorErr=1
OutputExs=1
PDFoutExs=1
[ $# -gt 2 ] && CondorErr=$3
[ $# -gt 3 ] && OutputExs=$4
[ $# -gt 4 ] && PDFoutExs=$5

cd $1
DIRNAME=${1%/}
DIRNAME=${DIRNAME##*/}
ERRCOUNT=0
for SubDir in $(ls) ; do
    if [ -d $SubDir ] ; then
	cd $SubDir
        if [ $CondorErr -eq 1 ] && [ -s condor_error.err ]; then
	    ERRCOUNT=$(($ERRCOUNT+1))
            echo "Found condor error message in ${1%/}/$SubDir/"
        fi
	if [ $OutputExs -eq 1 ] || [ $PDFoutExs -eq 1 ] ; then
	    WhichSect=$(echo $SubDir | sed -e "s/$DIRNAME//gI")
	    WhichSect=$(($WhichSect+1))
        fi
        if [ $OutputExs -eq 1 ] && [ ! -s $(echo $2 | sed -e "s/LO./LO.sect$WhichSect./1") ] ; then
            ERRCOUNT=$(($ERRCOUNT+1))
            echo "No output file in ${1%/}/$SubDir/"
        fi
        if [ $PDFoutExs -eq 1 ] && [ ! -s $(echo $2 | sed -e "s/LO./LO.sect$WhichSect.pdf./1") ] ; then
            ERRCOUNT=$(($ERRCOUNT+1))
            echo "No pdf output file in ${1%/}/$SubDir/"
        fi
	cd ..
    fi
done

echo $ERRCOUNT "Errors or Missing Files Found"
