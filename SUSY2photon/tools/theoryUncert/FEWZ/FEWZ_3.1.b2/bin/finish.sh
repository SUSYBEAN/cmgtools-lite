#!/bin/bash -e

# Usage: finish.sh <run_dir> <output_file> [ <init_sect> <last_sect> <sectloop_step> ]
#    or  finish.sh <run_dir1> <out_file1> +-*/. <run_dir2> <out_file2> <new_out_file>
# Note: output filename <out_file1> <out_file2> must be the full name, i.e. including the order prefix
#       +,-,*,/ or . is the operator between two output files, in which . means average the two

finishusage(){
printf "Usage: `basename $0` <run_dir> <fullname_output_file> \n\
 OR `basename $0` <run_dir1> <out_file1> +-*/. <run_dir2> <out_file2> <new_out_file>\n"
exit 1
}

if [ $# -ne 2 -a $# -ne 5 -a $# -ne 6 ]; then
   finishusage
fi

### required for both single job and double job modes
   rundir1=${1%/}
   rundirname1=${rundir1##*/}
   infile1=$2
   prefix1=${infile1%%.*}
   pdfinfile1=${prefix1}.pdf.${infile1#*.}
   p_infile1=${prefix1}.p_${infile1#*.}
   m_infile1=${prefix1}.m_${infile1#*.}
   if [ $prefix1 = NNLO ]; then
      if [ $# -eq 5 ]; then
         python scripts/merge_parallel.py $rundir1 $infile1 $3 $4 $5
      else
         python scripts/merge_parallel.py $rundir1 $infile1
      fi 
      ### merge_parallel will generate file in current directory if it can't write to the run directory
      ### -nt returns true if the RHS file doesn't exist or is older than the LHS file
      if [ ./$infile1 -nt $rundir1/$infile1 ]; then
         infiledir1='.'
         echo "No/Older merged file(s) found in run directory, use the file(s) in current directory"
      else
         infiledir1=$rundir1
      fi
   else
      ### LO and NLO only have one sector and skip merging, but move scale variation file to one level up
      infiledir1=$rundir1/${rundirname1}0
   fi

### extra merging for the second job
if [ $# -eq 6 ]; then
   rundir2=${4%/}
   rundirname2=${rundir2##*/}
   infile2=$5
   prefix2=${infile2%%.*}
   pdfinfile2=${prefix2}.pdf.${infile2#*.}
   p_infile2=${prefix2}.p_${infile2#*.}
   m_infile2=${prefix2}.m_${infile2#*.}
   if [ $prefix2 = NNLO ]; then
      python scripts/merge_parallel.py $rundir2 $infile2
      if [ ./$infile2 -nt $rundir2/$infile2 ]; then 
         infiledir2='.'
         echo "No/Older merged file(s) found in run directory, use the file(s) in current directory"
      else
         infiledir2=$rundir2
      fi
   else
      ### LO and NLO only have one sector and skip merging, but move scale variation file to one level up
      infiledir2=$rundir2/${rundirname2}0
   fi
fi

### generate file to the current directory for further processing
if [ $# -eq 2 -o $# -eq 5 ]; then
   outfile=$infile1
   pdfoutfile=$pdfinfile1
   p_outfile=$p_infile1
   m_outfile=$m_infile1
   ### directly copy the file for single job mode
   if [ $infiledir1 != '.' ]; then
      cp $infiledir1/$infile1 $outfile
      [ -e $infiledir1/$pdfinfile1 ] && cp $infiledir1/$pdfinfile1 $pdfoutfile
      [ -e $infiledir1/$p_infile1 ] && cp $infiledir1/$p_infile1 $p_outfile
      [ -e $infiledir1/$m_infile1 ] && cp $infiledir1/$m_infile1 $m_outfile
   fi
else
   ### combine two jobs using the operator given in the argument list
   outfile=$6
   prefix=${outfile%%.*}
   ### if the new output file has the same order prefix as the one of the two jobs
   if [ $prefix == $prefix1 ] && [ $prefix1 == $prefix2 ]; then
      pdfoutfile=$prefix.pdf.${outfile#*.}
      p_outfile=$prefix.p_${outfile#*.}
      m_outfile=$prefix.m_${outfile#*.}
   else
      pdfoutfile=pdf.$outfile
      p_outfile=p_$outfile
      m_outfile=m_$outfile
   fi
   python scripts/combine.py $infiledir1/$infile1 $3 $infiledir2/$infile2 $outfile
   if [ -e $infiledir1/$pdfinfile1 ] && [ -e $infiledir2/$pdfinfile2 ]; then
      python scripts/combine.py $infiledir1/$pdfinfile1 $3 $infiledir2/$pdfinfile2 $pdfoutfile
   fi
   if [ -e $infiledir1/$p_infile1 ] && [ -e $infiledir2/$p_infile2 ]; then
      python scripts/combine.py $infiledir1/$p_infile1 $3 $infiledir2/$p_infile2 $p_outfile
   fi
   if [ -e $infiledir1/$m_infile1 ] && [ -e $infiledir2/$m_infile2 ]; then
      python scripts/combine.py $infiledir1/$m_infile1 $3 $infiledir2/$m_infile2 $m_outfile
   fi
fi

### final processing of the merged/combined file
python scripts/get_momentA.py $outfile
if [ -e $p_outfile ] && [ -e $m_outfile ]; then
   python scripts/get_momentA.py $p_outfile
   python scripts/get_momentA.py $m_outfile
   python scripts/do_scales.py $outfile $p_outfile $m_outfile
fi
if [ -e $pdfoutfile ]; then
   python scripts/get_momentA.py $pdfoutfile
   python scripts/do_pdfs.py $outfile $pdfoutfile
fi
rm -f $pdfoutfile $p_outfile $m_outfile


