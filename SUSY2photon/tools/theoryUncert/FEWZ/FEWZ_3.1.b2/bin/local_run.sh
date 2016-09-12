#!/bin/bash
# Usage: local_run.sh z/w <run_dir> <input_file> <histo_file> <output_file_extension> <pdf_dir> [ <num_proc> ]
#    or  local_run.sh z/w <run_dir> <input_file> <histo_file> <output_file_extension> <pdf_dir> [ <num_proc> <init_sect> <last_sect> <sectloop_step> ]
#    or  local_run.sh z/w <run_dir> <input_file> <histo_file> <output_file_extension> <pdf_dir> [ <num_proc> <which_sect> ]
#  Note: arguments in [] are optional and for advanced users only
#        <num_proc> defaults to 1 if not given
#        <init_sect> <last_sect> <sectloop_step> specifies the sectors desired to run through the loop:
#                                                FOR <sect_num> FROM <init_sect> TO <last_sect> STEP <sectloop_step>
#        <which_sect> specifies the one sector to user want to run

localusage(){
echo "Usage: `basename $0` w/z <run_dir> <input_file> <histo_file> <output_file_extension> <pdf_dir> [ <num_proc> ]"
echo " Note: <num_proc> defaults to 1"
exit 1
}
#[ $# -lt 6 ] && localusage \
#             || printf "Running Directory: $2\nInput Setting: $INFILE\nHistogram File: $4\nOutput File: $5\nPDF Directory $6\nNo. of Processors: $7\n"
[ $# -lt 6 ] && localusage

BOSON=$1
INFILE=$3
OUTFILE=$5
HISTFILE=$4
RUNDIR=${2%/}
PDFDIR=${6%/}
### Turn the directory into relative directory just in case
RUNDIR=`python scripts/get_relpath.py $RUNDIR ./`
PDFDIR=`python scripts/get_relpath.py $PDFDIR ./`
CRTDIR=$PWD

### Have hacked Cuba file $CUBADIR/src/common/Fork.c by requiring t->ncores = 0 
### Set environmental variable just in case User re-installed unhacked Cuba
export CUBACORES=0 #makes sure Cuba doesn't run in parallel

# Read in necessary info
NUM=0
QUEUE=""
#NUM_CORE=$(cat /proc/cpuinfo | grep 'model name' | wc -l) ### checking for cpu numbers
SECTORS=`python scripts/get_sects.py $INFILE $BOSON`
SCALEVAR=`python scripts/get_scalevar.py $INFILE`

# check first argument, to make sure supported, and set executable
if [ $BOSON = "z" ] || [ $BOSON = "w" ]; then
   EXEC=fewz$BOSON
else
   echo "Unrecognized argument; defaulting to neutral current."
   EXEC=fewzz
fi

# default values of optional arguments if not given
MAX_PROC=1
INIT_SECT=1
LAST_SECT=$SECTORS
SECT_STEP=1

# read in optional arguments if given
[ $# -ge 7 ] && MAX_PROC=$7
[ $# -ge 8 ] && INIT_SECT=$8
[ $# -eq 8 ] && LAST_SECT=$INIT_SECT
[ $# -ge 9 ] && LAST_SECT=$9
[ $# -ge 10 ] && SECT_STEP=${10}

# check the optional input arguments
#[ $MAX_PROC -gt $NUM_CORE ] && MAX_PROC=$NUM_CORE ### checking for cpu numbers
if [ $SECT_STEP -eq 0 ] || \
   [ $SECT_STEP -gt 0 -a $INIT_SECT -gt $LAST_SECT ] || \
   [ $SECT_STEP -lt 0 -a $LAST_SECT -gt $INIT_SECT ]; then
   exit 2
else
   if [ $SECT_STEP -gt 0 ]; then
      if [ $INIT_SECT -gt $SECTORS ] || [ $LAST_SECT -lt 1 ]; then
         exit 2
      fi
   else
      if [ $LAST_SECT -gt $SECTORS ] || [ $INIT_SECT -lt 1 ]; then
         exit 2
      fi
   fi
fi
[ $INIT_SECT -lt 1 ] && INIT_SECT=1
[ $INIT_SECT -gt $SECTORS ] && INIT_SECT=$SECTORS
[ $LAST_SECT -lt 1 ] && LAST_SECT=1
[ $LAST_SECT -gt $SECTORS ] && LAST_SECT=$SECTORS

# function to update info on how many sectors are running
regenqueue(){
   oldrequeue=$QUEUE
   QUEUE=""
   NUM=0
   for PID in $oldrequeue
   do
      numproc=$(ps -p $PID | wc -l)
      if [ $numproc -gt 1 ] ; then
         QUEUE="$QUEUE $PID"
         NUM=$(($NUM+1))
      fi
   done
}

# function to check if any one of the running sector finishes
checkqueue(){
   oldchqueue=$QUEUE
   for PID in $oldchqueue
   do
      numproc=$(ps -p $PID | wc -l)
      if [ $numproc -le 1 ] ; then
         regenqueue
         break
      fi
   done
}

# create sector running directory and set up the structure
# ps. last three arguments for the python script below are optional and leaving them out will always be correct
python scripts/create_parallel.py $BOSON $INFILE $RUNDIR $INIT_SECT $LAST_SECT $SECT_STEP
if ! [ -e $RUNDIR/$EXEC ] ; then
    cp $EXEC $RUNDIR
fi
if ! [ -e $RUNDIR/$HISTFILE ] ; then
    python scripts/get_bin_files.py $HISTFILE $RUNDIR
fi
#if ! [ -e $RUNDIR/$INFILE ] ; then ### always copy input file, in case changed
cp $INFILE $RUNDIR
#fi
cd $RUNDIR

echo "Using "$MAX_PROC" parallel threads"
STARTTIME=`date +%s`

# loop through all sectors to finish job
#i=0
#while [ $i -lt $SECTORS ]
#do
#   echo "Starting sector $(($i+1))"
#   cd ${RUNDIR##*/}$i
#   #../$EXEC -i $INFILE -h ../$HISTFILE -o $OUTFILE -p ../../$PDFDIR > screen.out 2>&1 & ### tcsh incompatible
#   ../$EXEC -i $INFILE -h ../$HISTFILE -o $OUTFILE -p ../../$PDFDIR > screen.out &
#   PID=$!
#   QUEUE="$QUEUE $PID"
#   NUM=$(($NUM+1))
#   cd ..
#   while [ $NUM -ge $MAX_PROC ]; do
#      checkqueue
#      sleep 1
#   done
#   i=$(($i+1))
#done

# loop through all till finishing the job
# loop of scale variation
j=1
while [ $j -le $SCALEVAR ]
do
# loop of specified sectors
i=$INIT_SECT
if [ $SECT_STEP -gt 0 ] ; then
   # loop from small sector to large sector
   # for (( i=$INIT_SECT ; i<=$LAST_SECT ; i=i+$SECT_STEP ))
   while [ $i -le $LAST_SECT ]
   do
      cd ${RUNDIR##*/}$(($i-1))
      if [ $j -eq 1 ] ; then
         echo "Starting sector $i ..."
         #../$EXEC -i $INFILE -h ../$HISTFILE -o $OUTFILE -p ../../$PDFDIR > screen.out 2>&1 & ### tcsh incompatible
         ../$EXEC -i ../$INFILE -h ../$HISTFILE -o $OUTFILE -p ../../$PDFDIR -s $(($i-1)) > screen.out &
      elif [ $j -eq 2 ] ; then
         cd pscale
         echo "Starting sector $i for up scale variation ..."
         ../../$EXEC -i ../../p_$INFILE -h ../../$HISTFILE -o p_$OUTFILE -p ../../../$PDFDIR -l .. -s $(($i-1)) > p_screen.out &
      elif [ $j -eq 3 ] ; then
         cd mscale
         echo "Starting sector $i for down scale variation ..."
         ../../$EXEC -i ../../m_$INFILE -h ../../$HISTFILE -o m_$OUTFILE -p ../../../$PDFDIR -l .. -s $(($i-1)) > m_screen.out &
      fi
      PID=$!
      QUEUE="$QUEUE $PID"
      NUM=$(($NUM+1))
      cd $CRTDIR/$RUNDIR
      while [ $NUM -ge $MAX_PROC ]; do
         checkqueue
         sleep 3
      done
      i=$(($i+$SECT_STEP))
   done
else
   # loop from large sector to small sector
   # for (( i=$INIT_SECT ; i>=$LAST_SECT ; i=i+$SECT_STEP ))
   while [ $i -ge $LAST_SECT ]
   do
      cd ${RUNDIR##*/}$(($i-1))
      if [ $j -eq 1 ] ; then
         echo "Starting sector $i ..."
         #../$EXEC -i $INFILE -h ../$HISTFILE -o $OUTFILE -p ../../$PDFDIR > screen.out 2>&1 & ### tcsh incompatible
         ../$EXEC -i ../$INFILE -h ../$HISTFILE -o $OUTFILE -p ../../$PDFDIR -s $(($i-1)) > screen.out &
      elif [ $j -eq 2 ] ; then
         cd pscale 
         echo "Starting sector $i for up scale variation ..."
         ../../$EXEC -i ../../p_$INFILE -h ../../$HISTFILE -o p_$OUTFILE -p ../../../$PDFDIR -l .. -s $(($i-1)) > p_screen.out &
      elif [ $j -eq 3 ] ; then
         echo "Starting sector $i for down scale variation ..."
         cd mscale
         ../../$EXEC -i ../../m_$INFILE -h ../../$HISTFILE -o m_$OUTFILE -p ../../../$PDFDIR -l .. -s $(($i-1)) > m_screen.out &
      fi
      PID=$!
      QUEUE="$QUEUE $PID"
      NUM=$(($NUM+1))
      cd $CRTDIR/$RUNDIR
      while [ $NUM -ge $MAX_PROC ]; do
         checkqueue
         sleep 0.5
      done
      i=$(($i+$SECT_STEP))
   done
fi
j=$(($j+1))
done
wait

ENDTIME=`date +%s`
DIFFTIME=$(( $ENDTIME - $STARTTIME ))
echo "Job finished in $DIFFTIME seconds"
echo "Run the following to post-process the output files: finish.sh $RUNDIR <order>.$OUTFILE"
