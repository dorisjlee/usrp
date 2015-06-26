#!/bin/bash

usage () {
  echo "Usage:  join_vtk.sh -i <INBASENAME> -o <OUTBASENAME> -v <VAR_ID> -c <START:SKIP:END>"
  echo "  e.g., join_vtk.sh -i RadParGrav -o RadParGrav_joined -v starpar -c 0:5:100"
  exit 1
}

# Specify the location of the Athena installation directory
#homedir="/Users/askinner/Documents/teaching/usrp/hpc/athena4.2/vis/vtk"
homedir="/Users/dorislee/Desktop/Princeton_USRP/usrp/hpc_bootcamp/day3/hpc"
# Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
   echo "You must specify at least 1 argument!"
   usage
fi

# Set defaults 
outbasename=""
inbasename=""
varid=""
cycle="0:0"

# Parse the arguments
while getopts hi:o:v:c: opt
do
   case "$opt" in
      h) usage;;
      i) inbasename=$OPTARG;;
      o) outbasename=$OPTARG;;
      v) varid=$OPTARG;;
      c) cycle=$OPTARG;;
      \?) usage;;
   esac
done

# Parse in/out basenames and variable ID
if [[ $inbasename == "" ]]; then
  echo "Please specify an input basename!"
  usage
elif [[ $outbasename == "" ]]; then
  outbasename=${inbasename}_joined
fi
echo $inbasename $outbasename

# Parse the cycle variable
if [[ ${cycle##*:*} != "" ]]; then
  usage
else
  s=${cycle%%:*}
  e=${cycle##*:}
  k=${cycle#*:}
  if [[ ${k##*:*} != "" ]]; then
    k="1"
  else
    k=${k%:*}
    if [[ ${k##*:*} == "" ]]; then
      usage
    fi
  fi
fi

if [[ $varid != "" ]]; then
  varid=.${varid}
fi
filepattern=${inbasename}*${varid}.????.vtk
#filepattern=${inbasename}*.????${varid}.vtk
echo "filepattern = $filepattern"

declare -a proclist
proclist=( `find id* -name "$filepattern" | cut -d '/' -f 1 | cut -d 'd' -f 2 | sort | uniq` )
nproc=${#proclist[@]}
echo "nproc = $nproc"
echo "proclist = $proclist"
#echo `find id* -name "$filepattern"`

gcc -o $homedir/join_vtk $homedir/join_vtk.c -lm

#nproc=`ls -ld id* | wc -l`
#nproc=$((nproc-1))
#echo $nproc
for c in `seq $s $k $e`
do
  num=$(printf "%04d" "$c")
  arglist="id0/$inbasename$varid.$num.vtk"
#  arglist="id0/$inbasename.$num$varid.vtk"

  for n in `seq 0 $((nproc-1))`
  do
    p=${proclist[$n]}
#    echo "loop $n for proc $p"
    if [[ $p != 0 ]]; then
      arglist="$arglist id$p/$inbasename-id$p$varid.$num.vtk"
#      arglist="$arglist id$p/$inbasename-id$p.$num$varid.vtk"
    fi
  done

  outfilename=$outbasename$varid.$num.vtk
#  outfilename=$outbasename.$num$varid.vtk
  echo "Cycle $c:"
  echo "  Stitching $nproc processors together to form $outfilename"
  echo $arglist
  $homedir/join_vtk -o $outfilename $arglist
done



