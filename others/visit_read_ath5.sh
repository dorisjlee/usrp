#!/bin/bash
#Bash Script used for converting .ath5 to .ath file so that it could be read by VisIT's SmartGrouping Feature
set LANG=C #Might need to set regex sed encoding for MacOSX 
sed -i --  's/ath5/ath/g' *.xdmf    # globally, substitute all instances of 'ath5' with 'ath' inside  the .xdmf  files
rm *.xdmf-- #Remove all the temp files
for file in *.ath5
do
     new=${file/.ath5/.ath}
     mv $file $new 
done

for file in *.ath5.xdmf
do
     new=${file/.ath5.xdmf/.ath.xdmf}
     mv $file $new
done
