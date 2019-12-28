#! /bin/bash

rm -f ./fort.71 ./ALL_fort.71
rm -f ./total_eddy.dat ./ALL_total_eddy.dat
rm -f ./fort.666 ./ALL_fort.666
rm -f ./flux_control.dat ./ALL_flux_control.dat
rm -f ./lcfs_plot.dat ./ALL_lcfs_plot.dat
rm -f ./fort.21 ./ALL_fort.21
rm -f ./C_control.txt ./ALL_C_control.txt
rm -f ./fort.777 ./ALL_fort.777

exit

#foreach file (./FORT70/fort.70.*)
for file in ./FORT70/fort.70.*
do
  echo $f

  cp $file ./fort.70
#  ./ccs_sa > /dev/null
  ./ccs_sa
  cat ./fort.71 >> ./ALL_fort.71
  cat ./total_eddy.dat >> ./ALL_total_eddy.dat
  cat ./fort.666 >> ./ALL_fort.666
  cat ./flux_control.dat >> ./ALL_flux_control.dat
  cat ./lcfs_plot.dat >> ./ALL_lcfs_plot.dat

done

