#rm art2siesta log.file.1 siesta2art
cp ethane refconfig.dat
sed -i "/Counter:/c\Counter:      1000" filecounter
rm art2siesta CLOCK fdf-*.log INPUT_TMP.* log.file.* output.log siesta2art Si+I.* test_*
rm INPUT_TMP*
rm BASIS* fdf-*
rm NON_TRIMMED_KP_LIST
rm Si.ion  Si.ion.nc  Si.ion.xml
rm FORCE_STRESS  OCCS
rm 0_NORMAL_EXIT
rm DM*
rm fdf.log
rm min1*
rm sad1*
rm art2gaussian.inp gaussian2art.log temp.chk log events.list test.log temp.xyz
