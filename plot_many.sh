#!/bin/bash
#

	echo
	echo " +++ Starting script to plot ROMS surface maps +++"
	echo

	year=$1

	romsfile="roms_brz0.05r_01d_${year}0101_surf_daily.nc"
	offset="roms_brz0.05r_01d_mean_zeta.nc"

	vint=5
	vscl=15.
	wesn="-46 -33 -30 -13"
	wesn="1 218 120 553" 

	for d in {4..364..5}; do
		echo $d
		python plot_roms_surf_maps.py $romsfile $offset $d $vint $vscl $wesn
	done

	echo
	echo " +++ End of script +++"
	echo
