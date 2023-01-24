#!/bin/bash
#

	echo
	echo " +++ Starting script to plot CMEMS surface maps +++"
	echo

	offset="roms_brz0.05r_01d_mean_zeta_remaped.nc"

	#### make a list of file names

	year=$1

	ls -1 data/cmems_sla+vels_${year}*.nc > list_cmems 

	vint=5
	vscl=10.
	wesn="-46 -33 -30 -13"

	while read line; do
		echo $line
		python plot_cmems_surf_maps.py $line $offset $vscl
	done < list_cmems

	echo
	echo " +++ End of script +++"
	echo
