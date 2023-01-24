#!/bin/bash
#

	echo
	echo " +++ Starting script to plot across shelf velocities +++"
	echo

	gridfile='grid_brz0.05r_01d.nc'
	romsdir="/home/ivans/roms/cases/brz/brz0.05r/expts/d-storage/monmean"

	for year in 2015 2016 2017 2018 2019; do

		romsfile="$romsdir/roms_monmean_brz0.05r_01d_${year}0101_regridded_zlevs.nc"

		for mon in 1 3 7 9 12; do

			echo; echo " ... doing year $year , month $mon"; echo
			python plot_roms_across_shelf_monmean_ax97.py $romsfile $gridfile $mon $year

		done
	done

	echo
	echo " ++= End of script +++"
	echo
