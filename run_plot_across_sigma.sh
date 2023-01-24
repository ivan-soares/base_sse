#!/bin/bash
#

	echo
	echo " +++ Starting script to plot across shelf velocities +++"
	echo

	romsdir="/home/ivans/roms/cases/brz/brz0.05r/expts/d-storage/monmean"
	gridfile='grid_brz0.05r_01d.nc'


	for year in 2016 2017 2018 2019; do

		romsfile="$romsdir/roms_monmean_brz0.05r_01d_${year}0101_regridded_sigma.nc"

		for mon in 1 7; do

			echo; echo " ... doing year $year , month $mon"; echo

			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year 607   5   72  'Sao Chico'
			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year 540   8   72  'Salvador'
			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year 365  75  167  'Cabo S. Tomé'
			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year 340  71  157  'Cabo Frio'
			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year 283  27  114  'Ilha Bela'
			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year 187  28  120  'Floripa'
			python plot_roms_across_shelf_monmean_two_layers.py $romsfile $mon $year  96  30  122  'RGS'

		done
	done

	echo
	echo " ++= End of script +++"
	echo
