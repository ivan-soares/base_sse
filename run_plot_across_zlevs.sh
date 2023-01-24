#!/bin/bash
#

	echo
	echo " +++ Starting script to plot across shelf velocities +++"
	echo

	romsdir="/home/ivans/roms/cases/brz/brz0.05r/expts/d-storage/monmean"

	for year in 2015 2016 2017 2018 2019; do

		romsfile="$romsdir/roms_monmean_brz0.05r_01d_${year}0101_regridded_zlevs.nc"

		for mon in 1 7; do

			echo; echo " ... doing year $year , month $mon"; echo

			python plot_roms_across_shelf_monmean_z_two_layers.py $romsfile $mon $year -11. -37.0 -33. -34.0 1500. 'Sao Chico'
			python plot_roms_across_shelf_monmean_z_two_layers.py $romsfile $mon $year -22. -41.5 -34. -39.5  500. 'C. S Tomé' 
			python plot_roms_across_shelf_monmean_z_two_layers.py $romsfile $mon $year -24. -46.5 -36. -41.0  500. 'B. Santos'
			python plot_roms_across_shelf_monmean_z_two_layers.py $romsfile $mon $year -26. -49.5 -37. -42.0  600. 'Floripa'
			python plot_roms_across_shelf_monmean_z_two_layers.py $romsfile $mon $year -28. -50.0 -39. -44.0  600. 'S. Marta'
			python plot_roms_across_shelf_monmean_z_two_layers.py $romsfile $mon $year -30. -51.0 -40. -44.0  600. 'Rio Grande'

		done
	done

	echo
	echo " ++= End of script +++"
	echo
