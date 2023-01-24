#!/bin/bash
#
	echo
	echo " +++ Starting scritp to remap ROMS ssh +++"
	echo
		infile=$1
		outfile=$2
		cmems=$3

		cdo="cdo -s --no_warnings"

		$cdo genbil,$cmems $infile weights_r.nc
		$cdo remap,$cmems,weights_r.nc $infile $outfile


	echo
	echo " +++ End of script +++"
	echo

