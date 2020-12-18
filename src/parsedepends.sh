#!/bin/bash
cp Makefile.depends Makefile.depends.orig
automakefile=Makefile.am

# Remove xx_MOD variable definitions. This file makes this obsolete
# and replacing text in a variable name to form the name of a source
# variable name does not work with automake.
sed -i '/^.._MOD:=/d' Makefile.depends

# Don't use := with automake
sed -i 's/:=/=/g' Makefile.depends

# Similarly, remove the stated .o dependencies. Once the sources are delcared as such, automake does the rest.
sed -i '/^.*\.o:\ \$/d' Makefile.depends

# Now replace all instances of these abbreviations with their definitions
# Also, replace .o with .f; automake knows what to do with .f files
sed -i 's/\$(al_MOD)/allocate_mod.f/g' Makefile.depends
sed -i 's/\$(bk_MOD)/background_mod.f/g' Makefile.depends
sed -i 's/\$(dx_MOD)/degas2_xgc_mod.f/g' Makefile.depends
sed -i 's/\$(de_MOD)/detector_mod.f/g' Makefile.depends
sed -i 's/\$(ef_MOD)/efititp_mod.f/g' Makefile.depends
sed -i 's/\$(el_MOD)/element_mod.f/g' Makefile.depends
sed -i 's/\$(ff_MOD)/flight_frag_mod.f/g' Makefile.depends
sed -i 's/\$(g2_MOD)/geometry2d_mod.f/g' Makefile.depends
sed -i 's/\$(gi_MOD)/geomint_mod.f/g' Makefile.depends
sed -i 's/\$(cm_MOD)/macros_mod.f/g' Makefile.depends
sed -i 's/\$(ma_MOD)/materials_mod.f/g' Makefile.depends
sed -i 's/\$(mp_MOD)/mpi_mod.f/g' Makefile.depends
sed -i 's/\$(ou_MOD)/output_mod.f/g' Makefile.depends
sed -i 's/\$(pm_MOD)/pmi_mod.f/g' Makefile.depends
sed -i 's/\$(pd_MOD)/pmidata_mod.f/g' Makefile.depends
sed -i 's/\$(pf_MOD)/pmiformat_mod.f/g' Makefile.depends
sed -i 's/\$(po_MOD)/postdetector_mod.f/g' Makefile.depends
sed -i 's/\$(pr_MOD)/problem_mod.f/g' Makefile.depends
sed -i 's/\$(ps_MOD)/problem_mod.f/g' Makefile.depends
sed -i 's/\$(ra_MOD)/ratecalc_mod.f/g' Makefile.depends
sed -i 's/\$(rc_MOD)/reaction_mod.f/g' Makefile.depends
sed -i 's/\$(rd_MOD)/reactiondata_mod.f/g' Makefile.depends
sed -i 's/\$(rf_MOD)/readfilenames_mod.f/g' Makefile.depends
sed -i 's/\$(sc_MOD)/sector_mod.f/g' Makefile.depends
sed -i 's/\$(sn_MOD)/snapshot_pdf_mod.f/g' Makefile.depends
sed -i 's/\$(so_MOD)/sources_mod.f/g' Makefile.depends
sed -i 's/\$(sp_MOD)/species_mod.f/g' Makefile.depends
sed -i 's/\$(sa_MOD)/stat_mod.f/g' Makefile.depends
sed -i 's/\$(tl_MOD)/tally_mod.f/g' Makefile.depends
sed -i 's/\$(xi_MOD)/xgc1interp_mod.f/g' Makefile.depends
sed -i 's/\$(xs_MOD)/xsection_mod.f/g' Makefile.depends
sed -i 's/\$(zn_MOD)/zone_mod.f/g' Makefile.depends

# Now , define *_deps varibles to store unique lists of sources
bins=(`grep bin_PROGRAMS $automakefile | sed 's/bin_PROGRAMS.*=\(.*\)/\1/g'`)

echo " " >> Makefile.depends

for bin in "${bins[@]}"; do
	echo "bin=$bin"
	predeps=(`grep "^${bin}_SOURCES\ +=\ $bin" $automakefile | sed "s/${bin}_SOURCES.*=\(.*\)/\1/g"`)
	echo predeps="${predeps[@]}"

	deps=()
	for dep in "${predeps[@]}"; do
		depstrip=`echo $dep|sed 's/\(.*\)\.f/\1/g'`
		echo dep=$dep
		echo depstrip=$depstrip
		mods=(`grep "^${depstrip}_mods=" Makefile.depends | sed "s/${depstrip}_mods=\ \(.*\)$/\1/g"`)
		echo mods="${mods[@]}"
		deps+=("${mods[@]}")
	done
	echo "${deps[@]}"
	IFS=$'\n'
	sorted=($(sort <<<"${deps[*]}" | uniq))
	unset IFS
	echo "${bin}_moddeps = ${sorted[*]}" >> Makefile.depends
done
