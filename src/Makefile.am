# The following section between the asterisks is parsed by parsedepends.sh (part of the "make depend" target)
# ************************************************************************************************************
# Particular syntax ought to be followed in this section: whitespace must be on either side of the equals sign
# *_SOURCES must be defined on a single line
# Do not include sysdep.f, string.f, or modules. These are added later by the $TARGET_deps variables defined in Makefile.depends
bin_PROGRAMS = randomtest nrandomtest flighttest datasetup degas2_xgc poly_to_tri ucd_plot adaswrite efit2dg2d geomtest geomtesta matchout matcheir acc_output
bin_PROGRAMS += reactionwrite pmiwrite problemsetup geomread cramdtest cramdtesta degasgeom test readgeometry readgeometry2d definegeometry2d defineback boxgen
bin_PROGRAMS += plotflight cramdproblem readfilenames readbackground updatebackground reactiontest pmitest sourcetest voigttest voigtspec ratecalc dataexam outputbrowser
bin_PROGRAMS += postdetector datamatch tallysetup allocate array teststring elementlookup speciesetup netcdf
randomtest_SOURCES = randomtest.f random.f
nrandomtest_SOURCES = nrandomtest.f nrandom.f
flighttest_SOURCES = flighttest.f doflights.f elementread.f pmimatread.f geometry.f random.f stat.f randother.f reaction.f interpolate.f geomfront.f flight.f sources.f collide.f readfilenames.f plate.f degasinit.f score.f conversions.f voigt.f writeback.f aux2dsetup.f geomsetup.f geomfaces.f geomutils.f allocate.f composite.f
datasetup_SOURCES = datasetup.f elementread.f pmimatread.f elementsetup.f speciesetup.f reactionsetup.f materialsetup.f pmisetup.f readfilenames.f 
degas2_xgc_SOURCES = degas2_xgc.f doflights.f elementread.f pmimatread.f geometry.f random.f stat.f randother.f reaction.f interpolate.f geomfront.f flight.f sources.f collide.f readfilenames.f plate.f degasinit.f score.f conversions.f voigt.f writeback.f aux2dsetup.f geomsetup.f geomfaces.f geomutils.f allocate.f composite.f
poly_to_tri_SOURCES = poly_to_tri.f geometry.f allocate.f readfilenames.f 
ucd_plot_SOURCES = ucd_plot.f geometry.f allocate.f readfilenames.f elementread.f pmimatread.f  reaction.f interpolate.f random.f degasinit.f
adaswrite_SOURCES = adaswrite.f adas_xxdata_11.f elementread.f reaction.f interpolate.f readfilenames.f random.f
efit2dg2d_SOURCES = bspline.mod bspline90_22.o efit2dg2d.f aux2dsetup.f psiinterp.f interpolate.f
geomtest_SOURCES = geomtest.f geometry.f geomsetup.f allocate.f geomutils.f geomfaces.f composite.f string.f readfilenames.f geomfront.f dediagsetup.f
geomtesta_SOURCES = geomtesta.f geometry.f plotgeom.f elementread.f pmimatread.f readfilenames.f reaction.f interpolate.f collide.f randother.f random.f conversions.f voigt.f degasinit.f 
matchout_SOURCES = matchout.f pmimatread.f elementread.f geometry.f readfilenames.f degasinit.f
matcheir_SOURCES = matcheir.f pmimatread.f elementread.f geometry.f readfilenames.f degasinit.f
acc_output_SOURCES = acc_output.f pmimatread.f elementread.f geometry.f readfilenames.f degasinit.f conversions.f
reactionwrite_SOURCES = reactionwrite.f
pmiwrite_SOURCES = pmiwrite.f interpolate.f random.f readfilenames.f elementread.f 
problemsetup_SOURCES = problemsetup.f elementread.f pmimatread.f reaction.f interpolate.f readfilenames.f plate.f geometry.f random.f randother.f
geomread_SOURCES = geomread.f geomsetup.f allocate.f geomutils.f geomfaces.f geometry.f composite.f string.f elementread.f pmimatread.f
cramdtest_SOURCES = cramdtest.f cramdadd.f readfilenames.f elementread.f pmimatread.f geometry.f degasinit.f
cramdtesta_SOURCES = cramdtesta.f cramdadd.f pmimatread.f
degasgeom_SOURCES = degasgeom.f geomsetup.f allocate.f geomutils.f geomfaces.f geometry.f composite.f random.f
test_SOURCES = test.f geomsetup.f allocate.f geomutils.f geomfaces.f geometry.f composite.f
readgeometry_SOURCES = readgeometry.f geomsetup.f allocate.f geomutils.f geomfaces.f geometry.f composite.f string.f readdegas.f readfilenames.f elementread.f pmimatread.f dediagsetup.f aux2dsetup.f def2ddetector.f
definegeometry2d_SOURCES = definegeometry2d.f string.f geomsetup.f allocate.f composite.f geomutils.f geomfaces.f geometry.f aux2dsetup.f pmimatread.f dediagsetup.f readfilenames.f plotgeom.f def2ddetector.f
defineback_SOURCES = defineback.f string.f writeback.f pmimatread.f geometry.f elementread.f geomsetup.f aux2dsetup.f allocate.f geomutils.f geomfaces.f readfilenames.f reaction.f interpolate.f random.f composite.f def2dplasma.f psiinterp.f bspline.mod bspline90_22.o 
boxgen_SOURCES = boxgen.f readfilenames.f geomsetup.f allocate.f composite.f geomutils.f geomfaces.f geometry.f pmimatread.f elementread.f writeback.f cramdadd.f dediagsetup.f reaction.f interpolate.f random.f aux2dsetup.f
plotflight_SOURCES = plotflight.f string.f
cramdproblem_SOURCES = cramdproblem.f pmimatread.f cramdadd.f readfilenames.f elementread.f geometry.f
readfilenames_SOURCES =  readfilenames.f readfilenames.f
updatebackground_SOURCES = updatebackground.f elementread.f readfilenames.f geomsetup.f allocate.f writeback.f pmimatread.f random.f geomutils.f geometry.f geomfaces.f reaction.f interpolate.f
reactiontest_SOURCES = reactiontest.f elementread.f reaction.f interpolate.f readfilenames.f random.f randother.f collide.f pmimatread.f geometry.f
pmitest_SOURCES = pmitest.f elementread.f pmimatread.f plate.f readfilenames.f random.f randother.f geometry.f interpolate.f
sourcetest_SOURCES = sourcetest.f sources.f readfilenames.f randother.f random.f pmimatread.f elementread.f geometry.f cramdadd.f interpolate.f score.f collide.f reaction.f conversions.f
voigttest_SOURCES = voigttest.f voigt.f
voigtspec_SOURCES = voigtspec.f voigt.f
ratecalc_SOURCES = ratecalc.f aladdin_subr.o readfilenames.f interpolate.f
dataexam_SOURCES = dataexam.f readfilenames.f elementread.f
outputbrowser_SOURCES = outputbrowser.f readfilenames.f degasinit.f pmimatread.f elementread.f geometry.f
postdetector_SOURCES = postdetector.f readfilenames.f degasinit.f pmimatread.f elementread.f geometry.f dediagsetup.f def2ddetector.f
datamatch_SOURCES = datamatch.f 
tallysetup_SOURCES = tallysetup.f elementread.f pmimatread.f readfilenames.f string.f geometry.f
allocate_SOURCES = allocate.f
array_SOURCES = 
teststring_SOURCES =
elementlookup_SOURCES =
readbackground_SOURCES = readbackground.f readdegas.f geometry.f elementread.f readfilenames.f geomsetup.f allocate.f reaction.f interpolate.f composite.f geomutils.f geomfaces.f geometry.f writeback.f pmimatread.f cramdadd.f random.f aux2dsetup.f
speciesetup_SOURCES = elementread.f 
netcdf_SOURCES = netcdf.f
# ************************************************************************************************************

include Makefile.depends

randomtest_SOURCES += $(randomtest_deps)
nrandomtest_SOURCES += $(nrandomtest_deps)
flighttest_SOURCES += $(flighttest_deps)
datasetup_SOURCES += $(datasetup_deps)
degas2_xgc_SOURCES += $(degas2_xgc_deps)
poly_to_tri_SOURCES += $(poly_to_tri_xgc_deps)
ucd_plot_SOURCES += $(ucd_plot_deps)
adaswrite_SOURCES += $(adaswrite_deps)
efit2dg2d_SOURCES += $(efit2dg2d_deps)
geomtest_SOURCES += $(geomtest_deps)
geomtesta_SOURCES += $(geomtesta_deps)
matcheir_SOURCES += $(matcheir_deps)
acc_output_SOURCES += $(acc_output_deps)
reactionwrite_SOURCES += $(reactionwrite_deps)
pmiwrite_SOURCES += $(pmiwrite_deps)
problemsetup_SOURCES += $(problemsetup_deps)
geomread_SOURCES += $(geomread_deps)
cramdtest_SOURCES += $(cramdtest_deps)
cramdtesta_SOURCES += $(cramdtesta_deps)
degasgeom_SOURCES += $(degasgeom_deps)
test_SOURCES += $(test_deps)
readgeometry_SOURCES += $(readgeometry_deps)
definegeometry2d_SOURCES += $(definegeometry2d_deps)
defineback_SOURCES += $(defineback_deps)
boxgen_SOURCES += $(boxgen_deps)
plotflight_SOURCES += $(plotflight_deps)
cramdproblem_SOURCES += $(cramdproblem_deps)
readfilenames_SOURCES += $(readfilenames_deps)
readbackground_SOURCES += $(readbackground_deps)
updatebackground_SOURCES += $(updatebackground_deps)
reactiontest_SOURCES += $(reactiontest_deps)
pmitest_SOURCES += $(pmitest_deps)
sourcetest_SOURCES += $(sourcetest_deps)
voigttest_SOURCES += $(voigttest_deps)
voigtspec_SOURCES += $(voigtspec_deps)
ratecalc_SOURCES += $(ratecalc_deps)
dataexam_SOURCES += $(dataexam_deps)
outputbrowser_SOURCES += $(outputbrowser_deps)
postdetector_SOURCES += $(postdetector_deps)
datamatch_SOURCES += $(datamatch_deps)
tallysetup_SOURCES += $(tallysetup_deps)
allocate_SOURCES += $(allocate_deps)
array_SOURCES += $(array_deps)
teststring_SOURCES += $(teststring_deps)
elementlookup_SOURCES += $(elementlookup_deps)
speciesetup_SOURCES += $(speciesetup_deps)
netcdf_SOURCES += $(netcdf_deps)

LDADD=${LIBNETCDF_LIBS}

SUFFIXES = .web .hweb .mweb _mod.f

AM_FFLAGS=-fallow-argument-mismatch $(LDFLAGS)
AM_FFLAGS+=-ffixed-form

#FTANGLE=$(PWD)/../fweb/Web/ftangle
FTANGLE=ftangle
TANGLEFLAGS=-mINLINE -W[ -n/ -n')' -Tv -j -yn8000 -yb200000 -ybs15000 -mFORTRAN90 -mMACOSX -mGFORTRAN


#FC=gfortran
#AM_FFLAGS+=-fPIC -g
#LD=$(FC)

.PRECIOUS: %.f %.web.f

.web.f: 
	$(FTANGLE) $(TANGLEFLAGS) $<

.hweb.mweb: makemodule
%.mweb: %.hweb makemodule
	./makemodule $< $@

makemodule: makemodule.sh
	cp $< $@
	chmod a+x $@


#.mweb_mod.f: array.hweb sysdep.hweb
%_mod.f: %.mweb array.hweb sysdep.hweb
	$(FTANGLE) $(TANGLEFLAGS) -=n=#_mod.f $<

makedepend: makedepend.sh
	cp $< $@
	chmod a+x $@

parsedepends: parsedepends.sh
	cp $< $@
	chmod a+x $@

depend: makedepend parsedepends
	./makedepend *.web *.hweb
	./parsedepends

psiinterp.o: bspline.mod

#FORTRAN90=yes
#O=o
#-include Makefile.depends

#randomtest_SOURCES += $(randomtest_deps) $(foreach arg,$(subst .f,_mods,$(randomtest_deps)),$($(arg)))
	
#datasetup_SOURCES += $(datasetup_deps) $(foreach arg,$(subst .f,_mods,$(datasetup_deps)),$($(arg)))
#randomtest$(EXEEXT): sysdep.o random.o string.o

#randomtest.f: array.hweb constants.hweb macros.hweb random.hweb sysdep.hweb vector.hweb sysdep.o

#.f.o:
#	$(FC) -c $(AM_FFLAGS) $<

#%: %.o
#	$(LD) $^ -o $@

