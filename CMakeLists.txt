cmake_minimum_required(VERSION 3.12)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project(DEGAS2 Fortran C)
#enable_language (Fortran C)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#cmake_policy(SET CMP0074 NEW)

option(GRAPH_FILE "One of: SILO (default), HDF, or SILO_HDF5" SILO)
option(USE_MPI "Whether to use MPI when compiling" OFF)
option(SILO_ENABLE_HDF5 "Whether SILO is compiled with HDF5. Off by default because multiple versions on cluster can conflict" OFF)
option(NETCDF_DIR "Directory in which to look for NetCDF C and Fortran libraries")

########################################
# Dependencies
#
set (NETCDF_F90 "YES")
find_package (NetCDF REQUIRED)
if(USE_MPI)
   find_package (MPI REQUIRED)
endif()
find_package (HDF4)
find_package (HDF5)
find_package (SILO)


if (NETCDF_FOUND)
   message(STATUS "NetCDF found: ${NETCDF_LIBRARIES}")
endif (NETCDF_FOUND)

if (HDF4_FOUND)
   message(STATUS "HDF4 found: ${HDF4_INCLUDE_DIR}")
endif (HDF4_FOUND)

if (MPI_FOUND)
   message(STATUS "MPI found: ${MPI_Fortran_LIBRARIES}")
   if (USE_MPI)
      list(APPEND TANGLEFLAGS -mMPI)
      message(STATUS "Using MPI flag for FWEB")
   else()
      message(STATUS "Ignoring MPI")
   endif()

endif (MPI_FOUND)

message(STATUS ${CMAKE_Fortran_COMPILER_ID})
if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")

   if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
   #if (CMAKE_COMPILER_ID STREQUAL "GNU")
   #      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
      list(APPEND CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")
   #      add_compile_options(-fallow-argument-mismatch)
   endif()
endif()

########################################
# git submodules
#

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

add_subdirectory(deps/silo)
#add_subdirectory(deps/triangle)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/deps/silo/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

########################################
# FWEB Flags
#

if (CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT WIN32)
   list(APPEND TANGLEFLAGS -mDEBUG)
endif()

list(APPEND TANGLEFLAGS -mINLINE)
list(APPEND TANGLEFLAGS -mSILO)
list(APPEND TANGLEFLAGS -n/)
list(APPEND TANGLEFLAGS -Tv)
list(APPEND TANGLEFLAGS -j)
# If you want line information output as comments to .f files, remove the following flag
# Without this flag, be mindful that the absolute file path can't be too long or
# fweb will try continuing a comment to the next line, which is not allowed.
list(APPEND TANGLEFLAGS "-#")
list(APPEND TANGLEFLAGS -yn8000)
list(APPEND TANGLEFLAGS -yb200000)
list(APPEND TANGLEFLAGS -ybs15000)
#list(APPEND TANGLEFLAGS -mDEBUG)
list(APPEND TANGLEFLAGS -mSILO)
list(APPEND TANGLEFLAGS -mFORTRAN90)
list(APPEND TANGLEFLAGS -mLINUX64)
list(APPEND TANGLEFLAGS -W\[)

########################################
# Build Silo at configure time if not found
#
#if (SILO_FOUND)
#   message(STATUS "SILO found: ${SILO_LIBRARIES}")
#   message(STATUS "SILO include: ${SILO_INCLUDE_DIRS}")
#else()
#   message(FATAL_ERROR "SILO not found.")
#endif (SILO_FOUND)

########################################
# Build FWEB at configure time and set flags for ftangle
#
if(EXISTS "${PROJECT_SOURCE_DIR}/deps/fweb/Web/ftangle")
   set(FTANGLE ${PROJECT_SOURCE_DIR}/deps/fweb/Web/ftangle)
   message(STATUS "FWEB found: ${FTANGLE}")
else()
   message(STATUS "FWEB not found. Building from source.")
   if (NOT "${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
      message(FATAL_ERROR "Building FWEB requires GNU compiler. When configuring with CMake for the first time, specify path with -DCMAKE_C_COMPILER=[path to gcc]")
   endif()
   execute_process(
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/deps/fweb/Web
      COMMAND ./configure
   )
   execute_process(
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/deps/fweb/Web
      COMMAND make
   )
   set(FTANGLE ${PROJECT_SOURCE_DIR}/deps/fweb/Web/ftangle)
endif()


########################################
# General rule for building .f files from .web
#
file(GLOB_RECURSE WEBSOURCES ${PROJECT_SOURCE_DIR}/src/*.web)
set(ffiles)
foreach(webfile ${WEBSOURCES})
  string(REPLACE ".web" ".f" ffile ${webfile})
  get_filename_component(ffile ${ffile} NAME)
  add_custom_command(
    OUTPUT ${ffile}
    # Putting in this flag manually because I can't get it to play nice with the TANGLEFLAGS list
    COMMAND ${FTANGLE} ${TANGLEFLAGS} "-n')'" ${webfile}
    DEPENDS ${webfile}
    COMMAND_EXPAND_LISTS
  )
   add_custom_target(${ffile} DEPENDS ${webfile})
endforeach()

########################################
# Special treatment for definegeometry2d 
# Also compiles a .c file also to call triangle routines
#
file(GLOB dg2d_web ${PROJECT_SOURCE_DIR}/src/definegeometry2d.web)
string(REPLACE ".web" "c.c" dg2d_c ${dg2d_web})
get_filename_component(dg2d_c ${dg2d_c} NAME)
#list(APPEND TANGLEFLAGS -=c=\#c\.c)
add_custom_command(
   OUTPUT ${dg2d_c}
   COMMAND ${FTANGLE} ${TANGLEFLAGS} -=c=\#c\.c "-n')'" ${dg2d_web}
   DEPENDS ${dg2d_web}
   COMMAND_EXPAND_LISTS
)
add_custom_target(dg2d_c DEPENDS ${dg2d_web})
add_library(dg2dc OBJECT ${dg2d_c})
target_include_directories(dg2dc PRIVATE ${PROJECT_SOURCE_DIR}/deps/triangle)

###################################
# Look for "usr" files that specify 
# non-default detector and plasma properties
#

if(EXISTS ${PROJECT_SOURCE_DIR}/src/usr2dplasma.web)
   set(PLASMAFILE usr2dplasma)
else()
   set(PLASMAFILE def2dplasma)
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/src/usr2ddetector.web)
   set(DETECTORFILE usr2ddetector)
else()
   set(DETECTORFILE def2ddetector)
endif()

########################################
# Copy template files to data directory
# Replace D2ROOT in those files with this
# source directory.
#
file(READ ${CMAKE_SOURCE_DIR}/data/templates/pmi.input filedata)
string(REPLACE D2ROOT ${CMAKE_SOURCE_DIR} filedata ${filedata})
file(WRITE ${CMAKE_SOURCE_DIR}/data/pmi.input ${filedata})

file(READ ${CMAKE_SOURCE_DIR}/data/templates/reactions.input filedata)
string(REPLACE D2ROOT ${CMAKE_SOURCE_DIR} filedata ${filedata})
file(WRITE ${CMAKE_SOURCE_DIR}/data/reactions.input ${filedata})

file(READ ${CMAKE_SOURCE_DIR}/data/templates/degas2.in filedata)
string(REPLACE D2ROOT ${CMAKE_SOURCE_DIR} filedata ${filedata})
file(WRITE data/degas2.in ${filedata})

file(READ ${CMAKE_SOURCE_DIR}/data/templates/tally.in filedata)
file(WRITE data/tally.in ${filedata})

########################################
# Specify "modules"
#
set(MODULES allocate)
list(APPEND MODULES background)
list(APPEND MODULES degas2_xgc)
list(APPEND MODULES detector)
list(APPEND MODULES efititp)
list(APPEND MODULES element)
list(APPEND MODULES flight_frag)
list(APPEND MODULES geometry2d)
list(APPEND MODULES geomint)
list(APPEND MODULES macros)
list(APPEND MODULES materials)
list(APPEND MODULES mpi)
list(APPEND MODULES output)
list(APPEND MODULES pmi)
list(APPEND MODULES pmidata)
list(APPEND MODULES pmiformat)
list(APPEND MODULES postdetector)
list(APPEND MODULES problem)
list(APPEND MODULES ratecalc)
list(APPEND MODULES reaction)
list(APPEND MODULES reactiondata)
list(APPEND MODULES readfilenames)
list(APPEND MODULES sector)
list(APPEND MODULES snapshot_pdf)
list(APPEND MODULES sources)
list(APPEND MODULES species)
list(APPEND MODULES stat)
list(APPEND MODULES tally)
list(APPEND MODULES xgc1interp)
list(APPEND MODULES xsection)
list(APPEND MODULES zone)

foreach(module ${MODULES})

  set(hfile "${PROJECT_SOURCE_DIR}/src/${module}.hweb" )
  set(mfile "${module}.mweb" )
  set(ffile "${module}_mod.f" )
  add_custom_command(
    OUTPUT ${mfile}
    COMMAND ${PROJECT_SOURCE_DIR}/src/makemodule.sh ${hfile} ${mfile}
    DEPENDS ${hfile}
  )
  add_custom_command(
    OUTPUT ${ffile}
    COMMAND ${FTANGLE} ${TANGLEFLAGS} "-n')'" "-=n=#_mod.f" -I${PROJECT_SOURCE_DIR}/src ${mfile}
    DEPENDS ${mfile}
    COMMAND_EXPAND_LISTS
  )
  add_custom_target(${ffile} DEPENDS ${mfile})
endforeach()

#################################
# Special dependencies
# 
add_library(bspline90_22 OBJECT ${PROJECT_SOURCE_DIR}/src/EFIT2dg2d/bspline90_22.f90)

add_library(triangle OBJECT ${PROJECT_SOURCE_DIR}/deps/triangle/triangle.c)
target_compile_options(triangle PRIVATE -DTRILIBRARY -O)
#add_executable(randomtest randomtest.f random.f sysdep.f string.f)

add_library(aladdin OBJECT ${PROJECT_SOURCE_DIR}/data/Aladdin/aladdin_subr.f)
target_compile_options(aladdin PRIVATE -ffixed-form)

list(APPEND efit2dg2d_deps $<TARGET_OBJECTS:bspline90_22>)
list(APPEND defineback_deps $<TARGET_OBJECTS:bspline90_22>)
list(APPEND definegeometry2d_deps $<TARGET_OBJECTS:triangle> $<TARGET_OBJECTS:dg2dc>)
list(APPEND ratecalc_deps $<TARGET_OBJECTS:aladdin>)

########################################
# Until the line #<<<<<<<<<<<<<, the following is
# automatically generated from "make depend"
# which creates cmakecommands.txt.
# Paste the contents of that file here
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

set(al_MOD allocate_mod.f)
set(bk_MOD background_mod.f)
set(dx_MOD degas2_xgc_mod.f)
set(de_MOD detector_mod.f)
set(ef_MOD efititp_mod.f)
set(el_MOD element_mod.f)
set(ff_MOD flight_frag_mod.f)
set(g2_MOD geometry2d_mod.f)
set(gi_MOD geomint_mod.f)
set(cm_MOD macros_mod.f)
set(ma_MOD materials_mod.f)
set(mp_MOD mpi_mod.f)
set(ou_MOD output_mod.f)
set(pm_MOD pmi_mod.f)
set(pd_MOD pmidata_mod.f)
set(pf_MOD pmiformat_mod.f)
set(po_MOD postdetector_mod.f)
set(pr_MOD problem_mod.f)
set(ps_MOD problem_mod.f)
set(ra_MOD ratecalc_mod.f)
set(rc_MOD reaction_mod.f)
set(rd_MOD reactiondata_mod.f)
set(rf_MOD readfilenames_mod.f)
set(sc_MOD sector_mod.f)
set(sn_MOD snapshot_pdf_mod.f)
set(so_MOD sources_mod.f)
set(sp_MOD species_mod.f)
set(sa_MOD stat_mod.f)
set(tl_MOD tally_mod.f)
set(xi_MOD xgc1interp_mod.f)
set(xs_MOD xsection_mod.f)
set(zn_MOD zone_mod.f)
list(APPEND acc_output_mods  ${ou_MOD} ${pr_MOD} ${rf_MOD} ${so_MOD} ${tl_MOD})
list(APPEND adaswrite_mods  ${rc_MOD} ${sp_MOD} ${xs_MOD})
list(APPEND allocate_mods  ${al_MOD})
list(APPEND aux2dsetup_mods  ${gi_MOD} ${zn_MOD})
list(APPEND boxgen_mods  ${bk_MOD} ${de_MOD} ${gi_MOD} ${ma_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND btopdetector_mods  ${de_MOD} ${zn_MOD})
list(APPEND collide_mods  ${bk_MOD} ${gi_MOD} ${pr_MOD} ${rc_MOD} ${rd_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND composite_mods  ${gi_MOD} ${zn_MOD})
list(APPEND conversions_mods  ${bk_MOD} ${gi_MOD} ${pr_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND cramdadd_mods  ${bk_MOD} ${pr_MOD} ${rc_MOD} ${rd_MOD} ${rf_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND dataexam_mods  ${rc_MOD} ${xs_MOD})
list(APPEND datamatch_mods  ${xs_MOD})
list(APPEND dediagsetup_mods  ${de_MOD} ${gi_MOD} ${pr_MOD} ${sc_MOD} ${zn_MOD})
list(APPEND def2ddetector_mods  ${de_MOD} ${zn_MOD})
list(APPEND def2dplasma_mods  ${bk_MOD} ${gi_MOD} ${zn_MOD})
list(APPEND defineback_mods  ${bk_MOD} ${gi_MOD} ${pr_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND definegeometry2d_mods  ${g2_MOD} ${gi_MOD} ${ma_MOD} ${sc_MOD} ${zn_MOD})
list(APPEND degas2_xgc_mods  ${bk_MOD} ${dx_MOD} ${mp_MOD} ${ou_MOD} ${pr_MOD} ${rc_MOD} ${rd_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND degasinit_mods  ${mp_MOD} ${pr_MOD} ${tl_MOD})
list(APPEND doflights_mods  ${bk_MOD} ${de_MOD} ${ff_MOD} ${mp_MOD} ${ou_MOD} ${pr_MOD} ${rf_MOD} ${sa_MOD} ${sn_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND efit2dg2d_mods  ${ef_MOD})
list(APPEND elementread_mods  ${bk_MOD} ${el_MOD} ${mp_MOD} ${rc_MOD} ${rf_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND elementsetup_mods  ${el_MOD} ${rf_MOD})
list(APPEND flight_mods  ${bk_MOD} ${ff_MOD} ${pm_MOD} ${pr_MOD} ${rc_MOD} ${sa_MOD} ${sc_MOD} ${sn_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND flighttest_mods  ${bk_MOD} ${ff_MOD} ${mp_MOD} ${ou_MOD} ${pr_MOD} ${rc_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND geometry_mods  ${de_MOD} ${gi_MOD} ${mp_MOD} ${rf_MOD} ${sc_MOD} ${zn_MOD})
list(APPEND geomfront_mods  ${gi_MOD} ${sc_MOD})
list(APPEND geomsetup_mods  ${al_MOD} ${de_MOD} ${gi_MOD} ${rf_MOD} ${sc_MOD} ${zn_MOD})
list(APPEND geomtest_mods  ${bk_MOD} ${zn_MOD})
list(APPEND geomtesta_mods  ${bk_MOD} ${de_MOD} ${ou_MOD} ${pr_MOD} ${rf_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND gpicamera_mods  ${de_MOD} ${zn_MOD})
list(APPEND ltxlymandetector_mods  ${de_MOD} ${zn_MOD})
list(APPEND matcheir_mods  ${ou_MOD} ${pr_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND matchout_mods  ${ou_MOD} ${pr_MOD} ${rf_MOD} ${so_MOD} ${tl_MOD})
list(APPEND materialsetup_mods  ${ma_MOD} ${rf_MOD})
list(APPEND ncsxplasma_mods  ${bk_MOD} ${g2_MOD} ${pr_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND nstxgpi2_mods  ${bk_MOD} ${gi_MOD} ${zn_MOD})
list(APPEND outputbrowser_mods  ${de_MOD} ${gi_MOD} ${ou_MOD} ${pm_MOD} ${pr_MOD} ${rc_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND plate_mods  ${pd_MOD} ${pf_MOD} ${pm_MOD} ${pr_MOD} ${sp_MOD})
list(APPEND plotgeom_mods  ${gi_MOD} ${zn_MOD})
list(APPEND pmimatread_mods  ${ma_MOD} ${mp_MOD} ${ou_MOD} ${pd_MOD} ${pm_MOD} ${pr_MOD} ${rd_MOD} ${rf_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD})
list(APPEND pmisetup_mods  ${ma_MOD} ${pm_MOD} ${rf_MOD} ${sp_MOD})
list(APPEND pmitest_mods  ${ma_MOD} ${pd_MOD} ${pm_MOD} ${pr_MOD} ${rc_MOD} ${sp_MOD})
list(APPEND pmiwrite_mods  ${el_MOD} ${pf_MOD})
list(APPEND poly_to_tri_mods  ${g2_MOD} ${sc_MOD} ${zn_MOD})
list(APPEND tri_to_sonnet_mods  ${g2_MOD} )
list(APPEND postdetector_mods  ${bk_MOD} ${de_MOD} ${mp_MOD} ${ou_MOD} ${po_MOD} ${pr_MOD} ${rf_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND problemsetup_mods  ${el_MOD} ${ma_MOD} ${mp_MOD} ${pm_MOD} ${pr_MOD} ${ps_MOD} ${rc_MOD} ${rf_MOD} ${sp_MOD})
list(APPEND psiinterp_mods  ${ef_MOD})
list(APPEND pvmmaster_mods  ${zn_MOD})
list(APPEND pvmslave_mods  ${zn_MOD})
list(APPEND ratecalc_mods  ${ra_MOD} ${rf_MOD} ${xs_MOD})
list(APPEND reaction_mods  ${bk_MOD} ${gi_MOD} ${pr_MOD} ${rc_MOD} ${rd_MOD} ${sp_MOD} ${xs_MOD} ${zn_MOD})
list(APPEND reactionsetup_mods  ${rc_MOD} ${rf_MOD} ${sp_MOD})
list(APPEND reactiontest_mods  ${el_MOD} ${pr_MOD} ${rc_MOD} ${rd_MOD} ${sp_MOD})
list(APPEND reactionwrite_mods  ${xs_MOD})
list(APPEND readbackground_mods  ${bk_MOD} ${el_MOD} ${gi_MOD} ${pr_MOD} ${rc_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND readfilenames_mods  ${rf_MOD})
list(APPEND readgeometry_mods  ${de_MOD} ${gi_MOD} ${ma_MOD} ${sc_MOD} ${zn_MOD})
list(APPEND score_mods  ${bk_MOD} ${de_MOD} ${gi_MOD} ${pr_MOD} ${rc_MOD} ${rd_MOD} ${sa_MOD} ${sc_MOD} ${sn_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND sources_mods  ${bk_MOD} ${gi_MOD} ${pr_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND sourcetest_mods  ${gi_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND speciesetup_mods  ${el_MOD} ${rf_MOD} ${sp_MOD})
list(APPEND stat_mods  ${ou_MOD} ${pr_MOD} ${sa_MOD} ${tl_MOD})
list(APPEND sysdep_mods  ${mem_MOD})
list(APPEND tallysetup_mods  ${de_MOD} ${gi_MOD} ${pr_MOD} ${rc_MOD} ${rf_MOD} ${sc_MOD} ${so_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND ucd_plot_mods  ${bk_MOD} ${de_MOD} ${g2_MOD} ${ou_MOD} ${pr_MOD} ${rf_MOD} ${sc_MOD} ${sp_MOD} ${tl_MOD} ${zn_MOD})
list(APPEND uedgeback_mods  ${so_MOD})
list(APPEND updatebackground_mods  ${bk_MOD} ${pr_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND writeback_mods  ${bk_MOD} ${el_MOD} ${g2_MOD} ${gi_MOD} ${mp_MOD} ${pr_MOD} ${rc_MOD} ${rf_MOD} ${sc_MOD} ${sn_MOD} ${so_MOD} ${sp_MOD} ${zn_MOD})
list(APPEND xgc1_gpi_mods  ${bk_MOD} ${gi_MOD} ${xi_MOD} ${zn_MOD})
list(APPEND sysdeptest_deps sysdeptest.f)
list(APPEND randomtest_deps randomtest.f random.f)
list(APPEND pvmmaster_deps pvmmaster.f random.f stat.f geometry.f readfilenames.f)
list(APPEND pvmslave_deps pvmslave.f random.f randother.f stat.f geometry.f elementread.f pmimatread.f flight.f reaction.f interpolate.f geomfront.f readfilenames.f)
list(APPEND flighttest_deps flighttest.f doflights.f elementread.f pmimatread.f geometry.f random.f stat.f randother.f reaction.f interpolate.f geomfront.f flight.f sources.f collide.f readfilenames.f plate.f degasinit.f score.f conversions.f voigt.f writeback.f aux2dsetup.f geomsetup.f geomfaces.f geomutils.f allocate.f composite.f)
list(APPEND degas2_xgc_deps degas2_xgc.f doflights.f elementread.f pmimatread.f geometry.f random.f stat.f randother.f reaction.f interpolate.f geomfront.f flight.f sources.f collide.f readfilenames.f plate.f degasinit.f score.f conversions.f voigt.f writeback.f aux2dsetup.f geomsetup.f geomfaces.f geomutils.f allocate.f composite.f)
list(APPEND poly_to_tri_deps poly_to_tri.f geometry.f allocate.f readfilenames.f )
list(APPEND ucd_plot_deps ucd_plot.f geometry.f allocate.f readfilenames.f elementread.f pmimatread.f  reaction.f interpolate.f random.f degasinit.f)
list(APPEND adaswrite_deps adaswrite.f adas_xxdata_11.f elementread.f reaction.f interpolate.f readfilenames.f random.f)
list(APPEND efit2dg2d_deps efit2dg2d.f aux2dsetup.f psiinterp.f interpolate.f)
list(APPEND geomtest_deps geomtest.f geometry.f geomsetup.f allocate.f geomutils.f geomfaces.f composite.f string.f readfilenames.f geomfront.f dediagsetup.f)
list(APPEND geomtesta_deps geomtesta.f geometry.f plotgeom.f elementread.f pmimatread.f readfilenames.f reaction.f interpolate.f collide.f randother.f random.f conversions.f voigt.f degasinit.f )
list(APPEND matchout_deps matchout.f pmimatread.f elementread.f geometry.f readfilenames.f degasinit.f)
list(APPEND matcheir_deps matcheir.f pmimatread.f elementread.f geometry.f readfilenames.f degasinit.f)
list(APPEND acc_output_deps acc_output.f pmimatread.f elementread.f geometry.f readfilenames.f degasinit.f conversions.f)
list(APPEND reactionwrite_deps reactionwrite.f)
list(APPEND pmiwrite_deps pmiwrite.f interpolate.f random.f readfilenames.f elementread.f )
list(APPEND datasetup_deps datasetup.f elementread.f pmimatread.f elementsetup.f speciesetup.f reactionsetup.f materialsetup.f pmisetup.f readfilenames.f)
list(APPEND problemsetup_deps problemsetup.f elementread.f pmimatread.f reaction.f interpolate.f readfilenames.f plate.f geometry.f random.f randother.f)
list(APPEND cramdtesta_deps cramdtesta.f cramdadd.f sysdep.f string.f pmimatread.f)
list(APPEND readgeometry_deps readgeometry.f geomsetup.f allocate.f geomutils.f geomfaces.f geometry.f composite.f string.f readdegas.f readfilenames.f elementread.f pmimatread.f dediagsetup.f aux2dsetup.f ${DETECTORFILE}.f)
list(APPEND definegeometry2d_deps  string.f geomsetup.f allocate.f composite.f geomutils.f geomfaces.f geometry.f aux2dsetup.f pmimatread.f dediagsetup.f readfilenames.f plotgeom.f ${DETECTORFILE}.f)
list(APPEND defineback_deps defineback.f string.f writeback.f pmimatread.f geometry.f elementread.f geomsetup.f aux2dsetup.f allocate.f geomutils.f geomfaces.f readfilenames.f reaction.f interpolate.f random.f composite.f ${PLASMAFILE}.f psiinterp.f )
list(APPEND boxgen_deps boxgen.f readfilenames.f geomsetup.f allocate.f composite.f geomutils.f geomfaces.f geometry.f pmimatread.f elementread.f writeback.f cramdadd.f dediagsetup.f reaction.f interpolate.f random.f aux2dsetup.f)
list(APPEND cramdproblem_deps cramdproblem.f pmimatread.f cramdadd.f readfilenames.f elementread.f geometry.f)
list(APPEND readbackground_deps readbackground.f readdegas.f geometry.f elementread.f readfilenames.f geomsetup.f allocate.f reaction.f interpolate.f composite.f geomutils.f geomfaces.f geometry.f writeback.f pmimatread.f cramdadd.f random.f aux2dsetup.f)
list(APPEND updatebackground_deps updatebackground.f elementread.f readfilenames.f geomsetup.f allocate.f writeback.f pmimatread.f random.f geomutils.f geometry.f geomfaces.f reaction.f interpolate.f)
list(APPEND reactiontest_deps reactiontest.f elementread.f reaction.f interpolate.f readfilenames.f random.f randother.f collide.f pmimatread.f geometry.f)
list(APPEND pmitest_deps pmitest.f elementread.f pmimatread.f plate.f readfilenames.f random.f randother.f geometry.f interpolate.f)
list(APPEND sourcetest_deps sourcetest.f sources.f readfilenames.f randother.f random.f pmimatread.f elementread.f geometry.f cramdadd.f interpolate.f score.f collide.f reaction.f conversions.f)
list(APPEND voigttest_deps voigttest.f voigt.f)
list(APPEND voigtspec_deps voigtspec.f voigt.f)
list(APPEND ratecalc_deps ratecalc.f readfilenames.f interpolate.f)
list(APPEND dataexam_deps dataexam.f readfilenames.f elementread.f)
list(APPEND outputbrowser_deps outputbrowser.f readfilenames.f degasinit.f pmimatread.f elementread.f geometry.f)
list(APPEND postdetector_deps postdetector.f readfilenames.f degasinit.f pmimatread.f elementread.f geometry.f dediagsetup.f ${DETECTORFILE}.f)
list(APPEND datamatch_deps datamatch.f )
list(APPEND tallysetup_deps tallysetup.f elementread.f pmimatread.f readfilenames.f string.f geometry.f)
list(APPEND allocate_deps allocate.f)
foreach(str ${sysdeptest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND sysdeptest_mods "${${modstr}}")
endforeach()
add_executable(sysdeptest sysdeptest.f string.f sysdep.f ${sysdeptest_mods} ${sysdeptest_deps})
target_link_libraries(sysdeptest ${NETCDF_LIBRARIES})
 
foreach(str ${randomtest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND randomtest_mods "${${modstr}}")
endforeach()
add_executable(randomtest randomtest.f string.f sysdep.f ${randomtest_mods} ${randomtest_deps})
target_link_libraries(randomtest ${NETCDF_LIBRARIES})
 
foreach(str ${pvmmaster_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND pvmmaster_mods "${${modstr}}")
endforeach()
add_executable(pvmmaster pvmmaster.f string.f sysdep.f ${pvmmaster_mods} ${pvmmaster_deps})
target_link_libraries(pvmmaster ${NETCDF_LIBRARIES})
 
foreach(str ${pvmslave_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND pvmslave_mods "${${modstr}}")
endforeach()
add_executable(pvmslave pvmslave.f string.f sysdep.f ${pvmslave_mods} ${pvmslave_deps})
target_link_libraries(pvmslave ${NETCDF_LIBRARIES})
 
foreach(str ${flighttest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND flighttest_mods "${${modstr}}")
endforeach()
add_executable(flighttest flighttest.f string.f sysdep.f ${flighttest_mods} ${flighttest_deps})
target_link_libraries(flighttest ${NETCDF_LIBRARIES})
 
foreach(str ${degas2_xgc_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND degas2_xgc_mods "${${modstr}}")
endforeach()
#add_executable(degas2_xgc degas2_xgc.f string.f sysdep.f ${degas2_xgc_mods} ${degas2_xgc_deps})
#target_link_libraries(degas2_xgc ${NETCDF_LIBRARIES})
 
foreach(str ${poly_to_tri_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND poly_to_tri_mods "${${modstr}}")
endforeach()
add_executable(poly_to_tri poly_to_tri.f string.f sysdep.f ${poly_to_tri_mods} ${poly_to_tri_deps})
target_link_libraries(poly_to_tri ${NETCDF_LIBRARIES})

foreach(str ${tri_to_sonnet_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND tri_to_sonnet_mods "${${modstr}}")
endforeach()
add_executable(tri_to_sonnet tri_to_sonnet.f string.f sysdep.f ${tri_to_sonnet_mods} ${tri_to_sonnet_deps})

target_link_libraries(poly_to_tri ${NETCDF_LIBRARIES})
 
foreach(str ${ucd_plot_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND ucd_plot_mods "${${modstr}}")
endforeach()
add_executable(ucd_plot ucd_plot.f string.f sysdep.f ${ucd_plot_mods} ${ucd_plot_deps})
target_link_libraries(ucd_plot ${NETCDF_LIBRARIES})
 
foreach(str ${adaswrite_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND adaswrite_mods "${${modstr}}")
endforeach()
add_executable(adaswrite adaswrite.f string.f sysdep.f ${adaswrite_mods} ${adaswrite_deps})
target_link_libraries(adaswrite ${NETCDF_LIBRARIES})
 
foreach(str ${efit2dg2d_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND efit2dg2d_mods "${${modstr}}")
endforeach()
add_executable(efit2dg2d efit2dg2d.f string.f sysdep.f ${efit2dg2d_mods} ${efit2dg2d_deps})
target_link_libraries(efit2dg2d ${NETCDF_LIBRARIES})
 
foreach(str ${geomtest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND geomtest_mods "${${modstr}}")
endforeach()
add_executable(geomtest geomtest.f string.f sysdep.f ${geomtest_mods} ${geomtest_deps})
target_link_libraries(geomtest ${NETCDF_LIBRARIES})
 
foreach(str ${geomtesta_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND geomtesta_mods "${${modstr}}")
endforeach()
add_executable(geomtesta geomtesta.f string.f sysdep.f ${geomtesta_mods} ${geomtesta_deps})
target_link_libraries(geomtesta ${NETCDF_LIBRARIES})
 
foreach(str ${matchout_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND matchout_mods "${${modstr}}")
endforeach()
add_executable(matchout matchout.f string.f sysdep.f ${matchout_mods} ${matchout_deps})
target_link_libraries(matchout ${NETCDF_LIBRARIES})
 
foreach(str ${matcheir_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND matcheir_mods "${${modstr}}")
endforeach()
add_executable(matcheir matcheir.f string.f sysdep.f ${matcheir_mods} ${matcheir_deps})
target_link_libraries(matcheir ${NETCDF_LIBRARIES})
 
foreach(str ${acc_output_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND acc_output_mods "${${modstr}}")
endforeach()
add_executable(acc_output acc_output.f string.f sysdep.f ${acc_output_mods} ${acc_output_deps})
target_link_libraries(acc_output ${NETCDF_LIBRARIES})
 
foreach(str ${reactionwrite_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND reactionwrite_mods "${${modstr}}")
endforeach()
add_executable(reactionwrite reactionwrite.f string.f sysdep.f ${reactionwrite_mods} ${reactionwrite_deps})
target_link_libraries(reactionwrite ${NETCDF_LIBRARIES})
 
foreach(str ${pmiwrite_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND pmiwrite_mods "${${modstr}}")
endforeach()
add_executable(pmiwrite pmiwrite.f string.f sysdep.f ${pmiwrite_mods} ${pmiwrite_deps})
target_link_libraries(pmiwrite ${NETCDF_LIBRARIES})
 
foreach(str ${datasetup_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND datasetup_mods "${${modstr}}")
endforeach()
add_executable(datasetup datasetup.f string.f sysdep.f ${datasetup_mods} ${datasetup_deps})
target_link_libraries(datasetup ${NETCDF_LIBRARIES})
 
foreach(str ${problemsetup_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND problemsetup_mods "${${modstr}}")
endforeach()
add_executable(problemsetup problemsetup.f string.f sysdep.f ${problemsetup_mods} ${problemsetup_deps})
target_link_libraries(problemsetup ${NETCDF_LIBRARIES})
 
foreach(str ${cramdtesta_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND cramdtesta_mods "${${modstr}}")
endforeach()
add_executable(cramdtesta cramdtesta.f string.f sysdep.f ${cramdtesta_mods} ${cramdtesta_deps})
target_link_libraries(cramdtesta ${NETCDF_LIBRARIES})
 
foreach(str ${readgeometry_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND readgeometry_mods "${${modstr}}")
endforeach()
add_executable(readgeometry readgeometry.f string.f sysdep.f ${readgeometry_mods} ${readgeometry_deps})
target_link_libraries(readgeometry ${NETCDF_LIBRARIES})
 
foreach(str ${definegeometry2d_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND definegeometry2d_mods "${${modstr}}")
endforeach()
add_executable(definegeometry2d definegeometry2d.f string.f sysdep.f ${definegeometry2d_mods} ${definegeometry2d_deps})
target_link_libraries(definegeometry2d ${NETCDF_LIBRARIES})
 
foreach(str ${defineback_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND defineback_mods "${${modstr}}")
endforeach()
add_executable(defineback defineback.f string.f sysdep.f ${defineback_mods} ${defineback_deps})
target_link_libraries(defineback ${NETCDF_LIBRARIES})
 
foreach(str ${boxgen_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND boxgen_mods "${${modstr}}")
endforeach()
add_executable(boxgen boxgen.f string.f sysdep.f ${boxgen_mods} ${boxgen_deps})
target_link_libraries(boxgen ${NETCDF_LIBRARIES})
 
foreach(str ${cramdproblem_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND cramdproblem_mods "${${modstr}}")
endforeach()
add_executable(cramdproblem cramdproblem.f string.f sysdep.f ${cramdproblem_mods} ${cramdproblem_deps})
target_link_libraries(cramdproblem ${NETCDF_LIBRARIES})
 
foreach(str ${readbackground_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND readbackground_mods "${${modstr}}")
endforeach()
add_executable(readbackground readbackground.f string.f sysdep.f ${readbackground_mods} ${readbackground_deps})
target_link_libraries(readbackground ${NETCDF_LIBRARIES})
 
foreach(str ${updatebackground_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND updatebackground_mods "${${modstr}}")
endforeach()
add_executable(updatebackground updatebackground.f string.f sysdep.f ${updatebackground_mods} ${updatebackground_deps})
target_link_libraries(updatebackground ${NETCDF_LIBRARIES})
 
foreach(str ${reactiontest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND reactiontest_mods "${${modstr}}")
endforeach()
add_executable(reactiontest reactiontest.f string.f sysdep.f ${reactiontest_mods} ${reactiontest_deps})
target_link_libraries(reactiontest ${NETCDF_LIBRARIES})
 
foreach(str ${pmitest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND pmitest_mods "${${modstr}}")
endforeach()
add_executable(pmitest pmitest.f string.f sysdep.f ${pmitest_mods} ${pmitest_deps})
target_link_libraries(pmitest ${NETCDF_LIBRARIES})
 
foreach(str ${sourcetest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND sourcetest_mods "${${modstr}}")
endforeach()
add_executable(sourcetest sourcetest.f string.f sysdep.f ${sourcetest_mods} ${sourcetest_deps})
target_link_libraries(sourcetest ${NETCDF_LIBRARIES})
 
foreach(str ${voigttest_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND voigttest_mods "${${modstr}}")
endforeach()
add_executable(voigttest voigttest.f string.f sysdep.f ${voigttest_mods} ${voigttest_deps})
target_link_libraries(voigttest ${NETCDF_LIBRARIES})
 
foreach(str ${voigtspec_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND voigtspec_mods "${${modstr}}")
endforeach()
add_executable(voigtspec voigtspec.f string.f sysdep.f ${voigtspec_mods} ${voigtspec_deps})
target_link_libraries(voigtspec ${NETCDF_LIBRARIES})
 
foreach(str ${ratecalc_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND ratecalc_mods "${${modstr}}")
endforeach()
add_executable(ratecalc ratecalc.f string.f sysdep.f ${ratecalc_mods} ${ratecalc_deps})
target_link_libraries(ratecalc ${NETCDF_LIBRARIES})
 
foreach(str ${dataexam_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND dataexam_mods "${${modstr}}")
endforeach()
add_executable(dataexam dataexam.f string.f sysdep.f ${dataexam_mods} ${dataexam_deps})
target_link_libraries(dataexam ${NETCDF_LIBRARIES})
 
foreach(str ${outputbrowser_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND outputbrowser_mods "${${modstr}}")
endforeach()
add_executable(outputbrowser outputbrowser.f string.f sysdep.f ${outputbrowser_mods} ${outputbrowser_deps})
target_link_libraries(outputbrowser ${NETCDF_LIBRARIES})
 
foreach(str ${postdetector_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND postdetector_mods "${${modstr}}")
endforeach()
add_executable(postdetector postdetector.f string.f sysdep.f ${postdetector_mods} ${postdetector_deps})
target_link_libraries(postdetector ${NETCDF_LIBRARIES})
 
foreach(str ${datamatch_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND datamatch_mods "${${modstr}}")
endforeach()
add_executable(datamatch datamatch.f string.f sysdep.f ${datamatch_mods} ${datamatch_deps})
target_link_libraries(datamatch ${NETCDF_LIBRARIES})
 
foreach(str ${tallysetup_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND tallysetup_mods "${${modstr}}")
endforeach()
add_executable(tallysetup tallysetup.f string.f sysdep.f ${tallysetup_mods} ${tallysetup_deps})
target_link_libraries(tallysetup ${NETCDF_LIBRARIES})
 
foreach(str ${allocate_deps})
  string(REPLACE .f _mods modstr "${str}")
  list(APPEND allocate_mods "${${modstr}}")
endforeach()
add_executable(allocate allocate.f string.f sysdep.f ${allocate_mods} ${allocate_deps})
target_link_libraries(allocate ${NETCDF_LIBRARIES})
 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# End automatically generated dependencies
#################################################

target_link_libraries(definegeometry2d ${NETCDF_LIBRARIES} silo)
target_include_directories(definegeometry2d PUBLIC deps/silo/src/silo PUBLIC deps/triangle)


target_link_libraries(ucd_plot ${NETCDF_LIBRARIES} silo)
target_include_directories(ucd_plot PUBLIC deps/silo/src/silo)

separate_arguments(MPI_Fortran_COMPILE_FLAGS)

if (USE_MPI)
   target_link_libraries(flighttest ${MPI_Fortran_LIBRARIES})
   target_include_directories(flighttest PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(flighttest PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   
   target_link_libraries(datasetup ${MPI_Fortran_LIBRARIES})
   target_include_directories(datasetup PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(datasetup PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(problemsetup ${MPI_Fortran_LIBRARIES})
   target_include_directories(problemsetup PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(problemsetup PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(definegeometry2d ${MPI_Fortran_LIBRARIES})
   target_include_directories(definegeometry2d PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(definegeometry2d PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(defineback ${MPI_Fortran_LIBRARIES})
   target_include_directories(defineback PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(defineback PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(tallysetup ${MPI_Fortran_LIBRARIES})
   target_include_directories(tallysetup PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(tallysetup PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(ucd_plot ${MPI_Fortran_LIBRARIES})
   target_include_directories(ucd_plot PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_link_libraries(poly_to_tri ${MPI_Fortran_LIBRARIES})
   target_include_directories(poly_to_tri PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(ucd_plot PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(adaswrite ${MPI_Fortran_LIBRARIES})
   target_include_directories(adaswrite PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(adaswrite PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(dataexam ${MPI_Fortran_LIBRARIES})
   target_include_directories(dataexam PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(dataexam PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(outputbrowser ${MPI_Fortran_LIBRARIES})
   target_include_directories(outputbrowser PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(outputbrowser PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(postdetector ${MPI_Fortran_LIBRARIES})
   target_include_directories(postdetector PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(postdetector PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(boxgen ${MPI_Fortran_LIBRARIES})
   target_include_directories(boxgen PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(boxgen PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(readbackground ${MPI_Fortran_LIBRARIES})
   target_include_directories(readbackground PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(readbackground PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(readgeometry ${MPI_Fortran_LIBRARIES})
   target_include_directories(readgeometry PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(readgeometry PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(pmitest ${MPI_Fortran_LIBRARIES})
   target_include_directories(pmitest PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(pmitest PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
   target_link_libraries(sourcetest ${MPI_Fortran_LIBRARIES})
   target_include_directories(sourcetest PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
   target_compile_options(sourcetest PUBLIC ${MPI_Fortran_COMPILE_FLAGS})
endif()

#add_library(degas2_xgc STATIC degas2_xgc.f string.f sysdep.f ${degas2_xgc_mods} ${degas2_xgc_deps})
add_library(degas2_xgc SHARED degas2_xgc.f string.f sysdep.f ${degas2_xgc_mods} ${degas2_xgc_deps})
target_link_libraries(degas2_xgc ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES})
#target_link_libraries(degas2_xgc.a ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES})
#target_include_directories(degas2_xgc.a PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
target_include_directories(degas2_xgc PUBLIC ${MPI_Fortran_INCLUDE_DIRS})
