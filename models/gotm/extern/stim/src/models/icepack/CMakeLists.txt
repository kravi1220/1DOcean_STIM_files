if(STIM_ICEPACK)
   # Icepack column physics + Icepack's own column driver (minus its main program),
   # built as a separate library with Icepack's compile-time sizes, then linked into STIM.
   # Sizes: one column, one thickness category (no ITD: no ridging/advection in 1D),
   # 7 ice layers, 1 snow layer; tracer switches as in Icepack's smoke test.
   file(GLOB ICEPACK_COLUMN ${CMAKE_CURRENT_LIST_DIR}/columnphysics/*.F90)
   file(GLOB ICEPACK_DRIVER ${CMAKE_CURRENT_LIST_DIR}/driver/*.F90)
   add_library(stim_icepack_lib STATIC ${ICEPACK_COLUMN} ${ICEPACK_DRIVER})
   target_compile_definitions(stim_icepack_lib PRIVATE
      FORTRANUNDERSCORE NO_SNICARHC
      NXGLOB=1 NICELYR=7 NSNWLYR=1 NICECAT=1 NFSDCAT=1
      TRAGE=1 TRFY=1 TRLVL=1 TRPND=1 TRSNOW=0 TRBRI=0 NTRISO=0 NTRAERO=1 NBGCLYR=1
      TRALG=0 TRBGCZ=0 TRDOC=0 TRDIC=0 TRDON=0 TRFED=0 TRFEP=0 TRZAERO=0 TRBGCS=0)
   target_compile_options(stim_icepack_lib PRIVATE -ffree-line-length-none)
   set_property(TARGET stim_icepack_lib PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/icepack_modules)
   target_include_directories(stim_icepack_lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/icepack_modules>)

   install(TARGETS stim_icepack_lib EXPORT stimConfig
           LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

   target_sources(stim PRIVATE ${CMAKE_CURRENT_LIST_DIR}/stim_icepack.F90)
   target_link_libraries(stim PUBLIC stim_icepack_lib)
   target_compile_definitions(stim PUBLIC STIM_ICEPACK)
endif(STIM_ICEPACK)
