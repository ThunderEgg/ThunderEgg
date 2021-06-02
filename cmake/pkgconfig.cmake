
if(TARGET FFTW::FFTW)
  string(APPEND pc_req_public " fftw3")
endif()
if(TARGET PETSC::PETSC)
  string(APPEND pc_req_public " PETSc")
endif()
if(TARGET P4EST::P4EST)
  string(APPEND pc_req_public " p4est libsc")
endif()
if(TARGET BLAS::BLAS)
  string(APPEND pc_libs_public " -lblas")
endif()
if(TARGET LAPACK::LAPACK)
  string(APPEND pc_libs_public " -llapack")
endif()

configure_file(cmake/ThunderEgg.pc.in ThunderEgg.pc @ONLY) 

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ThunderEgg.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)