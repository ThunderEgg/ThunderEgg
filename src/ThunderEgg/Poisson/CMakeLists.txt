# -- clear header list form parent scope

unset(ThunderEgg_HDRS)

# -- add sources

list(APPEND ThunderEgg_HDRS StarPatchOperator.h)
target_sources(ThunderEgg PRIVATE StarPatchOperator.cpp)

if(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)
  list(APPEND ThunderEgg_HDRS DFTPatchSolver.h)
  target_sources(ThunderEgg PRIVATE DFTPatchSolver.cpp)
endif(TARGET LAPACK::LAPACK AND TARGET BLAS::BLAS)

if(TARGET FFTW::FFTW)

  list(APPEND ThunderEgg_HDRS FFTWPatchSolver.h)
  target_sources(ThunderEgg PRIVATE FFTWPatchSolver.cpp)

endif(TARGET FFTW::FFTW)

if(TARGET PETSC::PETSC)
  if(TARGET FFTW::FFTW)
    list(APPEND ThunderEgg_HDRS FastSchurMatrixAssemble2D.h)
    target_sources(ThunderEgg PRIVATE
                FastSchurMatrixAssemble2D.cpp)

    list(APPEND ThunderEgg_HDRS FastSchurMatrixAssemble3D.h)
    target_sources(ThunderEgg PRIVATE
                FastSchurMatrixAssemble3D.cpp)
  endif(TARGET FFTW::FFTW)

  list(APPEND ThunderEgg_HDRS MatrixHelper.h)
  target_sources(ThunderEgg PRIVATE MatrixHelper.cpp)

  list(APPEND ThunderEgg_HDRS MatrixHelper2d.h)
  target_sources(ThunderEgg PRIVATE MatrixHelper2d.cpp)

endif(TARGET PETSC::PETSC)

# -- install public headers

install(
  FILES ${ThunderEgg_HDRS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ThunderEgg/Poisson
)

# -- add list of headers to all

list(TRANSFORM ThunderEgg_HDRS PREPEND "ThunderEgg/Poisson/")
list(APPEND ALL_ThunderEgg_HDRS ${ThunderEgg_HDRS})

# -- propogate to parent scope

set(ALL_ThunderEgg_HDRS ${ALL_ThunderEgg_HDRS} PARENT_SCOPE)