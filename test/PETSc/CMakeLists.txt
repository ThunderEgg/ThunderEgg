if(TARGET PETSc::PETSc)

    target_sources(unit_tests_mpi1 PRIVATE MatShellCreator_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE MatWrapper_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE PCShellCreator_MPI1.cpp)

endif(TARGET PETSc::PETSc)