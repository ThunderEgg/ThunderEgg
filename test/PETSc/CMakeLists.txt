if(TARGET PETSC::PETSC)

    target_sources(unit_tests_mpi1 PRIVATE MatShellCreator_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE MatWrapper_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE PCShellCreator_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE VecViewManager_MPI1.cpp)

    target_sources(unit_tests_mpi1 PRIVATE VecWrapper_MPI1.cpp)

endif(TARGET PETSC::PETSC)