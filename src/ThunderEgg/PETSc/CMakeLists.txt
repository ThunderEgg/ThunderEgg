if(TARGET PETSC::PETSC)

  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/MatWrapper.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/MatShellCreator.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/PCShellCreator.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/VecWrapper.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/VecLocalDataManager.h)

endif(TARGET PETSC::PETSC)
