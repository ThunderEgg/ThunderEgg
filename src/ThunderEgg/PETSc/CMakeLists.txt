if(PETSC_FOUND)

  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/MatWrapper.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/MatShellCreator.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/PCShellCreator.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/VecWrapper.h)
  list(APPEND ThunderEgg_HDRS ThunderEgg/PETSc/VecLocalDataManager.h)

endif(PETSC_FOUND)
