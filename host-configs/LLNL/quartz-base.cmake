###############################################################################
#
# Base configuration for LC Quartz builds
# Calling configuration file must define the following CMAKE variables:
#
# MPI_HOME
#
###############################################################################

# Fortran
set(ENABLE_FORTRAN OFF CACHE BOOL "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER ${MPI_HOME}/bin/mpicc CACHE PATH "")
set(MPI_CXX_COMPILER  ${MPI_HOME}/bin/mpicxx CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")
set(MPIEXEC /usr/bin/srun CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

# PAPI (For TPL caliper builds)

set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX /usr/tce/packages/papi/papi-5.4.3 CACHE PATH "")

# OpenMP
set(ENABLE_OPENMP ON CACHE BOOL "")

# GEOSX specific options
set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)
set(ENABLE_PETSC OFF CACHE BOOL "Enables PETSc." FORCE)
set(ENABLE_HYPRE ON CACHE BOOL "" FORCE )

#set( GEOSX_LOCALINDEX_TYPE "int" CACHE STRING "" FORCE )
#set( GEOSX_GLOBALINDEX_TYPE "long long int" CACHE STRING "" FORCE )
#set( GEOSX_LA_INTERFACE "Hypre" CACHE STRING "Linear algebra interface to use in solvers" FORCE )

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
