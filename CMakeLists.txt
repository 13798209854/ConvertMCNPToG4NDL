#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ConvertMCNPToG4NDL)

# Parallel processing using TOP-C
option(TOPC_USE "Use TOP-C parallel architecture" OFF)
option(TOPC_USE_SEQ "Use TOP-C in sequential mode" OFF)
if(TOPC_USE)
  if(TOPC_USE_SEQ)
    set(PARMODE "--seq")
  else()
    set(PARMODE "--mpi")
  endif()

  add_definitions(-DTOPC_USE=1)
  # Get TOP-C libraries
  execute_process(COMMAND "topc-config" ${PARMODE} "--libs"
          OUTPUT_VARIABLE TOPC_LIBS)

  # Remove trailing whitespace
  string(STRIP ${TOPC_LIBS} TOPC_LIBS)

  # Get TOP-C compiler flags
  execute_process(COMMAND "topc-config" ${PARMODE} "--cflags"
          OUTPUT_VARIABLE TOPC_CFLAGS)
endif()

# Output library with MCNP6 settings
option(NeutHPMod_Use "Use modified NeutronHP format" OFF)
if(MCNP6Set_Use)
  add_definitions(-DMCNP6Set_Use=1)
  add_definitions(-DNeutHPMod_Use=1)
endif()

# Output library in the modified NeutronHP format
option(NeutHPMod_Use "Use modified NeutronHP format" OFF)
if(NeutHPMod_Use)
  add_definitions(-DNeutHPMod_Use=1)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc
           ${PROJECT_SOURCE_DIR}/src/*.icc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh
           ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(convertMCNPToG4NDL ConvertMCNPtoG4NDL.cc ${sources} ${headers})
target_link_libraries(convertMCNPToG4NDL ${TOPC_LIBS})


#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(ConvertMCNPToG4NDL DEPENDS convertMCNPToG4NDL)


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS convertMCNPToG4NDL DESTINATION bin)
