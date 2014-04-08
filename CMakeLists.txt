cmake_minimum_required(VERSION 2.8)
project(chrono_models)
include(SBELUtils.cmake)

#set path for local find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

SET(ENABLE_UNIT_GPU FALSE   CACHE   BOOL    "Turn on for GPU models")
SET(ENABLE_UNIT_MPI FALSE   CACHE   BOOL    "Turn on for GPU models")

# Additional packages
enable_openmp_support()
#enable_mpi_support()
#FIND_PACKAGE(hdf5)
FIND_PACKAGE(GLM)
# User specified path to ChronoEngine
SET(CHRONO_INC       ""   CACHE PATH "Where is ChronoEngine Source located")
SET(CHRONO_PARALLEL_INC       ""   CACHE PATH "Where is ChronoEngine Source located")
SET(CHRONO_OPENGL_INC  ""   CACHE PATH "Where is CHRONOOPENGL located")






# Include directories based on user path to ChronoEngine
SET(INC_DIRS ${CHRONO_INC} "${CHRONO_INC}/collision" "${CHRONO_INC}/collision/bullet" ${CHRONO_PARALLEL_INC} ${CHRONO_OPENGL_INC} ${HDF5_INCLUDE_DIRS})

# User specified path to ChronoEngine Libraries
SET(CHRONO_LIB_PATH  ""   CACHE PATH "Where is CHRONOLIB located")
SET(CHRONO_PARALLEL_LIB_PATH  ""   CACHE PATH "Where is CHRONOPARALLEL located")
SET(CHRONO_OPENGL_LIB_PATH  ""   CACHE PATH "Where is CHRONOOPENGL located")
# Find required Libraries
FIND_LIBRARY(CHRONO_LIB NAMES ChronoEngine PATHS ${CHRONO_LIB_PATH} REQUIRED)


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(GCC_OPT "sse" CACHE STRING "gcc optimizations")
	SET_PROPERTY(CACHE GCC_OPT PROPERTY STRINGS "mmx" "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2" "avx")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m${GCC_OPT}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m${GCC_OPT}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	SET(INTEL_OPT "SSE2" CACHE STRING "intel optimizations")
	SET_PROPERTY(CACHE INTEL_OPT PROPERTY STRINGS "SSE2" "SSE3" "SSSE3" "SSE4.1" "SSE4.2" "AVX" "CORE-AVX-I" "CORE-AVX2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -x${INTEL_OPT}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x${INTEL_OPT}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()


# Set unit-specific options
IF(ENABLE_UNIT_GPU)

    enable_cuda_support()
    LIST(APPEND INC_DIRS "${CUDA_TOOLKIT_ROOT_DIR}/include" "${CUDA_SDK_ROOT_DIR}/common/inc" )
    FIND_LIBRARY(CHRONO_LIB_PARALLEL NAMES ChronoEngine_Parallel PATHS ${CHRONO_PARALLEL_LIB_PATH} REQUIRED)
    #FIND_LIBRARY(CHRONO_LIB_POSTPROCESS NAMES ChronoEngine_POSTPROCESS PATHS ${CHRONO_LIB_PATH} REQUIRED)
    FIND_LIBRARY(CHRONO_LIB_OPENGL NAMES ChronoEngine_OPENGL PATHS ${CHRONO_OPENGL_LIB_PATH} REQUIRED)

    SET(GPU_LIBS
        #${CHRONO_LIB_POSTPROCESS}
        ${CHRONO_LIB_PARALLEL}
        ${CHRONO_LIB_OPENGL}
        ${HDF5_LIBRARIES}
        ${MPI_LIBRARIES}
    )
ENDIF()

IF(ENABLE_UNIT_MPI)
    enable_mpi_support()

    FIND_LIBRARY(CHRONO_LIB_MPI NAMES ChronoEngine_MPI PATHS ${CHRONO_LIB_PATH} REQUIRED)

    SET(MPI_LIBS
        ${CHRONO_LIB_MPI}
    )
ENDIF()

include_directories(${INC_DIRS})

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

# Add executables
IF (ENABLE_UNIT_GPU)
    #ADD_EXECUTABLE(ball_drop models/gpu/ball_drop.cpp)		
    #TARGET_LINK_LIBRARIES(ball_drop	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(convergence models/gpu/convergence.cpp)		
    TARGET_LINK_LIBRARIES(convergence	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(anchor models/gpu/anchor.cpp)		
    TARGET_LINK_LIBRARIES(anchor	${CHRONO_LIB} ${GPU_LIBS})
    
    #ADD_EXECUTABLE(aor models/gpu/aor.cpp)		
    #TARGET_LINK_LIBRARIES(aor	${CHRONO_LIB} ${GPU_LIBS})
    
    #ADD_EXECUTABLE(diagonal_impact models/gpu/diagonal_impact.cpp)		
    #TARGET_LINK_LIBRARIES(diagonal_impact	${CHRONO_LIB} ${GPU_LIBS})
    
    #ADD_EXECUTABLE(shaker models/gpu/shaker.cpp)		
    #TARGET_LINK_LIBRARIES(shaker	${CHRONO_LIB} ${GPU_LIBS})
    
    #ADD_EXECUTABLE(character_walk models/gpu/character_walk.cpp)		
    #TARGET_LINK_LIBRARIES(character_walk	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(trough models/gpu/trough.cpp)		
    TARGET_LINK_LIBRARIES(trough	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(density models/gpu/density.cpp)		
    TARGET_LINK_LIBRARIES(density	${CHRONO_LIB} ${GPU_LIBS})

    #ADD_EXECUTABLE(fluid_tank models/gpu/fluid_tank.cpp)		
    #TARGET_LINK_LIBRARIES(fluid_tank	${CHRONO_LIB} ${GPU_LIBS})

    #ADD_EXECUTABLE(monkey models/gpu/monkey.cpp)		
    #TARGET_LINK_LIBRARIES(monkey	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(sticky models/gpu/sticky.cpp)		
    TARGET_LINK_LIBRARIES(sticky	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(benchmark models/gpu/benchmark.cpp)		
    TARGET_LINK_LIBRARIES(benchmark	${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(foam models/gpu/foam.cpp)		
    TARGET_LINK_LIBRARIES(foam	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(extruder models/gpu/extruder.cpp)		
    TARGET_LINK_LIBRARIES(extruder	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(dragon models/gpu/dragon.cpp)		
    TARGET_LINK_LIBRARIES(dragon	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(teddy models/gpu/teddy.cpp)		
    TARGET_LINK_LIBRARIES(teddy	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(waterwheel models/gpu/waterwheel.cpp)		
    TARGET_LINK_LIBRARIES(waterwheel	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(deform models/gpu/deform.cpp)		
    TARGET_LINK_LIBRARIES(deform	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(drag_anchor models/gpu/drag_anchor.cpp)		
    TARGET_LINK_LIBRARIES(drag_anchor	${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(mixer models/gpu/mixer.cpp)		
    TARGET_LINK_LIBRARIES(mixer	${CHRONO_LIB} ${GPU_LIBS})
    
    
    ADD_EXECUTABLE(math models/gpu/math.cpp)		
    TARGET_LINK_LIBRARIES(math	${CHRONO_LIB} ${GPU_LIBS})
    
    #ADD_EXECUTABLE(rollers models/gpu/rollers.cpp)		
    #TARGET_LINK_LIBRARIES(rollers	${CHRONO_LIB} ${GPU_LIBS})

    #ADD_EXECUTABLE(sls models/gpu/sls.cpp)		
    #TARGET_LINK_LIBRARIES(sls	${CHRONO_LIB} ${GPU_LIBS})

    #ADD_EXECUTABLE(fording models/gpu/fording.cpp)      
    #TARGET_LINK_LIBRARIES(fording   ${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(fibers models/gpu/fibers.cpp)      
    TARGET_LINK_LIBRARIES(fibers   ${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(string models/gpu/string.cpp)      
    TARGET_LINK_LIBRARIES(string   ${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(tippytop models/gpu/tippytop.cpp)      
    TARGET_LINK_LIBRARIES(tippytop   ${CHRONO_LIB} ${GPU_LIBS})

	ADD_EXECUTABLE(mixer_fast models/gpu/mixer_fast.cpp)      
    TARGET_LINK_LIBRARIES(mixer_fast   ${CHRONO_LIB} ${GPU_LIBS})
    
     ADD_EXECUTABLE(snowball_smash models/gpu/snowball_smash.cpp)      
    TARGET_LINK_LIBRARIES(snowball_smash   ${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(cones models/gpu/cones)      
    TARGET_LINK_LIBRARIES(cones   ${CHRONO_LIB} ${GPU_LIBS})
    
     ADD_EXECUTABLE(plow models/gpu/plow)      
    TARGET_LINK_LIBRARIES(plow   ${CHRONO_LIB} ${GPU_LIBS})

    ADD_EXECUTABLE(pipe models/gpu/pipe)      
    TARGET_LINK_LIBRARIES(pipe   ${CHRONO_LIB} ${GPU_LIBS})
    
    ADD_EXECUTABLE(convergence_serial models/generic/convergence_serial)      
    TARGET_LINK_LIBRARIES(convergence_serial   ${CHRONO_LIB} ${GPU_LIBS})
    
    
    #ADD_EXECUTABLE(walker models/generic/demo_walker_nonIrr.cpp)		
    #TARGET_LINK_LIBRARIES(walker	${CHRONO_LIB} ${GPU_LIBS})

    #ADD_EXECUTABLE(cylinder_drop models/gpu/cylinder_drop.cpp)		
    #TARGET_LINK_LIBRARIES(cylinder_drop	${CHRONO_LIB} ${GPU_LIBS})
     ADD_EXECUTABLE(adams_comparison models/gpu/adams_comparison.cpp)		
    TARGET_LINK_LIBRARIES(adams_comparison	${CHRONO_LIB} ${GPU_LIBS})
ENDIF(ENABLE_UNIT_GPU)

IF (ENABLE_UNIT_MPI)
    ADD_EXECUTABLE(rover_enormous models/mpi/rover_enormous.cpp)
    TARGET_LINK_LIBRARIES(rover_enormous ${CHRONO_LIB} ${MPI_LIBS})
ENDIF()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
