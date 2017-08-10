CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(ImatiSTL)

SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Os")
SET(CMAKE_LIBRARY_PATH /opt/local/lib ${CMAKE_LIBRARY_PATH})

# We want to test the hybrid kernel.
ADD_DEFINITIONS(-DUSE_HYBRID_KERNEL -DUSE_CGAL_LAZYNT)

# Add CGAL support.
SET(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL
    "Disable CGAL from overwriting my cmake flags")
FIND_PACKAGE(CGAL)
INCLUDE(${CGAL_USE_FILE})

SET(LIBRARY_OUTPUT_PATH    ${PROJECT_SOURCE_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
MAKE_DIRECTORY(${LIBRARY_OUTPUT_PATH})
MAKE_DIRECTORY(${EXECUTABLE_OUTPUT_PATH})

INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/include/Kernel/
    ${PROJECT_SOURCE_DIR}/include/TMesh/
    ${PROJECT_SOURCE_DIR}/include/ImatiSTL/)

SET(KERNEL_INC
    ${PROJECT_SOURCE_DIR}/include/Kernel/basics.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/coordinates.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/graph.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/heap.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/jqsort.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/list.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/matrix.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/point.h
    ${PROJECT_SOURCE_DIR}/include/Kernel/tmesh_kernel.h)

SET(TMESH_INC
    ${PROJECT_SOURCE_DIR}/include/TMesh/detectIntersections.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/edge.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/edgeHeap.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/faceClustering.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/simplification.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/tin.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/tmesh.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/triangle.h
    ${PROJECT_SOURCE_DIR}/include/TMesh/vertex.h)

SET(IMATISTL_INC
    ${PROJECT_SOURCE_DIR}/include/ImatiSTL/cutIntersections.h
    ${PROJECT_SOURCE_DIR}/include/ImatiSTL/dijkstraGraph.h
    ${PROJECT_SOURCE_DIR}/include/ImatiSTL/imatistl.h
    ${PROJECT_SOURCE_DIR}/include/ImatiSTL/marchIntersections.h
    ${PROJECT_SOURCE_DIR}/include/ImatiSTL/trimesh.h)

SET(ALGORITHM_SRC
    ${PROJECT_SOURCE_DIR}/src/Algorithms/checkAndRepair.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/delaunize.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/detectIntersections.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/faceClustering.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/featureRecover.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/holeFilling.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/simplification.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/smoothing.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/spherize.cpp
    ${PROJECT_SOURCE_DIR}/src/Algorithms/subdivision.cpp
    )

SET(IMATISTL_SRC
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/computeOuterHull.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/cutIntersections.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/dijkstraGraph.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/epsilonSampling.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/marchIntersections.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/offset.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/test.cpp
    ${PROJECT_SOURCE_DIR}/src/ImatiSTL/uniformRemesh.cpp
    )

SET(KERNEL_SRC
    ${PROJECT_SOURCE_DIR}/src/Kernel/orientation.c
    ${PROJECT_SOURCE_DIR}/src/Kernel/coordinates.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/graph.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/heap.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/jqsort.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/list.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/matrix.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/point.cpp
    ${PROJECT_SOURCE_DIR}/src/Kernel/tmesh.cpp
    )

SET(MESHFIX_SRC
    #${PROJECT_SOURCE_DIR}/src/MeshFix/meshfix.cpp 
    )

SET(TMESH_SRC
    ${PROJECT_SOURCE_DIR}/src/TMesh/edge.cpp
    ${PROJECT_SOURCE_DIR}/src/TMesh/edgeHeap.cpp
    ${PROJECT_SOURCE_DIR}/src/TMesh/io.cpp
    ${PROJECT_SOURCE_DIR}/src/TMesh/tin.cpp
    ${PROJECT_SOURCE_DIR}/src/TMesh/triangle.cpp
    ${PROJECT_SOURCE_DIR}/src/TMesh/vertex.cpp
    )

ADD_EXECUTABLE(ImatiSTL
    ${ALGORITHM_SRC}
    ${KERNEL_INC} ${KERNEL_SRC}
    ${TMESH_INC} ${TMESH_SRC}
    ${MESHFIX_INC} ${MESHFIX_SRC}
    ${IMATISTL_INC} ${IMATISTL_SRC}
    )
