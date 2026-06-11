cmake_minimum_required(VERSION 3.12)
project(ReadSTL VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# VTK 9: request ONLY the modules we use. A bare find_package(VTK) loads every
# module and tries to resolve all their external deps (OGG/Theora/...), which fails
# against ParaView's VTK. Classes are factory-registered via vtk_module_autoinit.
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    CommonMath
    CommonExecutionModel
    CommonTransforms
    FiltersCore
    FiltersGeneral
    FiltersGeometry
    FiltersModeling
    FiltersSources
    FiltersExtraction
    IOGeometry
    IOXML
    IOLegacy
    RenderingCore
    RenderingOpenGL2
    RenderingAnnotation
    RenderingFreeType
    InteractionStyle
    InteractionWidgets
)
message(STATUS "Using VTK ${VTK_VERSION}")

# OpenMP for the parallel pixel loop in visibilityTrace (optional but recommended).
find_package(OpenMP)

add_executable(ReadSTL
    main.cpp
    vtkRaytrace.cpp
    meshTools.cpp
)

target_link_libraries(ReadSTL PRIVATE ${VTK_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(ReadSTL PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled")
else()
    message(STATUS "OpenMP NOT found — building serial")
endif()

# Required for VTK >= 9 so object factories (IO, Rendering) initialise.
vtk_module_autoinit(
    TARGETS ReadSTL
    MODULES ${VTK_LIBRARIES}
)
