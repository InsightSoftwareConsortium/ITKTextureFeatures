cmake_minimum_required(VERSION 3.9.5)
set(module computeFeatures)
project(${module})

find_package(ITK REQUIRED COMPONENTS TextureFeatures ITKIOImageBase)
include(${ITK_USE_FILE})

set(${module}_SRC
    computeGLRLMFeatures.cxx
)

add_executable(${module} ${${module}_SRC} )
target_link_libraries(${module} ${ITK_LIBRARIES})

set(INSTALL_RUNTIME_DESTINATION bin CACHE STRING "Install destination")

install(TARGETS ${module} DESTINATION ${INSTALL_RUNTIME_DESTINATION})
