# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component(MY_CURRENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(READ "${MY_CURRENT_DIR}/README.rst" DOCUMENTATION)

# itk_module() defines the module dependencies in ModuleTemplate
# ModuleTemplate depends on ITKCommon
# The testing module in ModuleTemplate depends on ITKTestKernel
# and ITKMetaIO(besides ModuleTemplate and ITKCore)
# By convention those modules outside of ITK are not prefixed with
# ITK.

if( NOT "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}" VERSION_LESS "4.12" )
  set(_GoogleTest_DEPENDS ITKGoogleTest)
endif()

# define the dependencies of the include module and the tests
itk_module(TextureFeatures
  DEPENDS
    ITKCommon
    ITKStatistics
    ITKImageGrid
    ITKMathematicalMorphology
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
    ITKImageIntensity
    ITKImageNoise
    ${_GoogleTest_DEPENDS}
  DESCRIPTION
    "${DOCUMENTATION}"
  EXCLUDE_FROM_DEFAULT
)
