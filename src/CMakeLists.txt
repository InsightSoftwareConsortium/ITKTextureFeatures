set(ITKTextureFeatures_SRCS
  )

add_library(ITKTextureFeatures ${ITKTextureFeatures_SRCS})
itk_module_link_dependencies()
itk_module_target(ITKTextureFeatures)
