itk_module_test()

set(ITKTextureFeaturesTests
  )

CreateTestDriver(MITKTextureFeatures "${ITKTextureFeatures-Test_LIBRARIES}" "${ITKTextureFeaturesTests}")

