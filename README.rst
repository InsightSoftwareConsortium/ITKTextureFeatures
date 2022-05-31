ITKTextureFeatures
==================

.. image:: https://github.com/InsightSoftwareConsortium/ITKTextureFeatures/workflows/Build,%20test,%20package/badge.svg
    :alt:    Build Status

Overview
--------

This repository contains `ITK <https://itk.org>`_ filters to estimate
texture feature maps from N-dimensional images.

For more information, see the `Insight Journal article <https://hdl.handle.net/10380/3574>`_::

  Vimort J., McCormick M., Budin F., Paniagua B.
  Computing Textural Feature Maps for N-Dimensional images
  McCormick M.
  The Insight Journal. January-December. 2017.
  https://hdl.handle.net/10380/3574
  https://insight-journal.org/browse/publication/985

Installation
------------

Python
^^^^^^

Binary `Python packages <https://pypi.python.org/pypi/itk-texturefeatures>`_
are available for macOS, Linux, and Windows. They can be installed with::

  python -m pip install --upgrade pip
  python -m pip install itk-texturefeatures

3D Slicer Extension
^^^^^^^^^^^^^^^^^^^

The module functionality is also available in the `3D Slicer
<https://slicer.org>`_ desktop application. Install the *BoneTextureExtension*
in the `3D Slicer extension manager
<https://www.slicer.org/wiki/Documentation/Nightly/SlicerApplication/ExtensionsManager>`_.
`Additional documentation
<https://raw.githubusercontent.com/Kitware/BoneTextureExtension/master/Docs/BoneTextureExtensionTutorial_2017.pdf>`_
is available for `the extension
<https://github.com/Kitware/BoneTextureExtension>`_.

C++
^^^

Since ITK 4.13.0, this module is available in the ITK source tree as a Remote
module. To enable it, set::

  Module_TextureFeatures:BOOL=ON

in ITK's CMake build configuration.

License
-------

The source code is distributed under the Apache 2.0 License. Please see LICENSE file for details.

Acknowledgements
----------------

This work was supported by the National Institute of Health (NIH) National
Institute for Dental and Craniofacial Research (NIDCR) grant R21DE025306
(Textural Biomarkers of Arthritis for the Subchondral Bone in the
Temporomandibular Joint) and NIDCR grant R01DE024450 (Quantification of 3D
bony Changes in Temporomandibular Joint Osteoarthritis).
