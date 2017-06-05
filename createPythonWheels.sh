#!/bin/bash

wget -L https://data.kitware.com/api/v1/file/592dd8068d777f16d01e1a92/download -O zstd-1.2.0-linux.tar.gz
gunzip -d zstd-1.2.0-linux.tar.gz
tar xf zstd-1.2.0-linux.tar

curl -L https://github.com/InsightSoftwareConsortium/ITKPythonBuilds/releases/download/latest/ITKPythonBuilds-linux.tar.zst -O ITKPythonBuilds-linux.tar.zst
./zstd-1.2.0-linux/bin/unzstd ITKPythonBuilds-linux.tar.zst -o ITKPythonBuilds-linux.tar
tar xf ITKPythonBuilds-linux.tar

./ITKPythonPackage/scripts/dockcross-manylinux-build-module-wheels.sh
