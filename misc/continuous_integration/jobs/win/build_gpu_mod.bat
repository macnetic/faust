rem should be in faust project root dir
cd gpu_mod
rem avoid interference between CUDA11 and CUDA12 by removing/setting VS build files
set VS_PATH=C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\MSBuild\Microsoft\VC\v160\BuildCustomizations
copy /Y "%VS_PATH%\bakCUDA 11.4.xml"  "%VS_PATH%\CUDA 11.4.xml"
copy /Y "%VS_PATH%\bakCUDA 11.4.props"  "%VS_PATH%\CUDA 11.4.props"
copy /Y "%VS_PATH%\bakCUDA 11.4.targets"  "%VS_PATH%\CUDA 11.4.targets"
del /Q "%VS_PATH%\CUDA 12.1.xml" "%VS_PATH%\CUDA 12.1.props" "%VS_PATH%\CUDA 12.1.targets"
if NOT EXIST build-cu11.4 (mkdir build-cu11.4) else (rmdir /S /Q build-cu11.4 & mkdir build-cu11.4)
cd build-cu11.4
cmake -G "Visual Studio 16 2019" -DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4/bin/nvcc.exe" -DCUDA_TOOLKIT_INCLUDE="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4/include" -DCUDA_TOOLKIT_ROOT_DIR="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.4" ..
cmake --build . --config %BUILD_CONFIG%
move %BUILD_CONFIG%\gm.dll .
cd ..
copy /Y "%VS_PATH%\bakCUDA 12.1.xml"  "%VS_PATH%\CUDA 12.1.xml"
copy /Y "%VS_PATH%\bakCUDA 12.1.props"  "%VS_PATH%\CUDA 12.1.props"
copy /Y "%VS_PATH%\bakCUDA 12.1.targets"  "%VS_PATH%\CUDA 12.1.targets"
del /Q "%VS_PATH%\CUDA 11.4.xml" "%VS_PATH%\CUDA 11.4.props" "%VS_PATH%\CUDA 11.4.targets"
if NOT EXIST build-cu12.1 (mkdir build-cu12.1) else (rmdir /S /Q build-cu12.1 & mkdir build-cu12.1)
cd build-cu12.1
cmake -G "Visual Studio 16 2019" -DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1/bin/nvcc.exe" -DCUDA_TOOLKIT_INCLUDE="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1/include" -DCUDA_TOOLKIT_ROOT_DIR="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1" ..
cmake --build . --config %BUILD_CONFIG%
move %BUILD_CONFIG%\gm.dll .
cd ..\..
