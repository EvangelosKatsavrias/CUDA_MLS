#   CUDA_MLS Framework
#
#   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
#
#   This file is part of the CUDA_MLS Framework.
#
#   CUDA_MLS Framework is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License version 3 as published by
#   the Free Software Foundation.
#
#   CUDA_MLS Framework is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
#
#   Contact Info:
#   Evangelos D. Katsavrias
#   email/skype: vageng@gmail.com
# -----------------------------------------------------------------------

# ============ General variables ==============
cc       = gcc
cxx      = g++
cc_CUDA  = nvcc
jcc      = javac
fcc      = gfortran

cc_flags = -mavx
cxx_flags = -std=c++14 -mavx -O3
cuda_flags = -std=c++11 -arch=sm_30

include_folder = ../../include
include_path = -I $(include_folder) -I $(include_folder)/Math -I $(include_folder)/DataClasses -I $(include_folder)/LeastSquares -I $(include_folder)/misc -I lib
src_path = .
obj_path = ../../obj/misc
robj_path = ../../obj/misc/rObj
lib_path = ../../lib
exe_path = ../../exe

cpp_srcFolder       = $(src_path)
cpp_objFolder       = $(obj_path)
c_srcFolder         = $(src_path)
c_objFolder         = $(obj_path)
java_srcFolder      = $(src_path)
java_objFolder      = $(obj_path)
fortran_srcFolder   = $(src_path)
fortran_objFolder   = $(obj_path)
cuda_srcFolder      = $(src_path)
cuda_objFolder      = $(obj_path)
cuda_robjFolder     = $(robj_path)
opencl_srcFolder    = $(src_path)
opencl_objFolder    = $(obj_path)

all: makeFolders compileCFiles compileCppFiles compileCudaFiles compileRelocatableCudaFiles

makeFolders:
	mkdir -p $(obj_path)
	mkdir -p $(robj_path)


# ============ Main files ==============
srcType_mainProg = cpp
name_mainProg = test
cc_mainProg = $(cxx)
relocatedCudaFiles = ./atomicOperations_doubles.cu 


# ============= C files ==============
cSrcFiles = $(wildcard $(c_srcFolder)/*.c)
cObjFiles = $(patsubst $(c_srcFolder)/%.c, $(c_objFolder)/%.o, $(cSrcFiles))



# ============= Cpp files ==============
cppSrcFiles = $(wildcard $(cpp_srcFolder)/*.cpp)
cppObjFiles = $(patsubst $(cpp_srcFolder)/%.cpp, $(cpp_objFolder)/%.o, $(cppSrcFiles))



# ============= CUDA files ==============
cudaSrcFiles1 = $(wildcard $(cuda_srcFolder)/*.cu)
cudaSrcFiles  = $(filter-out $(relocatedCudaFiles), $(cudaSrcFiles1))
cudaObjFiles = $(patsubst $(cuda_srcFolder)/%.cu, $(cuda_objFolder)/%.o, $(cudaSrcFiles))

# ============= CUDA relocatable files ==============
cudaRelocatableSrcFiles = $(relocatedCudaFiles)
cudaRelocatableObjFiles = $(patsubst %.cu, $(cuda_robjFolder)/%.o, $(cudaRelocatableSrcFiles))




# ============ Compile rules =============
$(cuda_objFolder)/%.o: $(cuda_srcFolder)/%.cu
	$(cc_CUDA) $< -o $@ $(include_path) -c $(cuda_flags)
$(cuda_robjFolder)/%.o: $(cuda_srcFolder)/%.cu
	$(cc_CUDA) -x cu $< -o $@ $(include_path) -dc $(cuda_flags)
$(c_objFolder)/%.o: $(c_srcFolder)/%.c
	$(cc) $< -o $@ $(include_path) -c $(c_flags)
$(cpp_objFolder)/%.o: $(cpp_srcFolder)/%.cpp
	$(cxx) $< -o $@ $(include_path) -c $(cxx_flags)
compileMainFile: $(src_path)/$(name_mainProg).$(srcType_mainProg)
	$(cc_mainProg) -c $^ -o $(obj_path)/$(name_mainProg).o $(include_path)


compileCudaFiles: $(cudaObjFiles)
compileRelocatableCudaFiles: $(cudaRelocatableObjFiles)
compileCppFiles: $(cppObjFiles)
compileCFiles: $(cObjFiles)



# ============ Link and Build rules ===========
linkBuild: compileCFiles compileCppFiles compileCudaFiles compileMainFile
	$(cc_mainProg) -o $(exe_path)/$(name_mainProg) $(cppObjFiles) $(cudaObjFiles) $(obj_path)/$(name_mainProg).o $(cxx_flags)


# ============ Run rules ====================
run:
	$(exe_path)/$(name_mainProg)


profile: linkBuild
	nvprof -o $(exe_path)/profile $(exe_path)/$(name_mainProg)
profile2: linkBuild
	nvprof -f $(exe_path)/$(name_mainProg)


runGridTest: gridDataTest.cpp
	$(cxx) gridDataTest.cpp $(include_path) -c $(cxx_flags)
	$(cxx) gridDataTest.o testGridProg.cpp -o ../../exe/test $(include_path) $(cxx_flags)
	cd ../../exe; ./test

# ================= Cleaning rules ======================
clean:
	rm $(obj_path)/* -r



#===============================================================
echoTest:
	echo $(cudaSrcFiles)

