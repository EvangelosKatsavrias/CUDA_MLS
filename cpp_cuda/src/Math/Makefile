# ============ General variables ==============
cc       = gcc
cxx      = g++
cc_CUDA  = nvcc
jcc      = javac
fcc      = gfortran

cc_flags = 
cxx_flags = -std=c++14 -mavx -O3

include_folder = ../../include
include_path = -I $(include_folder) -I $(include_folder)/Math -I $(include_folder)/DataClasses -I $(include_folder)/LeastSquares -I lib
src_path = .
obj_path = ../../obj/Math
robj_path = ../../obj/Math/rObj
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
relocatedCudaFiles = ./weightingFunction_cuda.cu



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
	$(cc_CUDA) $< -o $@ $(include_path) -c
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



# ================= Cleaning rules ======================
clean:
	rm $(obj_path)/* -r



#===============================================================
echoTest:
	echo $(cudaSrcFiles)

