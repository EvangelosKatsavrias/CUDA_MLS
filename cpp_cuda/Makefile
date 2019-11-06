# ============ General variables ==============
cc       = gcc
cxx      = g++
cc_CUDA  = nvcc
jcc      = javac
fcc      = gfortran

cc_flags = 
cxx_flags = -Wall -O3 -std=c++14 -L/path_to_cuda_distribution/cuda/lib64 -lcudart -lcuda -fPIC


include_path = -I include -I include/Math -I include/DataClasses -I include/LeastSquares -I lib -I include/Mesh -I include/Topology -I include/MeshView -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64
src_path = src
obj_path = obj
lib_path = -L/usr/X11R6/lib64 -lQt5Widgets -lQt5Gui -lQt5Core -lGL -lpthread -lcublas -lcublas_device -lQt53DCore -lQt53DRender -lQt53DExtras -lQt53DInput
exe_path = exe


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
opencl_srcFolder    = $(src_path)
opencl_objFolder    = $(obj_path)

all: otherParts linkBuild 


otherParts:
	+$(MAKE) -C $(src_path)/misc
	+$(MAKE) -C $(src_path)/Math
	+$(MAKE) -C $(src_path)/DataClasses
	+$(MAKE) -C $(src_path)/LeastSquares
	+$(MAKE) -C $(src_path)/Mesh
	+$(MAKE) -C $(src_path)/MeshView


otherObjFiles = $(wildcard $(obj_path)/LeastSquares/*.o) $(wildcard $(obj_path)/LeastSquares/rObj/*.o) $(wildcard $(obj_path)/Math/*.o) $(wildcard $(obj_path)/Math/rObj/*.o) $(wildcard $(obj_path)/misc/*.o) $(wildcard $(obj_path)/misc/rObj/*.o) $(wildcard $(obj_path)/DataClasses/*.o) $(wildcard $(obj_path)/Mesh/*.o)  $(wildcard $(obj_path)/MeshView/*.o)


# ============ Main files ==============
srcType_mainProg = cpp
name_mainProg = test
cc_mainProg = $(cxx)



# ============= C files ==============
cSrcFiles = $(wildcard $(c_srcFolder)/*.c)
cObjFiles = $(patsubst $(c_srcFolder)/%.c, $(c_objFolder)/%.o, $(cSrcFiles))



# ============= Cpp files ==============
cppSrcFiles = $(wildcard $(cpp_srcFolder)/*.cpp)
cppObjFiles = $(patsubst $(cpp_srcFolder)/%.cpp, $(cpp_objFolder)/%.o, $(cppSrcFiles))



# ============= CUDA files ==============
cudaSrcFiles = $(wildcard $(cuda_srcFolder)/*.cu)
cudaObjFiles = $(patsubst $(cuda_srcFolder)/%.cu, $(cuda_objFolder)/%.o, $(cudaSrcFiles))



# ============ Compile rules =============
$(cuda_objFolder)/%.o: $(cuda_srcFolder)/%.cu
	$(cc_CUDA) $< -o $@ $(include_path) -c
$(c_objFolder)/%.o: $(c_srcFolder)/%.c
	$(cc) $< -o $@ $(include_path) -c $(c_flags)
$(cpp_objFolder)/%.o: $(cpp_srcFolder)/%.cpp
	$(cxx) $< -o $@ $(include_path) -c $(cxx_flags)
compileMainFile: $(src_path)/$(name_mainProg).$(srcType_mainProg)
	$(cc_mainProg) -c $^ -o $(obj_path)/$(name_mainProg).o $(include_path)


compileCudaFiles: $(cudaObjFiles)
compileCppFiles: $(cppObjFiles)
compileCFiles: $(cObjFiles)


# ============ Link and Build rules ===========
linkBuild: compileCFiles compileCppFiles compileCudaFiles compileMainFile
	$(cc_mainProg) -o $(exe_path)/$(name_mainProg) $(cppObjFiles)  $(otherObjFiles) $(cxx_flags) $(lib_path)



# ============ Run rules ====================
run: all
	cd $(exe_path); ./$(name_mainProg) Leaf/Leaf.msh 


profile_export: all
	cd $(exe_path); sudo nvprof -f -o ./profile.nvvp ./$(name_mainProg)
profile_summary: all
	cd $(exe_path); nvprof --print-summary ./$(name_mainProg)
profile_tracegpu: all
	cd $(exe_path); nvprof --print-gpu-trace ./$(name_mainProg)
profile_traceapi: all
	cd $(exe_path); nvprof --print-api-trace ./$(name_mainProg)
profile_metrics: all
	cd $(exe_path); nvprof --metrics all ./$(name_mainProg)
profile_analysisMetrics: all
	cd $(exe_path); nvprof --analysis-metrics -o ./analysis.nvprof ./$(name_mainProg)
profile_occupancy: all
	cd $(exe_path); nvprof --metrics achieved_occupancy ./$(name_mainProg)
profile_branchEfficiency: all
	cd $(exe_path); nvprof --metrics branch_efficiency ./$(name_mainProg)
profile_smEfficiency: all
	cd $(exe_path); nvprof --metrics sm_efficiency ./$(name_mainProg)
profile_events: all
	cd $(exe_path); nvprof --events all ./$(name_mainProg)
profile_bankConflicts: all
	cd $(exe_path); nvprof --events l1_shared_bank_conflict ./$(name_mainProg)





# ================= Cleaning rules ======================
clean:
	rm $(obj_path)/* -r



#===============================================================
echoTest:
	echo $(otherObjFiles)

#======== code archiver =============
archiveCode:
	tar -cf arch.tar src exe include lib

#======== code lines counter =============
countCodeLines:
	cloc ./ --exclude-lang=XMI,XML,INI,Assembly,Fortran\ 77
