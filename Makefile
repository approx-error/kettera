 # TODO: Figure out wheter the common_flags only need to be given at the compilation step
src_dir := src
obj_dir := obj
mod_dir := mod
inc_dir := /usr/local/include
build_dir := build

parameters_src := kinds.f90 exit_codes.f90 io_parameters.f90 sim_parameters.f90
functions_src := read_write.f90 alloc_dealloc.f90 fftw3.f90 calculate.f90 initialize.f90 iterate.f90
main_src := main.f90
src_files := $(parameters_src) $(functions_src) $(main_src)
exe_name := kettera

compiler := gfortran
common_flags := -Wall -Wextra -Wno-maybe-uninitialized -Wno-uninitialized -O2 
# XXX: The -Wno-maybe-uninitialized and -Wno-uninitialized flags are needed to suppress bug number 77504 on GCC Bugzilla that is resolved in
# gcc (and by extension gfortran) version 12.5 but since the default gcc/gfortran on my system is 11.4 it is not resolved for me
compilation_flags := -ffree-form -ffree-line-length-none -fimplicit-none -fbackslash -fmodule-private -std=f2018 -I "$(inc_dir)" -J "$(mod_dir)"
libraries := -lfftw3 -lm  

src_path := $(patsubst %.f90, $(src_dir)/%.f90, $(src_files))
obj_path := $(patsubst $(src_dir)/%.f90, $(obj_dir)/%.o, $(src_path))
mod_path := $(patsubst $(src_dir)/%.f90, $(mod_dir)/%.mod, $(src_path))
exe_path := $(build_dir)/$(exe_name)


.PHONY: all
all: $(mod_path) $(exe_path)
	
# Shorthand command to run the executable. Since the executable is a prerequisite,
# make will run the previous two block if it doesn't exist yet. thus running
# 'make run' is the same as running 'make' and then 'make run'
.PHONY: run
run: $(mod_path) $(exe_path)
	@echo Running the program
	./$(exe_path)

$(mod_path): | $(mod_dir)

$(mod_dir):
	mkdir -p $(mod_dir)

# The executable is the target and the object files are the prerequisites
# If the executable doesn't exist, this block will be executed
$(exe_path): $(obj_path) | $(build_dir)
	@echo Linking object files
	$(compiler) $(common_flags) $^ -o $@ $(libraries)

$(build_dir):
	mkdir -p $(build_dir)

$(obj_path): | $(obj_dir)

$(obj_dir):
	mkdir -p $(obj_dir)

# The object files are the target and the source files are the prerequisites
# If the object files don't exist, this block will be executed
$(obj_dir)/%.o: $(src_dir)/%.f90
	@echo Compiling source files
	$(compiler) $(common_flags) $(compilation_flags) -c $< -o $@

# Shorthand commands to delete the object files, modules and the final executable
.PHONY: clean deepclean
clean:
	-rm -f $(wildcard $(obj_dir)/*.o) $(wildcard $(mod_dir)/*.mod) $(exe_path)

deepclean:
	-rm -rf $(obj_dir) $(mod_dir) $(build_dir)
	
