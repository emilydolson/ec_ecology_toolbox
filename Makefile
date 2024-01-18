# Project-specific settings
PROJECT := evo_comp_ecology
EMP_DIR := third-party/Empirical/include/emp

# Flags to use regardless of compiler
CFLAGS_all := -Wall -std=c++20 -fPIC -Wno-unused-function -I$(EMP_DIR)/

# Native compiler information
CXX := clang++
CXX_nat := $(CXX)
CFLAGS_nat := -O3 -DNDEBUG $(CFLAGS_all)
CFLAGS_nat_debug := -DEMP_TRACK_MEM -fprofile-arcs -ftest-coverage -g $(CFLAGS_all)
CFLAGS_nat_profile := -fprofile-arcs -ftest-coverage -pg $(CFLAGS_all)

# Emscripten compiler information
# CXX_web := emcc
# OFLAGS_web_all := -s TOTAL_MEMORY=67108864 --js-library $(EMP_DIR)/web/library_emp.js --js-library $(EMP_DIR)/web/d3/library_d3.js -s EXPORTED_FUNCTIONS="['_main', '_empCppCallback']" -s DISABLE_EXCEPTION_CATCHING=1 -s NO_EXIT_RUNTIME=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['cwrap', 'stringToUTF8', 'writeStringToMemory']"#--embed-file configs
# OFLAGS_web := -Oz -DNDEBUG
# OFLAGS_web_debug := -Oz -pedantic -Wno-dollar-in-identifier-extension

# CFLAGS_web := $(CFLAGS_all) $(OFLAGS_web) $(OFLAGS_web_all)
# CFLAGS_web_debug := $(CFLAGS_all) $(OFLAGS_web_debug) $(OFLAGS_web_all)

default: test

test: tests/test_interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat_debug) tests/test_interaction_networks.cc -o test
	time ./test

test-cag: tests/test_cag.cpp
	$(CXX_nat) $(CFLAGS_nat_debug) tests/test_cag.cpp -o test-cag
	time ./test-cag


profile: tests/test_interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat_profile) tests/test_interaction_networks.cc -o test
	./test

benchmark: tests/test_interaction_networks.cc
	$(CXX_nat) $(CFLAGS_nat) tests/test_interaction_networks.cc -o test
	./test benchmarks

wrapper:
	$(CXX_nat) -shared $(shell python3 -m pybind11 --includes)/ $(CFLAGS_all) source/wrapper.cc -o lexicase$(shell python3-config --extension-suffix)

clean:
	rm -f *~ source/*.o test *.gcda *.gcno

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
