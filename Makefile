# Blacklight Makefile

# Basic usage: "make"
# Optional arguments:
#   Targets:
#     <empty>: build bin/blacklight
#     all: same as <empty>
#     clean: remove bin/blacklight, as well as .o and .d files in obj/
#   Compiler options:
#     CXX=g++: use GNU g++ (default)
#     CXX=icpc: use Intel icpc
#   Other options:
#     -j: many processes to work in parallel (recommended except on busy resources)
#     -j <n>: use n processes to work in parallel
#     -j<n>: same as -j <n>

# Set executable and directory names
BIN_NAME = blacklight
BIN_DIR = bin
SRC_DIR = src
OBJ_DIR = obj
SRC_EXT = cpp
OBJ_EXT = o
DEP_EXT = d

# Set default compiler
CXX := g++

# Set compiler options: g++
ifeq ($(CXX), g++)
DEPENDENCY_OPTIONS := -MMD -MP
INCLUDE_OPTIONS :=
DIALECT_OPTIONS := -std=c++17 -fopenmp
OPTIMIZATION_OPTIONS := -O3 -flto -fno-math-errno -fno-signed-zeros -fno-trapping-math
WARNING_OPTIONS := -Wpedantic -Wall -Wextra -Wdouble-promotion -Wformat=2 -Wformat-signedness \
	-Wnull-dereference -Wmissing-include-dirs -Wunused -Wuninitialized -Wstringop-truncation \
	-Walloc-zero -Wduplicated-branches -Wduplicated-cond -Wshadow -Wundef -Wunused-macros \
	-Wcast-qual -Wcast-align -Wconversion -Wsign-conversion -Wlogical-op -Wmissing-declarations \
	-Wredundant-decls -Wlto-type-mismatch
WARNING_OPTIONS += -Wctor-dtor-privacy -Wnoexcept -Wnon-virtual-dtor -Wstrict-null-sentinel \
	-Wold-style-cast -Woverloaded-virtual -Wsign-promo -Wzero-as-null-pointer-constant \
	-Wplacement-new=2 -Wextra-semi -Wuseless-cast
LINKER_OPTIONS :=
LIBRARY_OPTIONS :=

# Set compiler options: icpc
else ifeq ($(CXX), icpc)
DEPENDENCY_OPTIONS := -MMD -MP
INCLUDE_OPTIONS :=
DIALECT_OPTIONS := -std=c++17 -qopenmp
OPTIMIZATION_OPTIONS := -O3 -ipo -xhost
WARNING_OPTIONS :=
LINKER_OPTIONS :=
LIBRARY_OPTIONS :=

# Handle unknown compiler
else
$(error "$(CXX)" not supported)
endif

# Set compiler flags
CPPFLAGS := $(DEPENDENCY_OPTIONS) $(INCLUDE_OPTIONS)
CXXFLAGS := $(DIALECT_OPTIONS) $(OPTIMIZATION_OPTIONS) $(WARNING_OPTIONS)
LDFLAGS := $(LINKER_OPTIONS)
LDLIBS := $(LIBRARY_OPTIONS)

# Set lists of files to be considered
SRC_FILES := $(shell find $(SRC_DIR) -name "*.$(SRC_EXT)")
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC_FILES:.$(SRC_EXT)=.$(OBJ_EXT))))
DEP_FILES := $(OBJ_FILES:.$(OBJ_EXT)=.$(DEP_EXT))
VPATH := $(sort $(dir $(SRC_FILES)))

# Set default target
.PHONY : all
all : $(BIN_DIR)/$(BIN_NAME)

# Include dependency lists
-include $(DEP_FILES)

# Compile sources into objects
$(OBJ_DIR)/%.$(OBJ_EXT) : %.$(SRC_EXT)
	@echo compiling $(basename $(notdir $@))
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Link objects into executable
$(BIN_DIR)/$(BIN_NAME) : $(OBJ_FILES)
	@echo linking $(basename $(notdir $@))
	@$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@

# Clean up
.PHONY : clean
clean :
	rm -f $(OBJ_DIR)/*.$(OBJ_EXT)
	rm -f $(OBJ_DIR)/*.$(DEP_EXT)
	rm -f $(BIN_DIR)/$(BIN_NAME)
