# Executable and directory names
BIN_NAME = blacklight
BIN_DIR = bin
SRC_DIR = src
OBJ_DIR = obj
SRC_EXT = cpp
OBJ_EXT = o
DEP_EXT = d

# Compiler and flags
CC := g++
CPPFLAGS := -MMD -MP
CXXFLAGS := -std=c++17 -O3 -fopenmp -Wpedantic -Wall -Wextra -Wdouble-promotion -Wformat=2 \
	-Wformat-signedness -Wnull-dereference -Wmissing-include-dirs -Wunused -Wuninitialized \
	-Wstringop-truncation -Walloc-zero -Wduplicated-branches -Wduplicated-cond -Wshadow \
	-Wplacement-new=2 -Wundef -Wunused-macros -Wcast-qual -Wcast-align -Wconversion \
	-Wzero-as-null-pointer-constant -Wuseless-cast -Wextra-semi -Wsign-conversion -Wlogical-op \
	-Wmissing-declarations -Wredundant-decls -Wctor-dtor-privacy -Wnoexcept -Wnon-virtual-dtor \
	-Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wsign-promo
LDFLAGS :=
LDLIBS :=

# Lists of files to be considered
SRC_FILES := $(shell find $(SRC_DIR) -name "*.$(SRC_EXT)")
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(SRC_FILES:.$(SRC_EXT)=.$(OBJ_EXT))))
DEP_FILES := $(OBJ_FILES:.$(OBJ_EXT)=.$(DEP_EXT))
VPATH := $(sort $(dir $(SRC_FILES)))

# Default target
.PHONY : all
all : $(BIN_DIR)/$(BIN_NAME)

# Include dependency lists
-include $(DEP_FILES)

# Compile sources into objects
$(OBJ_DIR)/%.$(OBJ_EXT) : %.$(SRC_EXT)
	@echo compiling $(basename $(notdir $@))
	@$(CC) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Link objects into executable
$(BIN_DIR)/$(BIN_NAME) : $(OBJ_FILES)
	@echo linking $(basename $(notdir $@))
	@$(CC) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@

# Cleanup
.PHONY : clean
clean :
	rm -f $(OBJ_DIR)/*.$(OBJ_EXT)
	rm -f $(OBJ_DIR)/*.$(DEP_EXT)
	rm -f $(BIN_DIR)/$(BIN_NAME)
