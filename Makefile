# -*- Makefile -*-
#make BUILD_DIR=./build/darwin TARGET_EXEC=QDLC-3.3.1.out TARGET_DIR=/Users/davidbauch/bin -j12
TARGET_EXEC := final_program.exe
TARGET_DIR := ./build
BUILD_DIR := ./build
SRC_DIRS := ./source ./external/ALGLIB
COMPILER = $(CXX)

# Find all the C and C++ files we want to compile
# Note the single quotes around the * expressions. Make will incorrectly expand these otherwise.
SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')


# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := include/ external/ #$(shell find $(SRC_DIRS) -type d)
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

LIB_FLAGS = -Wno-volatile -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# ADjust os dependent stuff
ifeq ($(OS),Windows_NT)
	LIB_FLAGS += -std=c++2a
	#BUILD_DIR += /win
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LIB_FLAGS += -std=c++20
		#BUILD_DIR += /darwin
		COMPILER = g++-10
	endif
endif

# String substitution for every C/C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
CPPFLAGS := $(INC_FLAGS) $(LIB_FLAGS)

# The final build step.
$(TARGET_DIR)/$(TARGET_EXEC): $(OBJS)
	@echo Compiling Main Program into $@
	@$(COMPILER) main.cpp -o $@ $(OBJS) $(CPPFLAGS)

# Build step for C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	@echo Compiling $@
	@mkdir -p $(dir $@)
	@$(COMPILER) $(CPPFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -r $(BUILD_DIR)