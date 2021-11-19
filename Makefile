# -*- Makefile -*-
#make BUILD_DIR=./build/darwin TARGET_EXEC=QDLC-3.3.1.out TARGET_DIR=/Users/davidbauch/bin -j12
#make BUILD_DIR=./build/win TARGET_EXEC=QDLC-3.3.1.exe TARGET_DIR=../../Threadhandler -j6
TARGET_EXEC := a.out
TARGET_DIR := ./build
BUILD_DIR := ./build
SRC_DIRS := ./source ./external/ALGLIB
COMPILER = $(CXX)

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')

INC_DIRS := include/ external/ #$(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

LIB_FLAGS = -Wno-volatile -O3 -DFMT_HEADER_ONLY -fopenmp -lstdc++fs

# ADjust os dependent stuff
ifeq ($(OS),Windows_NT)
	LIB_FLAGS += -std=c++2a
	TARGET_DIR = ../../Threadhandler
	TARGET_EXEC := QDLC-3.3.1.exe
	BUILD_DIR = ./build/win
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LIB_FLAGS += -std=c++20
		TARGET_EXEC = QDLC-3.3.1.out
		TARGET_DIR = /Users/davidbauch/bin
		BUILD_DIR = ./build/darwin
		COMPILER = g++-10
	endif
	ifeq ($(UNAME_S),Linux)
		LIB_FLAGS += -std=c++20
		TARGET_EXEC = QDLC-3.3.1.out
		TARGET_DIR = /Users/davidbauch/bin
		BUILD_DIR = ./build/linux
		COMPILER = g++
	endif
endif

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

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