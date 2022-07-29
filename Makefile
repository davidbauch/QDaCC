# -*- Makefile -*-
#make BUILD_DIR=./build/darwin TARGET_EXEC=QDLC-3.3.1.out TARGET_DIR=/Users/davidbauch/bin -j12
#make BUILD_DIR=./build/win TARGET_EXEC=QDLC-3.3.1.exe TARGET_DIR=../../Threadhandler -j6
TARGET_EXEC := a.out
TARGET_DIR := ./build
BUILD_DIR := ./build
SRC_DIRS := ./source ./external/ALGLIB
COMPILER = $(CXX)

VERSION := 3.4.0

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')

INC_DIRS := include/ external/ #$(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

LIB_DIRS :=  #$(shell find $(SRC_DIRS) -type d)
LIB_FLAGS := $(addprefix -L,$(LIB_DIRS)) -O3 -DFMT_HEADER_ONLY -D_GLIBCXX_PARALLEL
LIB_LINKS := 

ifneq ($(COMPILER), /opt/intel/oneapi/compiler/latest/mac/bin/intel64/icpc)
	LIB_FLAGS += -fopenmp -Wno-volatile -lstdc++fs
else
	LIB_FLAGS += -qopenmp -w0 -static-intel
endif

# ADjust os dependent stuff
ifeq ($(OS),Windows_NT)
	LIB_FLAGS += -std=c++2a
	TARGET_DIR = ../../Threadhandler
	TARGET_EXEC := QDLC-$(VERSION).exe
	BUILD_DIR = ./build/win
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LIB_FLAGS += -std=c++20
		TARGET_EXEC = QDLC-$(VERSION).out
		TARGET_DIR = /Users/davidbauch/bin
		BUILD_DIR = ./build/darwin
		COMPILER = g++-12
	endif
	ifeq ($(UNAME_S),Linux)
		LIB_FLAGS += -std=c++2a
		TARGET_EXEC = QDLC-$(VERSION).out
		TARGET_DIR = ./
		BUILD_DIR = ./build/linux
		COMPILER = g++
	endif
endif

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# Compile with UFLAG="-DLOG_DISABLE_L3" to disable L3
CPPFLAGS := $(INC_FLAGS) $(LIB_FLAGS) $(UFLAG)

# The final build step.
$(TARGET_DIR)/$(TARGET_EXEC): $(OBJS)
	@echo Compiling Main Program into $@
	@echo User Flags: $(UFLAG)
	@$(COMPILER) main.cpp -o $@ $(OBJS) $(LIB_LINKS) $(CPPFLAGS)
	ln -sf $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/QDLC

# Build step for C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	@echo Compiling $@
	@mkdir -p $(dir $@)
	@$(COMPILER) $(CPPFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)
