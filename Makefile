# -*- Makefile -*-
#MAKEFLAGS+="j 5"
TARGET_EXEC := a.out
TARGET_EXEC_LATEST := a.out
TARGET_DIR := ./build
BUILD_DIR := ./build
SRC_DIRS := ./source ./external/ALGLIB
COMPILER = $(CXX)

VERSION := 3.6.0

SRCS := $(shell find $(SRC_DIRS) -name "*.cpp")

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
	TARGET_EXEC_LATEST := QDLC.exe
	BUILD_DIR = ./build/win
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LIB_FLAGS += -std=c++2b
		TARGET_EXEC = QDLC-$(VERSION).out
		TARGET_EXEC_LATEST = QDLC.out
		TARGET_DIR = /Users/davidbauch/bin
		BUILD_DIR = ./build/darwin
		COMPILER = g++-12
	endif
	ifeq ($(UNAME_S),Linux)
		LIB_FLAGS += -std=c++2a
		TARGET_EXEC = QDLC-$(VERSION).out
		TARGET_EXEC_LATEST = QDLC.out
		TARGET_DIR = ./
		BUILD_DIR = ./build/linux
		COMPILER = g++
	endif
endif

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# Compile with UFLAG="-DLOG_DISABLE_L3" to disable L3
CPPFLAGS := $(INC_FLAGS) $(LIB_FLAGS) $(UFLAG)

# The final build step.
#$(TARGET_DIR)/$(TARGET_EXEC): $(OBJS)
main: $(OBJS)
	@echo Compiling Main Program into $(TARGET_DIR)/$(TARGET_EXEC), compiler is $(COMPILER), libs are $(LIB_FLAGS)
	@echo User Flags: $(UFLAG)
	@$(COMPILER) main.cpp -o $(TARGET_DIR)/$(TARGET_EXEC) $(OBJS) $(LIB_LINKS) $(CPPFLAGS)
	# ln -sf $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/QDLC
	cp $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/$(TARGET_EXEC_LATEST)
maindebug: $(OBJS)g
	@echo Compiling Main Program into $(TARGET_DIR)/$(TARGET_EXEC), compiler is $(COMPILER), libs are $(LIB_FLAGS)
	@echo User Flags: $(UFLAG)
	@$(COMPILER) main.cpp -g -o $(TARGET_DIR)/$(TARGET_EXEC) $(OBJS) $(LIB_LINKS) $(CPPFLAGS)
	# ln -sf $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/QDLC
	cp $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/$(TARGET_EXEC_LATEST)

# Build step for C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	@echo Compiling $@
	@mkdir -p $(dir $@)
	@$(COMPILER) $(CPPFLAGS) -c $< -o $@
$(BUILD_DIR)/%.cpp.og: %.cpp
	@echo Compiling $@
	@mkdir -p $(dir $@)
	@$(COMPILER) $(CPPFLAGS) -g -c $< -o $@

# .PHONY: clean
clean:
	rm -frv $(BUILD_DIR)

build: main #$(TARGET_DIR)/$(TARGET_EXEC)
all: clean build

#bake_before: 
#	@echo Generating cla_settings.h...
#	@xxd -i .\.settings.cla include/cla_settings.h
#bake_after:
#	@echo Removing cla_settings.h
#	@rm include/cla_settings.h
#bake: bake_before all bake_after