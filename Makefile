# ==========================
# include stuff from Bridge++ Makefile

root = $(HOME)/LQCD/SC.bridge.HAL

include $(root)/Makefile

# libraries and flags from Bridge++ makefile
LDLIBS   := $(subst -Lextra,  -L$(root)/extra,  $(LDLIBS))

CPPFLAGS := $(subst -Isrc,    -I$(root)/src,    $(CPPFLAGS))
CPPFLAGS := $(subst -Iextra/, -I$(root)/extra/, $(CPPFLAGS))

CXXFLAGS := $(subst -Isrc,    -I$(root)/src,    $(CXXFLAGS))
CXXFLAGS := $(subst -Iextra/, -I$(root)/extra/, $(CXXFLAGS))

# Bridge++ objects
osrc_dirs= $(subst src, $(root)/build/src, $(filter-out src/Main, $(src_dirs)))
objects =  $(wildcard $(addsuffix /*.o, $(osrc_dirs)))


# ===========================
# makefile for my files

# source and build directories
SRC_DIR   = src
SRC_DIR_VOJTA = $(SRC_DIR)/src_VOJTA
SRC_DIR_HAL = $(SRC_DIR)/src_HAL
SRC_DIR_EXTERNAL = $(SRC_DIR)/src_EXTERNAL

INCL_DIR = -I$(SRC_DIR)
INCL_DIR += -I$(SRC_DIR_VOJTA)
INCL_DIR += -I$(SRC_DIR_HAL)
INCL_DIR += -I$(SRC_DIR_EXTERNAL)

BUILD_DIR = build


# HAL library
HAL_DIR = $(HOME)/LQCD/SC.hadron_force.2014-04-24

INCLUDES_HAL = -I$(HOME)/LQCD/SC/include #-I$(HAL_DIR) -I$(HAL_DIR)/include 
HAL_LIBS =  $(HOME)/LQCD/SC/lib/libfftw3.a #$(HAL_DIR)/hal.a 


# compiler commands if different from Bridge++
# CXX = /usr/bin/mpicxx


# compiler flags if different from Bridge++
#CFLAGS   = -O3 -fopenmp -w  # -O3 optimize, -fopenmp use openMP, -w suppress warnings
#CXXFLAGS += -O3 -fopenmp -w
#LDFLAGS  =  -O3 -fopenmp -w


# objects of my code
SRCS =
SRCS += $(wildcard $(addsuffix /*.cpp, $(SRC_DIR)))
SRCS += $(wildcard $(addsuffix /*.cpp, $(SRC_DIR_VOJTA)))
SRCS += $(wildcard $(addsuffix /*.cpp, $(SRC_DIR_HAL)))
SRCS += $(wildcard $(addsuffix /*.cpp, $(SRC_DIR_EXTERNAL)))

BIN_OBJS =
BIN_OBJS += $(subst $(SRC_DIR)/,$(BUILD_DIR)/, $(subst .cpp,.o, $(SRCS)))
 

BINS =
BINS += $(BUILD_DIR)/V_LQCD.x


All:  msg2 build_dir $(BINS) $(BIN_OBJS)

build_dir:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(subst $(SRC_DIR)/,$(BUILD_DIR)/, $(SRC_DIR_VOJTA))
	mkdir -p $(subst $(SRC_DIR)/,$(BUILD_DIR)/, $(SRC_DIR_HAL))
	mkdir -p $(subst $(SRC_DIR)/,$(BUILD_DIR)/, $(SRC_DIR_EXTERNAL))


# compilation of HAL code
#$(BUILD_DIR)/$(SRC_DIR_HAL)/hal_run_V.o : $(SRC_DIR_HAL)/hal_run_V.cpp
#	$(CXX)  $(CXXFLAGS)  -c $< -o $@

#compilation of the rest
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES_HAL) $(INCL_DIR) -c -o $@ $<


$(BINS): $(BIN_OBJS) 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ \
	$(objects) \
	$(LDLIBS) \
        $(HAL_LIBS)


Clean: 
	rm -fr build

msg2:
	@echo $(SRCS)
	@echo $(BIN_OBJS)
#	@echo $(objects)
#	@echo $(osrc_dirs)

