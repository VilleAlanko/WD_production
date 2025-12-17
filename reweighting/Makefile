# Compiler
INIH_DIR = /home/alankovh/Documents/inih
CXX = g++
CC  = gcc
CXXFLAGS = -O2 -std=c++17 -I$(INIH_DIR) -I$(INIH_DIR)/cpp $(shell lhapdf-config --cxxflags)
LDFLAGS  = $(shell lhapdf-config --ldflags)

# Common INI sources
SRCS_CPP := $(INIH_DIR)/cpp/INIReader.cpp
SRCS_C   := $(INIH_DIR)/ini.c

# Source files for executables
SRC1 = src/new_PDFs.cpp
SRC2 = src/strangeness_asymmetry_uncertainties.cpp

# Targets
TARGET1 = new_PDFs
TARGET2 = strangeness_asymmetry_uncertainties

# Object files
OBJS_COMMON = $(SRCS_CPP:.cpp=.o) $(SRCS_C:.c=.o)
OBJS1 = $(SRC1:.cpp=.o) $(OBJS_COMMON)
OBJS2 = $(SRC2:.cpp=.o) $(OBJS_COMMON)

# Default target
all: $(TARGET1) $(TARGET2)

# Build rules
$(TARGET1): $(OBJS1)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(TARGET2): $(OBJS2)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compilation rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) -O2 -c $< -o $@

# Clean rule
clean:
	rm -f $(TARGET1) $(TARGET2) $(OBJS_COMMON) $(SRC1:.cpp=.o) $(SRC2:.cpp=.o)

