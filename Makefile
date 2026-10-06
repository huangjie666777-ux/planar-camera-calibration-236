CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Iinclude -Ithird_party/eigen3

SRC := src/camcal.cpp src/validation.cpp src/homography.cpp src/geometry_init.cpp src/projection.cpp src/optimizer.cpp
OBJ := $(SRC:.cpp=.o)

all: libcamcal236.a demo

libcamcal236.a: $(OBJ)
	ar rcs $@ $(OBJ)

demo: examples/demo.o libcamcal236.a
	$(CXX) $(CXXFLAGS) -o $@ examples/demo.o libcamcal236.a

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: demo
	./demo

clean:
	rm -f $(OBJ) examples/demo.o libcamcal236.a demo

.PHONY: all test clean

