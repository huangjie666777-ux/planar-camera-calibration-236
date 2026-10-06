CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Ithird_party/eigen3 -Iinclude

LIB_SRCS := src/validation.cpp src/homography.cpp src/initialization.cpp src/projection.cpp src/optimizer.cpp src/calibration.cpp
LIB_OBJS := $(LIB_SRCS:src/%.cpp=build/%.o)

.PHONY: all test clean

all: build/libcamcal236.a build/multiview_demo build/selftest

build:
	mkdir -p build

build/%.o: src/%.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/libcamcal236.a: $(LIB_OBJS)
	ar rcs $@ $^

build/multiview_demo: examples/multiview_demo.cpp build/libcamcal236.a
	$(CXX) $(CXXFLAGS) $< -Lbuild -lcamcal236 -o $@

build/selftest: tests/selftest.cpp build/libcamcal236.a
	$(CXX) $(CXXFLAGS) $< -Lbuild -lcamcal236 -o $@

test: all
	./build/selftest
	./build/multiview_demo

clean:
	rm -rf build
