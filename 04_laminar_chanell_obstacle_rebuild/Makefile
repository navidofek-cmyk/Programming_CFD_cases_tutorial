CXX      := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra
TARGET   := cfd_channel

SRCS := main.cpp mesh.cpp boundary.cpp solver.cpp io.cpp
OBJS := $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS) $(TARGET) *.csv

.PHONY: all clean
