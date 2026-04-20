CXX := g++
CXXFLAGS := -std=c++20 -O2 -Wall -Wextra -pedantic

TARGET := nozzle_solver
SRCS := main.cpp mesh.cpp geometry.cpp state.cpp solver.cpp io.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) solution.csv
