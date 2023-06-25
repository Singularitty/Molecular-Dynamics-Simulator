TARGET = simulator

# Compiler
CC = g++

# Compilation flags.
CFLAGS = -Wall -Wextra -pedantic -march=native -O2 -ftree-vectorize -fopt-info-vec -ffast-math

# Source files.
SRC = src/simulator.cpp src/*/*.cpp
HEADERS = src/*/*.hpp

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC)

clean:
	$(RM) $(TARGET)
