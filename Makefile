# Compiler
CC := g++

# Compiler flags
CFLAGS := -std=c++11 -Wall -Werror -Wextra -pedantic -march=native -O2 -ftree-vectorize -fopt-info-vec -ffast-math
#CFLAGS := -std=c++11 -Wall -Werror -Wextra -pedantic

# Source files
SRCS := Verlet.cpp

# Object files
OBJS := $(SRCS:.cpp=.o)

# Executable name
EXEC := simulation

# Rule to build the executable
$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Rule to compile source files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to clean the project
clean:
	rm -f $(OBJS) $(EXEC) energy.txt positions.xyz
