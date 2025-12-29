# --- Molekel macOS Makefile (ARM64) ---

# Compilers
CC = gcc
CXX = g++
FC = gfortran

# Paths
BREW_PATH = /opt/homebrew
INCLUDES = -I$(BREW_PATH)/include -I. -I./mui
LIBS = -L$(BREW_PATH)/lib -lglui -ltiff -ljpeg -lglfw

# Frameworks for macOS
FRAMEWORKS = -framework OpenGL -framework GLUT -framework Cocoa

# Compiler Flags
CXXFLAGS = -Wno-deprecated-declarations -w -DLINUX $(INCLUDES)  -O3
CFLAGS = -Wno-deprecated-declarations -w $(INCLUDES) -std=gnu89  -O3
FFLAGS = -std=legacy -fno-range-check  -O3

# Source groups
MUI_SRCS = $(wildcard mui_source/*.c)
MUI_OBJS = $(MUI_SRCS:.c=.o)

CPP_SRCS = $(wildcard *.C)
CPP_OBJS = $(CPP_SRCS:.C=.o)

# Targets
all: libmui.a ms glutmolekel

# 1. Build the MUI library
libmui.a: $(MUI_OBJS)
	ar rc $@ $(MUI_OBJS)
	ranlib $@

# 2. Build the Fortran utility
ms: ms.f
	$(FC) $(FFLAGS) ms.f -o ms

# 3. Build the main executable
glutmolekel: $(CPP_OBJS) libmui.a
	$(CXX) $(CXXFLAGS) -o $@ $(CPP_OBJS) ./libmui.a $(LIBS) $(FRAMEWORKS)
	cp $@ ../bin/
	@if [ "$(shell uname)" = "Darwin" ]; then \
		echo "Signing binary for macOS..."; \
		codesign --force --deep --sign - ../bin/$@; \
	fi
	@echo "Build complete: ../bin/glutmolekel is ready."

# Rules for object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.C
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Cleanup
clean:
	rm -f *.o mui_source/*.o libmui.a ms glutmolekel

