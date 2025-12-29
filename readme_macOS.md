# Molekel 4.3 for Modern macOS (ARM64)

This repository contains a ported version of Molekel 4.3, updated to compile on macOS Tahoe/Sonoma using Apple Silicon.

## Quick Start
1. Install dependencies: `brew install glui libtiff jpeg glfw gfortran`
2. `cd source`
3. Run `make`
4. Launch via `./bin/molekel`

This document outlines the steps and modifications required to compile and run this legacy version of Molekel on modern macOS.
2. Key Code Modifications
browser.C (The Filesystem Fix)

The original pwd() function used a manual inode-climbing method that fails on modern macOS security. It was replaced with the standard getcwd():

C++
void pwd(void) {
    if (getcwd(actual_directory, 1024) == NULL) {
        strncpy(actual_directory, ".", 128); 
    }
}
snap.C (Header Fixes)

Legacy Linux paths were updated to Homebrew paths:

Changed #include <UnixImageIO/jpeglib.h> to #include <jpeglib.h>

Added #include <tiffio.h>

Added typedef uint32_t uint32; to fix architecture-specific integer errors.

libmui (The UI Library)

Because the provided libmui.a was built for Intel/SGI, the source files in the mui/ folder had to be recompiled for arm64:

Bash
gcc -c -I/opt/homebrew/include -Wno-deprecated-declarations -w mui/*.c
ar rc libmui.a *.o
ranlib libmui.a
3. Compilation Commands
Fortran Components

To compile the math routines (like ms.f), use the legacy flag to ignore 1970s-style syntax errors:

Bash
gfortran -std=legacy -fno-range-check ms.f -o ms
C++ Components

Compile individual .C files with the Homebrew include path:

Bash
g++ -I/opt/homebrew/include -Wno-deprecated-declarations -w -DLINUX -c *.C
4. Final Linking Command
This command links the object files against the macOS native OpenGL and GLUT frameworks instead of the old X11 libraries:

Bash
g++ -Wno-deprecated-declarations -w -DLINUX -o glutmolekel \
    *.o \
    -L/opt/homebrew/lib \
    -lglui -ltiff -ljpeg -lglfw \
    -framework OpenGL \
    -framework GLUT \
    -framework Cocoa \
    ./libmui.a -lm
5. Troubleshooting Runtime Issues
Instant Exit (Status 2): Ensure you are running the binary from the directory containing its configuration and data folders.

Menu Warnings: macOS will log "Internal inconsistency in menus." These are harmless and caused by the legacy GLUT implementation; they do not affect the program's functionality.

