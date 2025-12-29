# Molekel 4.3 for Modern macOS (ARM64)

This repository contains a ported version of Molekel 4.3, updated to compile on macOS Tahoe/Sonoma using Apple Silicon.

## Quick Start
1. Install dependencies: `brew install glui libtiff jpeg glfw gfortran`
2. `cd source`
3. Run `make`
4. Launch via `./bin/molekel`

## Log of modifications from the original source required to compile and run this legacy version of Molekel on modern macOS.
### browser.C (filesystem fix)

The original pwd() function used a manual inode-climbing method that fails on modern macOS security. 
It was replaced with the standard getcwd():

`void pwd(void) {
    if (getcwd(actual_directory, 1024) == NULL) {
        strncpy(actual_directory, ".", 128); 
    }
}`

### snap.C (Header Fixes)

Legacy Linux paths were updated to Homebrew paths:

Changed `#include <UnixImageIO/jpeglib.h>` to `#include <jpeglib.h>`

Added `#include <tiffio.h>`

Added `typedef uint32_t uint32;` to fix architecture-specific integer errors.

### libmui (The UI Library)
Because the provided libmui.a was built for Intel/SGI, the source files in the mui/ folder had to be recompiled for arm64.

### Compilation Commands (added/included in new Makefile)
Fortran Components

`gfortran -std=legacy -fno-range-check ms.f -o ms`

C++ Components

`g++ -I/opt/homebrew/include -Wno-deprecated-declarations -w -DLINUX -c *.C`

### Final Linking Command
Links against the macOS native OpenGL and GLUT frameworks instead of the old X11 libraries:

`g++ -Wno-deprecated-declarations -w -DLINUX -o glutmolekel \`
`    *.o \`
`    -L/opt/homebrew/lib \`
`    -lglui -ltiff -ljpeg -lglfw \`
`    -framework OpenGL \`
`    -framework GLUT \`
`    -framework Cocoa \`
`    ./libmui.a -lm`

## Note:
Menu Warnings: macOS will log "Internal inconsistency in menus." These are harmless and caused by the legacy GLUT implementation; they do not affect the program's functionality.

