#!/bin/bash

rm -f visibility.exe
x86_64-w64-mingw32-gcc -static-libgcc -Wl,-Bstatic -fno-exceptions \
    -Wall -Werror -o visibility.exe -Ofast \
    visibility.cc thirdparty/astro_demo_common.c thirdparty/astronomy.c
x86_64-w64-mingw32-strip visibility.exe
