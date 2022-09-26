#!/bin/bash

rm -f visibility.exe
x86_64-w64-mingw32-gcc -static-libgcc -Wl,-Bstatic -fno-exceptions -DPIXEL_PER_DEGREE=2 \
    -Wall -Werror -o visibility.exe -O3 visibility.cc thirdparty/astronomy.c
x86_64-w64-mingw32-strip visibility.exe
