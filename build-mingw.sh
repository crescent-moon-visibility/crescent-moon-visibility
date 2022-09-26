#!/bin/bash

rm -f visibility.exe
x86_64-w64-mingw32-gcc -static-libgcc -Wl,-Bstatic -fno-exceptions -DPIXEL_PER_DEGREE=4 \
    -Wall -Werror -o visibility.exe -fopenmp -O3 visibility.cc thirdparty/astronomy.c
strip visibility.exe
