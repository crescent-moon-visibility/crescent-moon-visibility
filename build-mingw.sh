#!/bin/bash

rm -f visibility.exe
for i in 1 2 4; do
    x86_64-w64-mingw32-gcc -static-libgcc -Wl,-Bstatic -fno-exceptions -DPIXEL_PER_DEGREE=$i \
        -Wall -Werror -o visibility-$i.exe -fopenmp -O3 visibility.cc thirdparty/astronomy.c
    strip visibility$i.exe
done
