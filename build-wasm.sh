#!/bin/bash

emcc -DPIXEL_PER_DEGREE=1 -o visibility.wasm -O3 -fno-exceptions -fno-rtti -flto -o index.wasm \
    visibility.cc thirdparty/astronomy.c
