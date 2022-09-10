#!/bin/bash

CC=gcc
if [[ "$OSTYPE" == "darwin"* ]]; then
    CC="/usr/local/opt/llvm/bin/clang -L/usr/local/opt/llvm/lib" # -mllvm -polly
fi

rm -f visibility.out
$CC -fopenmp -Wall -Werror -o visibility.out -Ofastest -fno-exceptions -fno-rtti \
    visibility.cc thirdparty/astro_demo_common.c thirdparty/astronomy.c -lm \
    || exit $?

#rm -f *.png
# time ./visibility.out 2022-06-29.png 2022-06-29T00:00:00Z
# composite -blend 80 2022-06-29.png map.png 2022-06-29.png
#time ./visibility.out 2022-08-28.png 2022-08-28T00:00:00Z
#time ./visibility.out 2022-08-29.png 2022-08-29T00:00:00Z

time ./visibility.out 2020-05-23.png 2020-05-23T00:00:00Z
composite -blend 60 2020-05-23.png map.png 2020-05-23.png
