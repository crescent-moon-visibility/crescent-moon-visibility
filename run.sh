#!/bin/bash

CC="gcc -O3"
if [[ "$OSTYPE" == "darwin"* ]]; then
    CC="/opt/homebrew/opt/llvm/bin/clang -O3 -L/opt/homebrew/opt/llvm/lib -fno-rtti" # -mllvm -polly
fi

rm -f visibility.out
$CC -fopenmp -Wall -Werror -o visibility.out -fno-exceptions \
    visibility.cc thirdparty/astro_demo_common.c thirdparty/astronomy.c -lm \
    || exit $?

#rm -f *.png
# time ./visibility.out 2022-06-29.png 2022-06-29T00:00:00Z
# composite -blend 80 2022-06-29.png map.png 2022-06-29.png
#time ./visibility.out 2022-08-28.png 2022-08-28T00:00:00Z
#time ./visibility.out 2022-08-29.png 2022-08-29T00:00:00Z

DATE=2022-09-25
time ./visibility.out ${DATE}T00:00:00Z evening yallop map $DATE.png || (echo Not successful && exit 1)
composite -blend 60 $DATE.png map.png $DATE.png
#open $DATE.png
