#!/bin/bash

CC=gcc
if [[ "$OSTYPE" == "darwin"* ]]; then
    CC="/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/llvm/lib -fno-rtti" # -mllvm -polly
fi

rm -f visibility.out
$CC -fopenmp -O3 -Wall -Werror -o visibility.out -fno-exceptions visibility.cc thirdparty/astronomy.c -lm || exit $?

#rm -f *.png
# time ./visibility.out 2022-06-29.png 2022-06-29
# composite -blend 80 2022-06-29.png map.png 2022-06-29.png
#time ./visibility.out 2022-08-28.png 2022-08-28
#time ./visibility.out 2022-08-29.png 2022-08-29

DATE=2022-09-25
time ./visibility.out $DATE map evening yallop $DATE.png || (echo Not successful && exit 1)
composite -blend 60 $DATE.png map.png $DATE.png
#open $DATE.png
