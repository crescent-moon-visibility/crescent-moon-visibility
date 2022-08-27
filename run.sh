#!/bin/bash

rm -f visibility.out
gcc -Wall -Werror -o visibility.out -Ofast \
    visibility.c thirdparty/astro_demo_common.c thirdparty/astronomy.c \
    || exit $?

rm -f *.png
time ./visibility.out 2022-08-27.png 2022-08-27T00:00:00Z
time ./visibility.out 2022-08-28.png 2022-08-28T00:00:00Z
time ./visibility.out 2022-08-29.png 2022-08-29T00:00:00Z
