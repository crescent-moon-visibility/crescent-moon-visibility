#!/bin/bash

CC=gcc
OPEN=xdg-open
if [[ "$OSTYPE" == "darwin"* ]]; then
    CC="/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/llvm/lib -fno-rtti" # -mllvm -polly
    OPEN=open
fi

rm -f visibility.out
$CC -fopenmp -O3 -Wall -Werror -o visibility.out -fno-exceptions -DPIXEL_PER_DEGREE=4 visibility.cc thirdparty/astronomy.c -lm || exit $?

echo "Compiliation is completed, now let's run the code."

DATE=${DATE:-2022-08-27}
TYPE=${TYPE:-evening}
METHOD=${METHOD:-yallop}
time ./visibility.out $DATE map $TYPE $METHOD $DATE.png || (echo Not successful && exit 1)

FONT_ARGS=()
if [[ "$OSTYPE" == "darwin"* && -f /System/Library/Fonts/Supplemental/Arial.ttf ]]; then
    FONT_ARGS=(-font /System/Library/Fonts/Supplemental/Arial.ttf)
elif [[ -f /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf ]]; then
    FONT_ARGS=(-font /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf)
fi

magick composite -blend 60 "$DATE.png" map.png "$DATE.png"
TYPE="$(tr '[:lower:]' '[:upper:]' <<< ${TYPE:0:1})${TYPE:1}"
METHOD="$(tr '[:lower:]' '[:upper:]' <<< ${METHOD:0:1})${METHOD:1}"
magick "$DATE.png" "${FONT_ARGS[@]}" -pointsize 20 -fill black -draw "gravity south text 0,0 '$TYPE, $METHOD, $DATE'" "$DATE.png"
$OPEN $DATE.png
