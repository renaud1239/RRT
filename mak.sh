#!/bin/sh

echo '----------------------------------------'

echo "Compiling rrt-sse.c -> ./a-sse"
gcc -D__TEXTURE__ -O3 -Wno-attributes -mfpmath=sse -msse3 rrt-sse.c -lm -lpthread -o a-sse

echo '----------------------------------------'
echo "Happy crunching!"

