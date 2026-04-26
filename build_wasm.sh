#!/bin/bash
# Build the WASM module using Emscripten
#
# Prerequisites:
#   1. Install Emscripten: https://emscripten.org/docs/getting_started/downloads.html
#   2. Activate: source emsdk_env.sh
#
# Usage:
#   ./build_wasm.sh

set -e

echo "Building WASM module..."

emcc -O3 \
    -s WASM=1 \
    -s EXPORTED_FUNCTIONS='["_init_table","_query_pi","_get_table_end","_get_checkpoint_interval","_get_num_checkpoints","_malloc","_free"]' \
    -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap","HEAPU8"]' \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=16777216 \
    -s ENVIRONMENT=web \
    -s MODULARIZE=1 \
    -s EXPORT_NAME='PiModule' \
    -s NO_EXIT_RUNTIME=1 \
    --no-entry \
    -o docs/pi_query.js \
    pi_query_wasm.cpp

echo "Done! Files created:"
echo "  docs/pi_query.js    (JS glue)"
echo "  docs/pi_query.wasm  (WASM binary)"
ls -lh docs/pi_query.js docs/pi_query.wasm
