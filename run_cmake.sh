#!/bin/bash

if [ ! -d "build" ]; then
  mkdir build
fi

pushd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Ninja" ./..
# cmake -DCMAKE_BUILD_TYPE=Debug -G "Ninja" ./..
popd

mv ./build/compile_commands.json ./compile_commands.json
