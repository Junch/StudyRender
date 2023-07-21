cmake -B_build64 -H. -A x64
cmake --build _build64
::ctest --test-dir _build64 -C Debug