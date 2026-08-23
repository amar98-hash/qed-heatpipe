cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release

cmake --build build --target cppqed_smoke_test -j$(nproc)


if the submodule cppqed is mistakenly cleaned and updated, run this patch 
git -C qed-submodule/cppqed apply ../../patches/cppqed-gitrevision-cmake-fix.patch