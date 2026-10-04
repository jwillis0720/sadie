#!/bin/bash
# Build IgBLAST 1.22.0 igblastn + makeblastdb natively for arm64 macOS with Apple's toolchain only
# (Xcode or Command Line Tools; no Homebrew/conda). Static NCBI libs, no SRA/VDB, macOS 11+.
# Usage: build_igblast.sh DEST_DIR [WORK_DIR]   (igblast-arm64.patch must sit next to this script;
# WORK_DIR defaults to /tmp/igblast-1.22.0-build, so the ~500 MB build tree stays out of the repo)
set -euo pipefail
[ "$(uname -sm)" = "Darwin arm64" ] || { echo "$0: needs arm64 macOS, found $(uname -sm)" >&2; exit 1; }
# Re-run once in a clean environment so no conda/Homebrew tools or inherited CFLAGS leak in.
[ -n "${IGB_CLEAN:-}" ] || exec /usr/bin/env -i IGB_CLEAN=1 HOME="$HOME" PATH=/usr/bin:/bin:/usr/sbin:/sbin /bin/bash "$0" "$@"

ver=1.22.0
sha256=a04eb195087447cf42e673401ae2d9f6ef8414b34f1534724a6953f41224df68 # NCBI posts no .md5 for IgBLAST; same as bioconda
here=$(cd "$(dirname "$0")" && pwd)
dest=$(mkdir -p "${1:?usage: $0 DEST_DIR [WORK_DIR]}" && cd "$1" && pwd)
work=$(mkdir -p "${2:-/tmp/igblast-$ver-build}" && cd "${2:-/tmp/igblast-$ver-build}" && pwd -P)
export CC=/usr/bin/clang CXX=/usr/bin/clang++ MACOSX_DEPLOYMENT_TARGET=11.0
# Drop the build directory from __FILE__ strings so the binaries embed no local paths. configure writes physical
# paths (hence pwd -P); CPPFLAGS keeps its default -O2, which a user CFLAGS/CXXFLAGS would turn into -O.
export CPPFLAGS="-fmacro-prefix-map=$work/="
# NCBI's default fast flags minus -fomit-frame-pointer, which NCBI already drops for AppleClang 15 on arm
# (breaks C++ exceptions) but not for newer AppleClang.
export DEF_FAST_FLAGS="-O3 -finline-functions -fstrict-aliasing"
# pipe2 is a macOS 27 API: found on a newer build host, it gets weak-linked and is NULL on macOS 11-26.
export ac_cv_func_pipe2=no

cd "$work"
tgz=ncbi-igblast-$ver-src.tar.gz
[ -f "$tgz" ] || curl -fsSLO "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/$ver/$tgz"
echo "$sha256  $tgz" | shasum -a 256 -c -
rm -rf "ncbi-igblast-$ver-src" && tar xzf "$tgz"
cd "ncbi-igblast-$ver-src"
patch -p1 < "$here/igblast-arm64.patch"
cd c++
# configure writes `hostname` into the build signature that -version-full prints; give it a neutral name.
mkdir -p "$work/shim" && printf '#!/bin/sh\necho sadie-build\n' > "$work/shim/hostname" && chmod +x "$work/shim/hostname"
echo "configure + make (log: $work/build.log)"
{ PATH="$work/shim:$PATH" ./configure.orig --with-build-root=ReleaseMT --without-debug --with-mt --without-dll \
    --with-bin-release --without-openmp --without-vdb --without-python --with-flat-makefile --without-autodep \
    --without-makefile-auto-update --without-caution &&
  cd ReleaseMT/build &&
  make -j"$(sysctl -n hw.ncpu)" -f Makefile.flat igblastn.exe makeblastdb.exe; } > "$work/build.log" 2>&1 ||
  { tail -40 "$work/build.log"; exit 1; }

for b in igblastn makeblastdb; do
  cp "../bin/$b" "$dest/$b"
  /usr/bin/strip "$dest/$b"
  /usr/bin/codesign --verify "$dest/$b" 2>/dev/null || /usr/bin/codesign -s - -f "$dest/$b"
  /usr/bin/codesign --verify --strict "$dest/$b"
  if /usr/bin/otool -L "$dest/$b" | tail -n +2 | grep -vE '^[[:space:]]+/(usr/lib|System/Library)/'; then
    echo "$b links a library outside /usr/lib and /System" >&2; exit 1
  fi
  if /usr/bin/nm -mu "$dest/$b" | grep -F ') weak external '; then
    echo "$b weak-imports symbols newer than macOS 11 (NULL on older systems)" >&2; exit 1
  fi
  "$dest/$b" -version
done
