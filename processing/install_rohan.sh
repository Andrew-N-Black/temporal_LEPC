#!/bin/bash
# =============================================================================
# One-time build of ROHan (Renaud et al. 2019) from source.
#
# ROHan works directly from a BAM/CRAM + reference FASTA
# =============================================================================
set -euo pipefail

PROJ=/scratch/gautschi/blackan/GROUSE/old_vs_new
ROHAN_DIR=$PROJ/tools/ROHan
GSL_PREFIX=$PROJ/tools/gsl

mkdir -p "$(dirname "$ROHAN_DIR")"

module --force purge
module load gcc/14.1.0

# GSL is a real link-time dependency (unlike bzip2 below, this can't just be
# disabled) and Gautschi has no 'gsl' module at all (confirmed). Build it
# from source into the project directory if it isn't there already, then
# point the compiler/linker at it via CPATH/LIBRARY_PATH/LD_LIBRARY_PATH --
# this lets ROHan's own Makefile find it with a plain `-lgsl` exactly as it
# would a system-installed copy, with no need to edit ROHan's build files.
if [ ! -f "$GSL_PREFIX/lib/libgsl.a" ] && [ ! -f "$GSL_PREFIX/lib/libgsl.so" ]; then
    echo "Building GSL from source into $GSL_PREFIX..."
    GSL_BUILD_DIR=$(mktemp -d)
    ( cd "$GSL_BUILD_DIR" \
      && wget -q https://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz \
      && tar xzf gsl-latest.tar.gz \
      && cd gsl-*/ \
      && ./configure --prefix="$GSL_PREFIX" \
      && make -j4 \
      && make install )
    rm -rf "$GSL_BUILD_DIR"
    if [ ! -f "$GSL_PREFIX/lib/libgsl.a" ] && [ ! -f "$GSL_PREFIX/lib/libgsl.so" ]; then
        echo "ERROR: GSL build did not produce libgsl in $GSL_PREFIX/lib" >&2
        echo "  Check the configure/make output above for errors." >&2
        exit 1
    fi
    echo "GSL built successfully at $GSL_PREFIX"
else
    echo "GSL already built at $GSL_PREFIX -- skipping."
fi

export CPATH="$GSL_PREFIX/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="$GSL_PREFIX/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="$GSL_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# rohan (the binary, once built) will also need LD_LIBRARY_PATH set to find
# libgsl.so at runtime, not just at build time -- make sure it's set in
# run_rohan.sh too (see the note printed at the end of this script).

# htslib (bundled by ROHan as a submodule) needs libbzip2 dev headers to
# build CRAM support. Gautschi has no 'bzip2' module (confirmed), so don't
# bother attempting it -- the build-with-retry logic below falls back to
# disabling bzip2 support in htslib automatically if needed.

ROHAN_BIN="$ROHAN_DIR/bin/rohan"

if [ -x "$ROHAN_BIN" ]; then
    echo "ROHan already built at $ROHAN_BIN -- skipping."
    echo "Delete $ROHAN_DIR first if you want to rebuild from scratch."
    exit 0
fi

if [ -d "$ROHAN_DIR" ]; then
    # A directory here but no binary just means a previous build attempt
    # didn't finish (e.g. the bzip2/GSL issues from earlier) -- reuse the
    # existing clone and retry rather than treating this as fatal. Only
    # bail out if it doesn't actually look like a ROHan checkout at all.
    if [ ! -f "$ROHAN_DIR/Makefile" ]; then
        echo "ERROR: $ROHAN_DIR exists but doesn't look like a ROHan checkout" >&2
        echo "  (no top-level Makefile found). Remove it and rerun, or" >&2
        echo "  investigate what's actually there." >&2
        exit 1
    fi
    echo "$ROHAN_DIR already exists from a previous attempt -- reusing it and retrying the build."
else
    echo "Cloning ROHan (with submodules -- it bundles its own htslib copy)..."
    git clone --recursive https://github.com/grenaud/ROHan.git "$ROHAN_DIR"
fi

cd "$ROHAN_DIR"
echo "Building..."
BUILD_LOG=$(mktemp)
if ! make 2>&1 | tee "$BUILD_LOG"; then
    if grep -qE "libbzip2 development files not found|cannot find -lbz2" "$BUILD_LOG" \
       && [ ! -x "$ROHAN_BIN" ]; then
        # Only worth patching+retrying if rohan itself is what's missing.
        # ROHan's build also fetches a separate companion tool, bam2prof
        # (ancient-DNA damage profiling -- not used here since this is
        # modern high-coverage Illumina data), as part of the same `make`
        # invocation. If IT fails on -lbz2/-lcurl after rohan has already
        # linked successfully, that's fine to ignore entirely -- see the
        # final check below, which looks for rohan specifically rather
        # than trusting make's overall exit status.
        echo ""
        echo ">>> Build failed on a bzip2 dependency, and no working bzip2 dev"
        echo ">>> headers/library exist on this system. Removing bzip2 from"
        echo ">>> ROHan's own build (both htslib's CRAM support and ROHan's"
        echo ">>> final link step reference it) and retrying. This only"
        echo ">>> affects CRAM files whose blocks specifically use bzip2"
        echo ">>> compression, which is uncommon -- most GATK/sarek-produced"
        echo ">>> CRAMs don't use it."

        # Patch 1: htslib's own configure -- don't require bzlib.h to build
        # (idempotent: only add --disable-bz2 if it isn't already there).
        if grep -q -- "--disable-libcurl" src/Makefile && ! grep -q -- "--disable-bz2" src/Makefile; then
            sed -i 's/--disable-libcurl/--disable-libcurl --disable-bz2/' src/Makefile
        fi

        # Patch 2: ROHan's own final link command hardcodes -lbz2 separately
        # from htslib's configuration -- strip it. Safe to run even if
        # already stripped (no-op), so no existence guard needed here.
        sed -i 's/-lbz2 //g; s/-lbz2//g' src/Makefile

        # htslib's own leftover build state from the failed attempt needs to
        # go, or the retry will skip straight back to the cached configure
        # failure instead of re-running configure with the new flag.
        rm -rf lib/htslib
        make 2>&1 | tee "$BUILD_LOG" || true
    fi
fi
rm -f "$BUILD_LOG"

# Check specifically for the rohan binary, NOT make's overall exit status --
# make also builds the separate bam2prof companion tool (ancient-DNA damage
# profiling) as part of the same invocation, and its failure is irrelevant
# to running rohan itself on modern samples.
if [ ! -x "$ROHAN_BIN" ]; then
    echo "ERROR: build finished but $ROHAN_BIN is not there or not executable." >&2
    echo "  Check the make output above for errors -- a missing GSL/htslib" >&2
    echo "  dev package is the most common cause." >&2
    exit 1
fi

echo ""
echo "Build succeeded: $ROHAN_BIN"
echo ""
echo "(If you saw errors above about bam2prof/-lcurl/-lbz2: that's a separate"
echo " companion tool for ancient-DNA damage profiling, not needed for"
echo " modern samples, and safe to ignore -- rohan itself is what matters"
echo " and it built fine, per the check above.)"
echo ""
echo "IMPORTANT: rohan was linked against a custom GSL build, not a system"
echo "one, so it needs LD_LIBRARY_PATH set to find libgsl.so at RUNTIME too,"
echo "not just here at build time. run_rohan.sh has this added already --"
echo "if you run rohan any other way, remember:"
echo "  export LD_LIBRARY_PATH=$GSL_PREFIX/lib:\$LD_LIBRARY_PATH"
echo ""
echo "Before running the batch script, confirm the exact flag names with:"
echo "  export LD_LIBRARY_PATH=$GSL_PREFIX/lib:\$LD_LIBRARY_PATH"
echo "  $ROHAN_BIN --help"
echo "(flag names in run_rohan.sh were written from ROHan's documented usage,"
echo " but it's worth a 10-second sanity check against your actual build)."
