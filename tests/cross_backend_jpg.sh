#!/usr/bin/env bash

# Copyright 2025
#
# This file is under BSD-3-Clause license.

set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
bundled="${OUTGUESS_BUNDLED:-}"
system="${OUTGUESS_SYSTEM:-}"
build_binaries="${CROSS_BACKEND_BUILD:-0}"

if [ -z "$bundled" ] || [ -z "$system" ]; then
	if [ "$build_binaries" != "1" ]; then
		echo "SKIP: set OUTGUESS_BUNDLED and OUTGUESS_SYSTEM or CROSS_BACKEND_BUILD=1" >&2
		exit 77
	fi

	build_dir="$root/tests/.cross-backend"
	mkdir -p "$build_dir"

	(
		cd "$root"
		./autogen.sh
		./configure --with-generic-jconfig
		make -C src
	)
	cp "$root/src/outguess" "$build_dir/outguess-bundled"
	bundled="$build_dir/outguess-bundled"

	(
		cd "$root"
		./configure --with-system-libjpeg
		make -C src
	)
	cp "$root/src/outguess" "$build_dir/outguess-system"
	system="$build_dir/outguess-system"
fi

if [ ! -x "$bundled" ] || [ ! -x "$system" ]; then
	echo "ERROR: outguess binaries not executable" >&2
	exit 1
fi

cd "$root/tests"

outfile_a="test-cross-bundled.jpg"
outfile_b="test-cross-system.jpg"
txt_a="text-cross-bundled.txt"
txt_b="text-cross-system.txt"

cleanup() {
	rm -f "$outfile_a" "$outfile_b" "$txt_a" "$txt_b"
}
trap cleanup EXIT

echo -e "\nEmbedding with bundled, extracting with system..."
"$bundled" -k "secret-key-001" -d message.txt test.jpg "$outfile_a"
"$system" -k "secret-key-001" -r "$outfile_a" "$txt_a"
grep "inside of the image" "$txt_a" || { echo ERROR; exit 1; }

echo -e "\nEmbedding with system, extracting with bundled..."
"$system" -k "secret-key-001" -d message.txt test.jpg "$outfile_b"
"$bundled" -k "secret-key-001" -r "$outfile_b" "$txt_b"
grep "inside of the image" "$txt_b" || { echo ERROR; exit 1; }
