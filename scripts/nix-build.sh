#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"

configure_flags=()
if [ -n "${CONFIGURE_FLAGS-}" ]; then
	read -r -a configure_flags <<< "${CONFIGURE_FLAGS}"
fi
run_tests="${RUN_TESTS:-1}"

run_build() {
	./autogen.sh
	./configure "${configure_flags[@]}"
	make -C src
	if [ "$run_tests" != "0" ]; then
		make -C tests check
	fi
}

if [ -n "${IN_NIX_SHELL-}" ]; then
	run_build
else
	nix develop -c bash -c "cd '$root' && $(declare -f run_build); run_build"
fi
