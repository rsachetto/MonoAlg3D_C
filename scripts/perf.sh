#!/bin/bash
perf record -F 99 -g -- $@
perf script | ../FlameGraph/stackcollapse-perf.pl > out.perf-folded
../FlameGraph/flamegraph.pl out.perf-folded > perf-kernel.svg

