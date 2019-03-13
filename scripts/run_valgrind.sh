#!/usr/bin/env bash
valgrind --xml=yes --xml-file=val.xml --suppressions=scripts/valgrind.supp  --leak-check=full --leak-check=yes $@
