set makeprg=./build.sh
set guioptions+=!

let g:ale_linters_explicit = 1

let b:ale_linters = ["gcc", "cppcheck", "flawfinder"]
let b:ale_fixers = ['clangtidy', 'clang-format']

