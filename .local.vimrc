set makeprg=./build.sh
set guioptions+=!

let g:CommandTWildIgnore=&wildignore . ",*.o,*.so,*.a,build_*,.lib*"

