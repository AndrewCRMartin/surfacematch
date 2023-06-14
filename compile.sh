compile="/usr/bin/cc -I${HOME}/include -L${HOME}/lib -ansi -pedantic -Wall"


$compile -o surface surface.c -lbiop -lgen -lm -lxml2
$compile -o match match.c -lbiop -lgen -lm -lxml2

