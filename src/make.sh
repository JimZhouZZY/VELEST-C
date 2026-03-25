cc -fsanitize=address -I/opt/homebrew/opt/lapack/include -L/opt/homebrew/opt/lapack/lib -llapacke -llapack -lblas -std=c11 -O0 -g vel_globals.c velest.c sub/*.c -lm -o velest_c
