# cc -fsanitize=address -I/opt/homebrew/opt/lapack/include -L/opt/homebrew/opt/lapack/lib -llapacke -llapack -lblas -std=c11 -O0 -g velest.c globals.c sub/*.c -lm -o velest_c
cc -I/opt/homebrew/opt/lapack/include -L/opt/homebrew/opt/lapack/lib -llapacke -llapack -lblas -std=c11 -O3  velest.c globals.c sub/*.c -lm -o velest_c
