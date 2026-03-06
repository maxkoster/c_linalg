linalgmake: main.c lib/linalg.c
	gcc -o linalg main.c lib/linalg.c -I. -Wall