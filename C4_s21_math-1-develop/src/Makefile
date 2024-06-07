.PHONY: all test s21_math.a clang clean cppcheck rebuild
WWW = -Wall -Wextra -Werror -std=c11
OS=$(shell uname -s)
LINUX=-lsubunit -lrt
LIBS=-lcheck -lm -lpthread

ifeq ($(OS), Linux)
	OPEN = xdg-open
else
	OPEN = open
endif

all:s21_math.a

test: clean_test s21_math.a
ifeq ($(OS), Darwin)
	gcc $(WWW) -fprofile-arcs -ftest-coverage -pedantic -std=c11 -g s21_math.c tests/*.c s21_math.a -o unit_test $(LIBS) 
else
	gcc  $(WWW) -fprofile-arcs -ftest-coverage -pedantic -std=c11 -g s21_math.c tests/*.c s21_math.a -o unit_test $(LIBS) $(LINUX)
endif
	./unit_test

gcov_report: clean_test s21_math.a
ifeq ($(OS), Darwin)
	gcc $(WWW) -fprofile-arcs -ftest-coverage -pedantic -std=c11 -g s21_math.c tests/*.c s21_math.a -o unit_test $(LIBS) 
else
	gcc  $(WWW) -fprofile-arcs -ftest-coverage -pedantic -std=c11 -g s21_math.c tests/*.c s21_math.a -o unit_test $(LIBS) $(LINUX)
endif
	./unit_test
	lcov -t "unit_test" -o test.info -c -d .
	genhtml -o report test.info
	$(OPEN) report/index.html

s21_math.a:
	gcc -c s21_*.c
	ar rcs s21_math.a s21_*.o
	ranlib s21_math.a

clang:
	clang-format -i *.c

clean_obj:
	rm -rf *.o

clean_lib:
	rm -rf *.a

clean_test:
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -rf test
	rm -rf report
	rm -rf *.dSYM

clean: clean_lib clean_lib clean_test clean_obj
	rm -rf unit_test

rebuild: clean all
