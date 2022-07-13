GCC = gcc -Wall -Werror -Wextra -g
SOURCE = s21_*.c
OBJ = $(SOURCE:.c=.o)
GCOVFLAGS=-L. --coverage
HTML=lcov -t "test" -o rep.info -c -d .
OS = $(shell uname)

ifeq ($(OS), Darwin)
	CHECK_H=-lcheck
else
	CHECK_H=-lcheck -pthread -lrt -lm
endif

all: gcov_report

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt s21_decimal *.cfg *a.out *.dSYM

test: s21_decimal.a
	$(GCC) decimal_test.c $(SOURCE) -o test $(CHECK_H) $(GCOVFLAGS)
	./test

s21_decimal.a: clean $(OBJ)
	ar rcs s21_decimal.a $(OBJ)

gcov_report: test
	$(HTML)
	genhtml -o report rep.info

check:
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config $(SOURCE) *.h
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extensions=c $(SOURCE) *.h
	rm -rf CPPLINT.cfg
	make test
ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif