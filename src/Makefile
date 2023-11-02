all: clean test

.PHONY: test
test:
	g++ -std=c++17 -Wall -Werror -Wextra -pthread -fPIC matrix_oop_tests.cc matrix_oop.cc -o test -lgtest -lgtest_main -lm
	./test

.PHONY: matrix_oop.a
matrix_oop.a : matrix_oop.o
	ar rc libmatrix_oop.a *.o
	ranlib libmatrix_oop.a
	cp libmatrix_oop.a matrix_oop.a

.PHONY: matrix_oop.o
matrix_oop.o:
	g++ -std=c++17 -Wall -Werror -Wextra -c matrix_oop.cc

clean:
	rm -rf *.o *.out *.gch *.dSYM *.gcov *.gcda *.gcno *.a matrix_oop_tests *.css *.html vgcore* report *.info *.gz *.log test