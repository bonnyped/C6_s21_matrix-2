CC         = gcc 
CFLAGS     = -fno-stack-protector -g -Wall -Werror -Wextra -std=c11
TARGET     = s21_matrix.a
CLEAN      = rm -rf
CLANG      = clang-format --style=google
TO_CLANG   = *.c *.h
SRC        = $(wildcard *.c)
OBJ        = $(patsubst %.c, %.o, $(SRC))
TEST_FLAGS = $(shell pkg-config --cflags --libs check)
OS        := $(shell uname -s)

PREF_TEST  = ./tests/
PREF_GCOV  = ./gcov/

all: clean $(TARGET) test clang gcov_report open

$(TARGET) : $(OBJ)
	ar rc $(TARGET) *.o
	$(CLEAN) *.o

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(PREF_TEST)test.o $(TARGET)
	$(CC) $(PREF_TEST)test.o $(TARGET) -o $(PREF_TEST)test -L. $(TEST_FLAGS)
	$(CLEAN) $(PREF_TEST)*.o
	$(PREF_TEST)test

$(PREF_TEST)test.o: $(PREF_TEST)test.c
	$(CC) -c $(PREF_TEST)test.c -o $(PREF_TEST)test.o

$(PREF_TEST)test.c:
	$(CLEAN) $(PREF_TEST)test $(PREF_TEST)*.c
	checkmk clean_mode=1 $(PREF_TEST)include.file $(PREF_TEST)*.check > $(PREF_TEST)test.c
	$(CLANG) $(PREF_TEST)$(TO_CLANG) -i

gcov_report: clean $(PREF_TEST)test.c
	mkdir $(PREF_GCOV)
	@$(CC) -fprofile-arcs -ftest-coverage $(PREF_TEST)test.c $(SRC) -o $(PREF_GCOV)gcov_report $(TEST_FLAGS)
	@$(PREF_GCOV)gcov_report
	@lcov -t "test" -o $(PREF_GCOV)test.info -c -d .
	@genhtml -o report $(PREF_GCOV)test.info
	$(CLEAN) ./*.gcda  ./*.gcno $(PREF_GCOV)*.gcda  $(PREF_GCOV)*.gcno

open:
ifeq ($(OS),Darwin)
	open -a "Google Chrome" report/index.html
else
	xdg-open report/index.html
endif

leak: test
ifeq ($(OS),Darwin)
	leaks -atExit -- $(PREF_TEST)test
else
	valgrind --track-origins=yes --leak-check=full  --log-file=leak.log $(PREF_TEST)test 
endif 

clang:
	$(CLANG) $(TO_CLANG) -n
	$(CLANG) $(TO_CLANG) -i

clean:
	$(CLEAN) *.o *.out *.a $(PREF_TEST)test.o $(PREF_TEST)test.c report *.log $(PREF_TEST)*.info $(PREF_TEST)test $(PREF_GCOV)