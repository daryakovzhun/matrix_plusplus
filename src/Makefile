CC             =    g++
CFLAGS         =    -Wall -Werror -Wextra
GCOV_FLAGS     =    -fprofile-arcs -ftest-coverage -lgcov
CHECK_FLAGS    =    -lgtest -lgmock -pthread
ALL_FLAGS      =    $(CFLAGS) $(GCOV_FLAGS) $(CHECK_FLAGS)

UNAME := $(shell uname)

SRCSTESTS      =    $(wildcard tests/*.cpp)
OBJSTESTS      =    $(SRCSTESTS:%.cpp=objs_test/%.o)
SRCS           =    $(wildcard source/*.cpp)
OBJS           =    $(SRCS:%.cpp=objs_test/%.o)
ALL            =     $(OBJS) $(OBJSTESTS)

OBJTESTDIR 	   = 	objs_test/tests objs_test/source
NAME           =    test
REPORT_NAME    = 	report

objs_test/%.o: %.cpp
	@mkdir -p $(OBJTESTDIR)
	$(CC) -o $@ -c $< $(ALL_FLAGS)

all: s21_matrix_oop.a

test: $(ALL) 
	$(CC) $(ALL) -o $(NAME) $(ALL_FLAGS)
	./$(NAME)

clean:
	rm -rf $(OBJS) *.o *.gcno *.gcda *.info *.a $(NAME) check_res finish objs_test

s21_matrix_oop.a: $(OBJS)
	ar rc s21_matrix_oop.a $(OBJS)
	ranlib s21_matrix_oop.a

gcov_report: test
	lcov -t "$(REPORT_NAME)" -o objs_test/$(REPORT_NAME).info -c -d .
	genhtml objs_test/$(REPORT_NAME).info -o finish

open_html:
ifeq ($(UNAME), Linux)
	firefox finish/index.html
endif
ifeq ($(UNAME), Darwin)
	open finish/index.html
endif


.PHONY: all test gcov_report clean s21_matrix_oop.a open_html