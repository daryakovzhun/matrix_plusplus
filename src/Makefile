CC             =    g++
CFLAGS         =    -Wall -Werror -Wextra
GCOV_FLAGS     =    -fprofile-arcs -ftest-coverage -lgcov
CHECK_FLAGS    =    -lgtest -lgmock -pthread
ALL_FLAGS      =    $(CFLAGS) $(GCOV_FLAGS) $(CHECK_FLAGS)

SRCS           =    $(wildcard s21_matrix_oop.cpp)
OBJS           =    $(SRCS:%.cpp=objs_test/%.o)
OBJTESTDIR 	   = 	objs_test/ 
NAME           =    objs_test/tests
REPORT_NAME    = 	report

objs_test/%.o: %.cpp
	@mkdir -p $(OBJTESTDIR)
	$(CC) -MMD -o $@ -c $< $(ALL_FLAGS)

all: s21_matrix_oop.a test

test: $(OBJS)
# $(CC) $(CFLAGS) $(CHECK_FLAGS) s21_matrix_test.c $(OBJS) -o $(NAME)
	$(CC) $(OBJS) tests/s21_matrix_test.cpp -o $(NAME) $(ALL_FLAGS)
	./objs_test/tests

clean:
	rm -rf $(OBJS) *.o *.gcno *.gcda *.info *.a $(NAME) check_res finish

s21_matrix_oop.a: $(OBJS)
	ar rc s21_matrix_oop.a $(OBJS)
	ranlib s21_matrix_oop.a

gcov_report: test
# $(CC) -o $(NAME) $(SRCS) s21_matrix_test.c $(ALL_FLAGS) 
	./$(NAME)
	lcov -t "$(REPORT_NAME)" -o objs_test/$(REPORT_NAME).info -c -d .
	genhtml objs_test/$(REPORT_NAME).info -o finish

# geninfo *.gcda *.gcno -b . -o check_res 
# genhtml check_res -o finish

.PHONY: all test gcov_report clean s21_matrix_oop.a