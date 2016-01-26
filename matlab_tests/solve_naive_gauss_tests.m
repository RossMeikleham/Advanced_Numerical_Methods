addpath('../matlab_src/matrix_factorisations')


%% Solve Naive Gauss 2x2
expected = [-8; 3];
actual = solve_naive_gauss([1,5,7;-2,-7,-5]);
assert(isequal(expected, actual));



