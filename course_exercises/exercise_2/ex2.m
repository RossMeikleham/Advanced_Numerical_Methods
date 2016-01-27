a1 = [2, -1; 1, 3];
b1 = [3, 5];

a2 = [2, 3, -1; 1, 0, 1; -1, 2, -2];
b2 = [1, 0, 0];

a3 = [2, 1, 5; 5, 1, -5; 1, 2, 0];
b3 = [7, 9, 4];

a4 = [0, 6, 4; -2, 2, 2; 2, 1, -1];
b4 = [4, 2, 1];

tol = eps(10);

expected1 = a1 \ b1';
expected2 = a2 \ b2';
expected3 = a3 \ b3';
expected4 = a4 \ b4';

s1 = solve_naive_gauss(a1, b1);
s2 = solve_naive_gauss(a2, b2);
s3 = solve_naive_gauss(a3, b3);
%s4 = solve_naive_gauss(a4, b4);

assert(iseqtol(s1', expected1, tol));
assert(iseqtol(s2', expected2, tol));
assert(iseqtol(s3', expected3, tol));

s1p = solve_pivot_gauss(a1, b1);
s2p = solve_pivot_gauss(a2, b2);
s3p = solve_pivot_gauss(a3, b3);
s4p = solve_pivot_gauss(a4, b4);

assert(iseqtol(s1p', expected1, tol));
assert(iseqtol(s2p', expected2, tol));
assert(iseqtol(s3p', expected3, tol));
assert(iseqtol(s4p', expected4, tol));

cf = [2, -1, 0; -1, 2, -1; 0, -1, 2];
l = cholesky_factorisation(cf);
actual = chol(cf)';

assert(isequal(l, actual));
