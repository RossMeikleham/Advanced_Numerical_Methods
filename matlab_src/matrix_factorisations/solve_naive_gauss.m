function x = solve_naive_gauss(a)
    a = naive_gauss_elimination(a);
    a = back_substitution_gauss(a);
    x = a(:,size(a, 2));
end
