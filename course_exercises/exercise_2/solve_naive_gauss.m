function x = solve_naive_gauss(a, b)
    a = naive_gauss_elimination(a);
    x = back_substitution_gauss(a, b);
end
