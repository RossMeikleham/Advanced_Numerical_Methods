function x = solve_naive_gauss(a, b)
    a = naive_gauss_elimination(a);
    b = forward_substitution_gauss(a, b);   
    x = back_substitution_gauss(a, b);
end
