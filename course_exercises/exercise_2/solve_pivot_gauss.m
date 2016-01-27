function x = solve_pivot_gauss(a, b)
    
    [u, l, p] = pivoting_gauss_elimination(a);
    b = p * b'; 
    b = forward_substitution_gauss(l, b);
    x = back_substitution_gauss(u, b)';
end
