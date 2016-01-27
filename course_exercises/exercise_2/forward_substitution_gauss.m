 function b = forward_substitution_gauss(a, b)
    [m, n] = size(a);
    for k=1:m-1
        for i=k+1:m
            b(i) = b(i) - (a(i, k) * b(k));
        end
    end
 
end
