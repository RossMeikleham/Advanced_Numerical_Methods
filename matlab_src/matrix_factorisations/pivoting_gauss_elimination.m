% Pivoting method of Gaussian Elimination of a Given Matrix
function a = naive_gauss_elimination(a) 
    m = size(a, 1);
    n = size(a, 2);

    for k=1:(m-1)
        for j=(k+1):m
            a(j,k) = a(j, k) / a(k, k);
            a(j,k+1:m) = a(j,k+1:m) - (a(j,k) * a(k,k+1:m));
        end
    end
end
