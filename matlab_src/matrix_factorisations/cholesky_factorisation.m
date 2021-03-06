% Pivoting method of Gaussian Elimination of a Given Matrix
function l = cholesky_factorisation(a) 
    rows = size(a, 1);
    cols= size(a, 2);
    l = zeros(size(a));

    for i=1:cols
        for j=i:rows
            if (i == j) 
                s = a(i, i);
                for k=1:(i-1)
                    s = s - (l(i, k) ^ 2);
                end
                l(i, i) = sqrt(s);
             
            % i != j
            else
                s = a(j, i);
                for k=1:(i-1)
                    s = s - (l(i, k) * l(j, k));
                end
                l(j, i) = s / l(i, i);
            end
        end
    end
end

cholesky_factorisation([25, 15, -5; 15, 18, 0; -5, 0, 11])
%expected 5 0 0; 3 3 0; -1 1 3;
