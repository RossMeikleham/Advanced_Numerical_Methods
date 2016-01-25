% Modified Gram Schmidt Process of Factorising a Given Matrix
% for numerical stability
function [q,r] = modified_gram_schmidt(a) 
    rows = size(a, 1);
    cols = size(a, 2);

    q = zeros(rows, cols); %q: mxn matrix
    r = zeros(cols, cols); %r: nxn matrix
    v = a; 
    
    for i = 1:cols
        r(i, i) = norm(v(:, i));
        q(:, i) = v(:, i) / r(i, i);
        for j=i+1:cols
            r(i, j) = q(:, i)' * v(:, j);
            v(:, j) = v(:, j) - r(i, j) * q(:, i);
        end
    end
end  


