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

[q, r] = modified_gram_schmidt([1 1 0; 1 0 1; 0 1 1])

q_expected = [1/sqrt(2) 1/sqrt(6) -1/sqrt(3); 
              1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 
              0 2/sqrt(6) 1/sqrt(3)];

r_expected = [sqrt(2), 1/sqrt(2), 1/sqrt(2);
              0, 3/sqrt(6), 1/sqrt(6);
              0, 0, 2/sqrt(3)];

tol = eps(1000)
q - q_expected < tol
r - r_expected < tol

