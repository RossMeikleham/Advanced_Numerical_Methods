% Gram Schmidt Process of Factorising a Given Matrix
function [q,r] = gram_schmidt(a) 
    rows = size(a, 1);
    cols = size(a, 2);

    q = zeros(rows, cols); %q: mxn matrix
    r = zeros(cols, cols); %r: nxn matrix
    
    for n=1:cols
        v = a(:, n);
        % Calc R-Values for n-th column
        for i=1:n-1
            r(i, n) = q(:, i)' * a(:, n);
            v = v - r(i, n) * q(:, i); 
        end
  
        r(n, n) = norm(v);

        % Calc Q-Values for n-th column
        q(:, n) = v / r(n, n);
    end
end  

[q, r] = gram_schmidt([1 1 0; 1 0 1; 0 1 1]);
q
r

q_expected = [1/sqrt(2) 1/sqrt(6) -1/sqrt(3); 
              1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 
              0 2/sqrt(6) 1/sqrt(3)];

r_expected = [sqrt(2), 1/sqrt(2), 1/sqrt(2);
              0, 3/sqrt(6), 1/sqrt(6);
              0, 0, 2/sqrt(3)];

tol = eps(1)
q - q_expected < tol
r - r_expected < tol

