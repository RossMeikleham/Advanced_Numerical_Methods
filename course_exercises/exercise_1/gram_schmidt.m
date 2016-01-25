% Gram Schmidt Process of Factorising a Given Matrix
function [q,r] = gram_schmidt(a) 
    rows = size(a, 1);
    cols = size(a, 2);

    q = zeros(rows, cols); %q: mxn matrix
    r = zeros(cols, cols); %r: nxn matrix
    
    for n=1:cols
        v = a(:, n);
        % Calc r-Values for n-th column
        for i=1:n-1
            r(i, n) = q(:, i)' * a(:, n);
            v = v - (r(i, n) * q(:, i)); 
        end
  
        r(n, n) = norm(v);
        % Calc q-Values for n-th column
        q(:, n) = v / r(n, n);
    end 
    
    
