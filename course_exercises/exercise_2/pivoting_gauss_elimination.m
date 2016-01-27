% Pivoting method of Gaussian Elimination of a Given Matrix
function [a, l, p] = pivoting_gauss_elimination(a) 
    rows = size(a, 1);
    cols= size(a, 2);
    p = eye(rows);
    l = eye(rows);

    for i=1:cols-1
        %Locate the row which contains the largest magnitude
        %value for the current column, and swap it with the
        %topmost row which hasn't been considered 
        [_, r] = max(arrayfun(@abs, a(i:rows,i))); 
        r = r + (i - 1);
        a([i r], :) = a([r i], :); 
        l([i r], 1:i-1) = l([r i], 1:i-1);
        p([i r], :) = p([r i], :);       
        
        %Subtract from first row, should gives 0s
        %in the current column below the given row
        for j=(i+1):rows
            l(j, i) = a(j, i) / a(i, i);
            a(j, i:cols) = a(j, i:cols) - (l(j, i) * a(i, i:cols));
        end
    end
end
