% Pivoting method of Gaussian Elimination of a Given Matrix
function a = pivoting_gauss_elimination(a) 
    rows = size(a, 1);
    cols= size(a, 2);

    for i=1:cols 
        %Locate the row which contains the largest magnitude
        %value for the current column, and swap it with the
        %topmost row which hasn't been considered 
        [_, r] = max(arrayfun(@abs, a(i:rows,i))); 
        r = r + (i - 1);
        a([i r], :) = a([r i], :);        

        %Subtract from first row, should gives 0s
        %in the current column below the given row
        for j=(i+1):rows
            a(j, :) = a(j,:) - ((a(j, i) * a(i, :)) / a(i, i));
        end
    end
end
