% Pivoting method of Gaussian Elimination of a Given Matrix
function a = pivoting_gauss_elimination(a) 
    rows = size(a, 1);
    cols= size(a, 2);

    for i=1:(cols-1) 
        %Locate the row which contains the largest magnitude
        %value for the current column, and swap it with the
        %topmost row which hasn't been considered 
        [_, r] = max(arrayfun(@abs, a(i:rows,i))); 
        r = r + (i - 1);
        temp = a(r, :);
        a(r, :) = a(i, :);
        a(i, :) = temp;

        %Subtract from first row, should gives 0s
        %in the current column below the given row
        for j=(i+1):rows
            a(j, :) = a(j,:) - ((a(j, i) * a(i, :)) / a(i, i));
        end
    end
end

pivoting_gauss_elimination([4, 2, 5; 5, 6, 1])
pivoting_gauss_elimination([3, 2, -4, 3; 2, 3, 3, 15; 5, -3, 1, 14])
pivoting_gauss_elimination([2, 1, 1, 0, 1; 4, 3, 3, 1, 2; 8, 7, 9, 5, 3; 6, 7, 9, 8, 4])
