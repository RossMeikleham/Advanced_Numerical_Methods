% Naive method of Gaussian Elimination of a Given Matrix
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

function a = back_substitution_gauss(a)
    m = size(a, 1);
    n = size(a, 2);   
 
    % Row operations on result
    for k=1:m-1
        for i=k+1:m
            a(i, n) = a(i, n) - (a(i, k) * a(k, n));
        end
    end

    % Back substitution
	a(m, n) = a(m, n) / a(m, m);
	for i=m-1:-1:1
		s = a(i, n);
		for j=i+1:m
			s = s - (a(i, j) * a(j, n));
		end
		a(i, n) = s / a(i, i);
	end 
end
