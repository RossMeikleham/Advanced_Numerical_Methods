function b = back_substitution_gauss(a, b)
    [m, n] = size(a); 
 
    % Row operations on result
    for k=1:m-1
        for i=k+1:m
            b(i) = b(i) - (a(i, k) * b(k));
        end
    end   
 
    % Back substitution
	b(m) = b(m) / a(m, m);
	for i=m-1:-1:1
		s = b(i);
		for j=i+1:m
			s = s - (a(i, j) * b(j));
		end
		b(i) = s / a(i, i);
	end 
end
