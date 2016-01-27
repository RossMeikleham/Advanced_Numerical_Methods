function b = back_substitution_gauss(a, b)
    [m, n] = size(a); 
    
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
