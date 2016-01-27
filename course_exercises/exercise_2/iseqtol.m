function t = iseqtol(a, b, tol)
    tolMat = repmat(tol, size(a));
    t = isequal(arrayfun(@abs, a - b) < tolMat, ones(size(a)));

end
