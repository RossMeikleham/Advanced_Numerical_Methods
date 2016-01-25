
function exercise1
format long 

q_expected = [1/sqrt(2) 1/sqrt(6) -1/sqrt(3); 
              1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 
              0 2/sqrt(6) 1/sqrt(3)];

r_expected = [sqrt(2), 1/sqrt(2), 1/sqrt(2);
              0, 3/sqrt(6), 1/sqrt(6);
              0, 0, 2/sqrt(3)];

disp('Testing Naive Gram Schmidt:');
[q1, r1] = gram_schmidt([1 1 0; 1 0 1; 0 1 1]);
q1
r1


tol = eps(1000)
q1 - q_expected < tol
r1 - r_expected < tol

disp('Testing Modified Gram Schmidt:');
[q2, r2] = modified_gram_schmidt([1 1 0; 1 0 1; 0 1 1]);
q2
r2

q2 - q_expected < tol
r2 - r_expected < tol

t5 = [0 -1 2; 1 0 2; 1 -1 0];
disp('Original Gram Shmidt T5:');
gram_schmidt(t5)
disp('Modified Gram Schmidt T5:');
modified_gram_schmidt(t5)

% Calculate Sigma
n = 80;
s =  zeros(n, n);
for j = 1:n
    s(j, j) = 2^(-j);
end

[u, r] = qr(randn(n, n));
[v, r] = qr(randn(n, n));
a = u * s * v'; 

disp('C5 Gram Schmidt:');
[q1, r1] = gram_schmidt(a);

disp('C5 Modified Gram Schmidt:');
[q2, r2] = modified_gram_schmidt(a);
 
semilogy(1:80, diag(r1), 'ro');
hold on
semilogy(1:80, diag(r2), 'bx');

end



