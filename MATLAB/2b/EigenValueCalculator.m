function [] = EigenValueCalculator()
%This function is designed to calculate the eigenvalues of a 2x2 matrix

m1 = 100/1000000000;
m2 = 50/1000000000;
k = 1;

A = [((2*k)/m1) (-k/m1); (-k/m2) ((2*k)/m2)];

[v, d] = eig(A);

d1 = d(1,1)
d2 = d(2,2)



end

