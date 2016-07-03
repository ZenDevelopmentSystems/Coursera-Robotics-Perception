function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

[n,m] = size(x1);
A = zeros(n,9);

for i=1:n
    A(i,1)= x1(i,1)*x2(i,1);
    A(i,2)= x1(i,1)*x2(i,2);
    A(i,3)= x1(i,1);
    A(i,4)= x1(i,2)*x2(i,1);
    A(i,5)= x1(i,2)*x2(i,2);
    A(i,6)= x1(i,2);
    A(i,7)= x2(i,1);
    A(i,8)= x2(i,2);
    A(i,9)= 1;
end

[U, S, V] = svd(A);

h = V(:,9);

H = [h(1) h(2) h(3);h(4) h(5) h(6);h(7) h(8) h(9)];

[U2,S2,V2] = svd(H);

S2(3,3)=0;

F_b = U2 * S2 * V2.';

F = F_b/norm(F_b);

end
