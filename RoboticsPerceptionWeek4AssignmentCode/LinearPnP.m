function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 3) pose rotation

[n,ccc] = size(x);
A=zeros(2*n,12);

for i=1:n
    A(i*2-1,1:4)= zeros(1,4);
    A(i*2-1,5:8)= -1*[X(i,:) 1];
    A(i*2-1,9:12)= x(i,2)*[X(i,:) 1];
    
%     P(i*3-2,7)= x2(i,1);
%     P(i*3-2,8)= x2(i,2);
%     P(i*3-2,9)= 1;
%     P(i*3-2,10)= x2(i,1);
%     P(i*3-2,11= x2(i,2);
%     P(i*3-2,12)= 1;
    
    A(i*2,1:4) = [X(i,:) 1];
    A(i*2,5:8) = zeros(1,4);
    A(i*2,9:12)= -1*x(i,1)*[X(i,:) 1];
    
%     P(i*3-1,4)= x1(i,2)*x2(i,1);
%     P(i*3-1,5) x1(i,2)*x2(i,2);
%     P(i*3-1,6)= x1(i,2);
%     P(i*3-1,7)= x2(i,1);
%     P(i*3-1,8)= x2(i,2);
%     P(i*3-1,9)= 1;
%     P(i*3-1,10) = x2(i,1);
%     P(i*3-1,11) = x2(i,2);
%     P(i*3-1,12)= 1;
        
%     A(i*3,1:4)= -1*x(i,2)*[X(i,:) 1];
%     A(i*3,5:8)= x(i,1)*[X(i,:) 1];
%     A(i*3,9:12)= zeros(1,4);
%     P(i*3,4)= x1(i,2)*x2(i,1);
%     P(i*3,5) x1(i,2)*x2(i,2);
%     P(i*3,6)= x1(i,2);
%     P(i*3,7)= x2(i,1);
%     P(i*3,8)= x2(i,2);
%     P(i*3,9)= 1;
%     P(i*3,10)= x2(i,1);
%     P(i*3,11)= x2(i,2);
%     P(i*3,12)= 1;
end
[U,S,V]=svd(A);

p = V(:,12)/V(12,12);

P = [transpose(p(1:4,:)) ; transpose(p(5:8,:)) ; transpose(p(9:12,:))];

P_a = inv(K)*P;

R_b = P_a(1:3,1:3);

[U2,S2,V2] = svd(R_b);

if det(U2*(V2.'))>0
    R =  U2*transpose(V2);
    t = P_a(:,4)/S2(1,1);
    C = -1* R.' * t;
end

if det(U2*(V2.'))<0
    R = -1* U2* transpose(V2);
    t = -1* P_a(:,4)/S2(1,1);
    C = -1* R.' * t;
end

end
