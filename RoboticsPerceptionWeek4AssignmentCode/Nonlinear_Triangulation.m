function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix for both
%     cameras
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 


[n,m] = size(X0);
X = zeros(n,3);
for i =1:n
    re = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), transpose(X0(i,:)));  
    re2 = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), re); 
    re3 = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1(i,:), x2(i,:), x3(i,:), re2); 
    X(i,:) = transpose(re3);
end

end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
    J = [Jacobian_Triangulation(C1, R1, K, X0);
         Jacobian_Triangulation(C2, R2, K, X0);
         Jacobian_Triangulation(C3, R3, K, X0);];
   
    x_1 = K*R1*(X0-C1);
    x_2 = K*R2*(X0-C2);
    x_3 = K*R3*(X0-C3);
   
    f = [x_1(1)/x_1(3); x_1(2)/x_1(3); x_2(1)/x_2(3); x_2(2)/x_2(3); x_3(1)/x_3(3); x_3(2)/x_3(3) ];
    b = [x1(1,1)        ; x1(1,2)        ; x2(1,1)        ; x2(1,2)        ; x3(1,1)        ; x3(1,2)         ];
    X= X0 + inv(J.' * J ) * (J.') * (b-f);
    
end
    
function J = Jacobian_Triangulation(C, R, K, X)
    x = K*R*(X-C);
    f=K(1,1);
    p_x = K(1,3);
    p_y = K(2,3);
    
    dev_w_X = R(3,:); 
    dev_u_X = [f*R(1,1)+p_x*R(3,1) f*R(1,2)+p_x*R(3,2) f*R(1,3)+p_x*R(3,3)];
    dev_v_X = [f*R(2,1)+p_y*R(3,1) f*R(2,2)+p_y*R(3,2) f*R(2,3)+p_x*R(3,3)];
    J = [ (x(3)*dev_u_X-x(1)*dev_w_X)/(x(3)*x(3)) ; (x(3)*dev_v_X-x(2)*dev_w_X)/(x(3)*x(3)) ];
end
