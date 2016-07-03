function E = EssentialMatrixFromFundamentalMatrix(F,K)
%% EssentialMatrixFromFundamentalMatrix
% Use the camera calibration matrix to esimate the Essential matrix
% Inputs:
%     K - size (3 x 3) camera calibration (intrinsics) matrix for both
%     cameras
%     F - size (3 x 3) fundamental matrix from EstimateFundamentalMatrix
% Outputs:
%     E - size (3 x 3) Essential matrix with singular values (1,1,0)

E_b = K.' * F * K;

[U,S,V] = svd(E_b);

S = S/S(1,1);

S(3,3) = 0;

S(2,2) = 1;

E = U * S * V.';

end

