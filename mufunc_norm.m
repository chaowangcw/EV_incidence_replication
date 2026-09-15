function f = mufunc_norm(x2,theta2w)
% MUFUNC_NORM Construct heterogeneous utility in the random-coefficients model.
% INPUTS
%   x2      - [N x K] product characteristics with random coefficients.
%   theta2w - [K x (D+1)] parameter matrix. Column 1 contains coefficients
%             on random draws; columns 2:D+1 contain demographic interactions.
%
% OUTPUT
%   f or mu      - [N x ns] heterogeneous utility component by product and draw.
%
% REQUIRED GLOBAL STATE
%   ns    - number of simulation draws.
%   vfull - [N x (K*ns)] random-coefficient draws, stored in interleaved
%           draw-by-characteristic blocks expected by the indexing below.
%   dfull - [N x (D*ns)] log demographic/income draws in the corresponding
%           interleaved layout. EXP(dfull) enters the utility calculation.

global ns vfull dfull
[n k] = size(x2);
j = size(theta2w,2)-1;
mu = zeros(n,ns);
for i = 1:ns%
    	v_i = vfull(:,i:ns:k*ns);
       d_i = exp(dfull(:,i:ns:j*ns));
       mu(:,i) = (x2.*v_i*theta2w(:,1))+x2.*(d_i*theta2w(:,2:j+1)')*ones(k,1);
end
f = mu;



