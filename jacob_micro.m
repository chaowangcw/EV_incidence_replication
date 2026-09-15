function f = jacob_micro(mval,theta2)
% JACOB_MICRO Compute the mean-utility Jacobian for nonlinear parameters.
%
% INPUTS
%   mval   - [N x 1] exponentiated mean utilities.
%   theta2 - [P x 1] free nonlinear demand parameters.
%
% OUTPUT
%   f      - [N x P] Jacobian d(delta)/d(theta2).
%
% REQUIRED GLOBAL STATE
%   ns, theti, thetj - simulation count and theta2 index mapping.
%   cdid, cdindex    - market membership and cumulative market endpoints.
%   x2, dfull        - random-coefficient characteristics and log-income draws.
%
% FILE INPUT
%   Loads ps2_<ns>.mat, which must contain the simulation draws used in the
%   demand estimation, including v. The file name depends on ns.
%
% SOURCE
%   Modified from Aviv Nevo's random-coefficients demand code.

global ns theti thetj cdid cdindex x2 dfull
s=['load ps2_' num2str(ns)];
eval(s)
theta2w = full(sparse(theti,thetj,theta2));
expmu = exp(mufunc_norm(x2,theta2w));
shares = ind_sh(mval,expmu);

clear expmu

[n,K] = size(x2);
J = size(theta2w,2) - 1;
f1 = zeros(size(cdid,1),K*(J + 1));

% computing (partial share)/(partial sigma)
for i = 1:K
	xv = (x2(:,i)*ones(1,ns)).*v(cdid,ns*(i-1)+1:ns*i); 
    temp = xv.*shares;
    for ii = 1:length(cdindex)
        if ii == 1
            sum1(ii,:) = sum(temp(1:cdindex(ii),:));
        else
            sum1(ii,:) = sum(temp(cdindex(ii-1) + 1: cdindex(ii),:)); 
        end
    end

	f1(:,i) = mean((shares.*(xv-sum1(cdid,:)))')';
	clear xv temp sum1
end

% If no demogr comment out the next para
% computing (partial share)/(partial pi)
for j = 1:J
d = exp(dfull(cdid,ns*(j-1)+1:ns*j));
	temp1 = zeros(size(cdid,1),K);
	for i = 1:K
		xd=(x2(:,i)*ones(1,ns)).*d;  
        temp = xd.*shares;
        for ii = 1:length(cdindex)
            if ii == 1
                sum1(ii,:) = sum(temp(1:cdindex(ii),:));
            else
                sum1(ii,:) = sum(temp(cdindex(ii-1) + 1: cdindex(ii),:)); 
            end
        end
		temp1(:,i) = mean((shares.*(xd-sum1(cdid,:)))')';
		clear xd temp sum1
	end
	f1(:,K*j+1:K*(j+1)) = temp1;
	clear temp1
end

rel = theti + (thetj - 1) * max(theti) ;

% computing (partial delta)/(partial theta2)

f = zeros(size(cdid,1),size(rel,1));
n = 1;
for i = 1:size(cdindex,1)
	temp = shares(n:cdindex(i),:);
	H1 = temp*temp';
	H = (diag(sum(temp')) - H1)/ns;
	f(n:cdindex(i),:) = - inv(H)*f1(n:cdindex(i),rel);
	n = cdindex(i) + 1;
end
