function f = var_cov_micro_2(theta2)
% VAR_COV_MICRO_2 Compute the joint GMM covariance matrix of demand estimates.
%
% INPUT
%   theta2 - second-stage nonlinear demand estimates.
%
% OUTPUT
%   f / vcov  - [(K+1+P) x (K+1+P)] covariance matrix, where K=size(x1,2)
%               and the additional linear parameter is the price coefficient.
%
% REQUIRED GLOBAL STATE
%   NewW, inc_c_W - estimated covariance matrices for macro and micro moments.
%   IV, inc_IV    - aggregate and micro instruments.
%   x1, x2, op    - demand characteristics and observed price in levels.
%   theti, thetj, ns, dfull - parameter mapping and simulation data.
%   inc_r, s_inc, GLOBAL_SCALE_FACTOR - income micro-moment inputs.
%
% FILE INPUTS
%   Loads mvalold from mvalold_norm.mat. Loads theta1.mat and relies on that
%   file containing both theta1 and gmmresid, as written by gmmobj_microw_2.
%
% SOURCE
%   Modified from Aviv Nevo's May 1998 covariance calculation to include
%   income micro moments.

global NewW IV x1 x2 op theti thetj ns dfull inc_r s_inc
global inc_c_W inc_b_W brand_matrix city_matrix inc_IV GLOBAL_SCALE_FACTOR nobs nobs_c
load mvalold_norm
load theta1

%% (partial mom2)/(partial theta)
[N, K] = size(x1);
Z = size(IV,2);
theta2w = full(sparse(theti,thetj,theta2));
expmu = exp(mufunc_norm(x2,theta2w));
mktshr_ind = ind_sh(mvalold, expmu);
avginc = sum(exp(dfull).*mktshr_ind./(sum(mktshr_ind,2)*ones(1,ns)),2);
inc_resid0 = (s_inc*GLOBAL_SCALE_FACTOR - avginc(inc_r,:));
ee = eye(K+1);
nth2 = length(theta2);
ee2 = eye(nth2);
tao_c = zeros(size(inc_IV,1), K+1+nth2);
for j = 1: K+ 1 + nth2
    if j <=K + 1
        tolx = 10^(-5);
        theta1nd = theta1.*(ones(K+1,1)+ee(:,j)*tolx);
        deltacoef = theta1(j)*tolx;
        dmval = exp([op x1]*theta1nd);
        mktshr_ind = ind_sh(dmval, expmu);
        avginc = sum(exp(dfull).*mktshr_ind./(sum(mktshr_ind,2)*ones(1,ns)),2);
        inc_resid = (s_inc*GLOBAL_SCALE_FACTOR - avginc(inc_r,:));

    else
        tolx = 10^(-5);
        theta2n=theta2.*(ones(nth2,1) + ee2(:,j-K - 1)*tolx);
        theta2w = full(sparse(theti,thetj,theta2n));
        deltacoef = theta2(j-K-1)*tolx;
        dmu = exp(mufunc_norm(x2,theta2w));
        mktshr_ind = ind_sh(mvalold, dmu);
        avginc = sum(exp(dfull).*mktshr_ind./(sum(mktshr_ind,2)*ones(1,ns)),2);
        inc_resid = (s_inc*GLOBAL_SCALE_FACTOR - avginc(inc_r,:));
    end

    tao_c(:,j) = (inc_resid - inc_resid0)/deltacoef;

end
N = size(x1,1);
Z = size(IV,2);
temp = jacob_micro(mvalold,theta2);
a =  [[op x1 temp]; tao_c]'*blkdiag(IV, inc_IV);
IVres = IV.*(gmmresid*ones(1,Z));
shatc = inc_IV.*(inc_resid0*ones(1,cols(inc_IV)));
b = blkdiag(IVres'*IVres, shatc'*shatc);
winvA = blkdiag(inv(NewW), inv(inc_c_W));

f = inv(a*winvA*a')*a*winvA*b*winvA*a'*inv(a*winvA*a');

