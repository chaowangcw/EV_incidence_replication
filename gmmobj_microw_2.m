function f = gmmobj_microw_2(theta2)
% GMMOBJ_MICROW_2 Optimal-weight GMM objective with macro and micro moments.
%
% INPUT
%   theta2    - free nonlinear demand parameters.
%
% OUTPUT
%   f - scalar GMM criterion.
%
% REQUIRED GLOBAL STATE
%   x1, op        - linear characteristics and observed product prices.
%   IV, NewW      - aggregate instruments and estimated moment covariance.
%   inc_IV, inc_c_W - micro instruments and estimated moment covariance.
%   ns, dfull     - simulation count and log-income draws.
%   inc_r, s_inc  - micro-sample row indices and observed incomes.
%   GLOBAL_SCALE_FACTOR - scale conversion applied to observed incomes.
%   theta1, niter - updated here for later use and progress reporting.
%
% SOURCE
%   Modified from Aviv Nevo's code to implement optimal weighting with income
%   micro moments.

global theta1 x1 IV op  niter ns NewW dfull inc_r s_inc NewincW 
global inc_c_W inc_b_W brand_matrix city_matrix inc_IV GLOBAL_SCALE_FACTOR nobs nobs_c
tic
niter=niter+1;

[delta, expmu] = meanval_micro_PP_norm(theta2);


% the following deals with cases were the min algorithm drifts into region
% where the objective is not defined
if max(isnan(delta)) == 1
    f = 1e+10;
else
    
    
    y = [delta]; %
    X1 = [log(op) x1]; %,xs;

    iv2=IV;
    temp1 = X1'*iv2;
    temp2 = iv2'*y;
    theta1 = inv(temp1*inv(NewW)*temp1')*(temp1*inv(NewW)*temp2);
    clear temp1 temp2
    gmmresid = y - X1*theta1;
    
    save theta1 theta1 y  gmmresid;
    
    [theta1(1:10)]
    temp1 = gmmresid'*iv2;
    mktshr_ind = ind_sh(exp(delta), expmu);
    avginc = sum(exp(dfull).*mktshr_ind./(sum(mktshr_ind,2)*ones(1,ns)),2);
    inc_resid = (s_inc*GLOBAL_SCALE_FACTOR - avginc(inc_r,:));
    temp2 = inc_resid'*inc_IV;

    %% Moments
    g= [temp1, temp2];
    Wt = blkdiag(inv(NewW), inv(inc_c_W));
    f=g*Wt*g';
    
    save incresid_sep  gmmresid

end

disp(['GMM objective:  ' num2str(f)])
if floor(niter/10)==niter/10
    display(['Parameter estimates']);
    num2str([theta2])
end


toc

