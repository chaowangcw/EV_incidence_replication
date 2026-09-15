function f = gmmobj_micro_2(theta2)
% GMMOBJ_MICRO_2 Initial-weight GMM objective with macro and micro moments.
%
% INPUT
%   theta2    - free nonlinear demand parameters.
%
% OUTPUT
%   f / objective - scalar GMM criterion.
%
% REQUIRED GLOBAL STATE
%   x1, op        - linear characteristics and observed product prices.
%   IV, invA      - instruments and initial aggregate-moment weight.
%   ns, dfull     - simulation count and log-income draws.
%   inc_r, s_inc  - micro-sample row indices and observed incomes.
%   inc_IV        - instruments for the income micro moments.
%   GLOBAL_SCALE_FACTOR - scale conversion applied to observed incomes.
%   theta1, niter - updated here for later use and progress reporting.
%
% SOURCE
%   Based on Aviv Nevo's May 1998 code and modified to add income micro moments.

global theta1 x1 IV  op niter invA ns dfull inc_r s_inc inc_IV GLOBAL_SCALE_FACTOR

tic
niter=niter+1;

[delta, expmu] = meanval_micro_PP_norm(theta2);

% the following deals with cases were the min algorithm drifts into region where the objective is not defined
if max(isnan(delta)) == 1
    f = 1e+10;
else

    y = [delta]; %
    X1 = [log(op) x1];

    temp1 = X1'*IV;
    temp2 = IV'*y;
    theta1 = inv(temp1*invA*temp1')*(temp1*invA*temp2);
    %save theta1 theta1;
    clear temp1 temp2
    gmmresid = y - X1*theta1;

    theta1(1:9)
    temp1 = gmmresid'*IV;

    %%%% micro moments %%%%
    mktshr_ind = ind_sh(exp(delta), expmu);
    avginc = sum(exp(dfull).*mktshr_ind./(sum(mktshr_ind,2)*ones(1,ns)),2);
    inc_resid = (s_inc*GLOBAL_SCALE_FACTOR - avginc(inc_r,:));
    temp2 = inc_resid'*inc_IV;

    %% Moments
    g= [temp1,temp2];
    Wt = blkdiag(invA, inv([inc_IV'*inc_IV]));
    %% objective function
    f=g*Wt*g';
    
    save gmmresid  gmmresid inc_resid;

end

disp(['GMM objective:  ' num2str(f)])
if floor(niter/10)==niter/10
    display(['Parameter estimates']);
    num2str([theta2])
end

toc
