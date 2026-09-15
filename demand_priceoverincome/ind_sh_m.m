function f = ind_sh_m(expmval,expmu, nprod_m, ns)
% IND_SH_M Compute simulated individual choice probabilities in one market.
%
% INPUTS
%   expmval - [J x 1] exponentiated product mean utilities in the market.
%   expmu   - [J x ns] exponentiated heterogeneous utility components.
%   nprod_m - scalar product count J for this market.
%   ns      - scalar number of simulation draws.
%
% OUTPUT
%   f  - [J x ns] simulated individual product-choice probabilities.
%
% REPRODUCIBILITY AND SOURCE
%   Uses no globals, files, or random draws. Based on Aviv Nevo's May 1998
%   random-coefficients demand code.

eg = expmu.*(expmval*ones(1,ns));

sum1=ones(nprod_m,1)*sum(eg);
f = eg./(1+sum1);
