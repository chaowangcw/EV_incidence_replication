function f = ind_sh(expmval,expmu)
% ins_sh Compute simulated individual product-choice probabilities in all markets.
%
% INPUTS
%   expmval - [N x 1] exponentiated product mean utilities.
%   expmu   - [N x ns] exponentiated heterogeneous utility components.
%
% OUTPUT
%   shares  - [N x ns] simulated individual product-choice probabilities.
%
% SOURCE
%   Based on Aviv Nevo's May 1998 random-coefficients demand code.

global ns cdindex cdid
eg = expmu.*kron(ones(1,ns),expmval);
for i = 1:length(cdindex)
    if i == 1
        sum1(i,:) = sum(eg(1:cdindex(i),:));
    else
        sum1(i,:) = sum(eg(cdindex(i-1) + 1: cdindex(i),:)); 
    end
end

denom1 = 1./(1+sum1);
denom = denom1(cdid,:);
f = eg.*denom;
