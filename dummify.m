function f = dummify(var_to_dum)
% DUMMIFY Convert a categorical vector into a full set of indicator columns.
%
% PURPOSE
%   Construct dummy-variable matrices for firm, product, city, time, brand,
%   and body-style identifiers in the data-preparation stage of the main
%   replication script.
%
% INPUT
%   var_to_dum  - [N x 1] numeric vector or cell array of character vectors.
%
% OUTPUT
%   dummy_matrix - [N x K] full indicator matrix, where K is the number of
%                  levels returned by UNIQUE. No reference category is
%                  removed inside this function; the caller removes one when
%                  required for estimation.


temp=unique(var_to_dum);
dum=zeros(length(var_to_dum), length(temp));



if iscell(var_to_dum)==1  % cell array
    var_to_dum=upper(var_to_dum);
    for i=1:length(temp)
        temp1=find(strcmp(var_to_dum, temp(i)));
        dum(temp1,i)=1;
    end
    
else
    for i=1:length(temp)
        temp1=find(var_to_dum==temp(i));
        dum(temp1,i)=1;
    end
end

f=dum(:,1:length(temp));
clear temp temp0 temp1

