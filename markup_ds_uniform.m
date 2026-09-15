function f = markup_ds_uniform(theta2, delta, expmu, hyd, mid, mktsize) 
% PURPOSE: Calculate unified markup at model level using Bertrand competition model
% INPUTS:
%   - theta2: random coefficients vector
%   - delta: mean utility from BLP contraction mapping
%   - expmu: heterogeneous utility component
%   - hyd: half-year indicator (time dimension)
%   - mid: model ID (car model identifier)
%   - mktsize: market size (for weighting)
% OUTPUTS:
%   - f: unified markup for each product (same markup for same model across all cities)
% REQUIRED GLOBAL STATE
%   ns, cdid, fnmb, theta1, theti, thetj, op, dfull, and vfull must have been
%   initialized by the demand-estimation section of the main script.


global ns cdid fnmb theta1 theti thetj op dfull vfull

expmval = exp(delta); 
theta2w = full(sparse(theti,thetj,theta2));
shares = ind_sh(expmval,expmu);
s_j = mean(shares')';
nj = size(theta2w,2)-1;
sigmap = theta2w(1, 1);
betapi = theta2w(1, 2:nj+1)';


% ============================================================
% UNIFIED MARKUP CALCULATION - MODEL LEVEL AGGREGATION
% ============================================================
% This section aggregates market-level data to model level to compute
% unified markup for each car model across all cities.

% Initialize output variable
f = zeros(length(op),1);

% ============================================================
% MAIN LOOP: Iterate over each time period (half-year)
% ============================================================
for hy = min(hyd):max(hyd)
    
    % Initialize model-level variables for aggregation across citie
    hy_omega = zeros(max(unique(mid)));  % aggregated elasticity matrix at model level
    hy_s = zeros(max(unique(mid)),1);     % aggregated market share at model level
    
    % Get unique city IDs for the current time period
    cd_range = unique(cdid(find(hyd==hy)));
    
    % ============================================================
    % INNER LOOP: Iterate over each city market
    % ============================================================
    for t = min(cd_range):max(cd_range)
        trows = find(cdid == t); tempn = size(trows,1);
        
        % Check if this city has data for current time period
        if isempty(trows)
            continue;
        end
        
        % Get row indices for products in this city and time period
        time_trows = trows(hyd(trows) == hy);
        if isempty(time_trows)
            continue;
        end
        
        tempn = size(time_trows,1);
        
        vi = vfull(time_trows,1:ns); di = exp(dfull(time_trows,:)); 
        demo_i = zeros(tempn, ns);
        for n = 1: nj
            demo_i = di(:, (n-1)*ns + 1:n*ns)*betapi(n) + demo_i;
        end
        s_i = shares(time_trows,:);
        
        % ============================================================
        % Compute price elasticity matrix for this city
        % ============================================================
        temp = (s_i.*(vi*sigmap + theta1(1)+ demo_i))*s_i'/ns;
        temp1 = sum((s_i.*(vi*sigmap + theta1(1)+ demo_i))')'/ns;
        
        % Compute price elasticity matrix H
        H = (diag(temp1) - temp)./(op(trows)*ones(1,tempn));
        
        clear temp1 temp
        
        % ============================================================
        % Construct firm identity matrix
        % ============================================================
        omega_sta = zeros(tempn);
        f_i = fnmb(time_trows,:);
        for j = 1:tempn
            for h = 1:tempn
                omega_sta(j,h) = (f_i(j) == f_i(h));
            end
        end
        
        omega = H.*omega_sta;
        
        % ============================================================
        % Aggregate city-level data to model level
        % ============================================================

        modelid = mid(time_trows);
        city_mktsize = unique(mktsize(time_trows));
        if length(city_mktsize) ~= 1, error('Market size not constant within city'); end
        s_j_city = s_j(time_trows);
        
        for row = 1:length(time_trows)
            for col = 1:length(time_trows)
                % Aggregate elasticity matrix weighted by market size
                hy_omega(modelid(row),modelid(col)) = hy_omega(modelid(row),modelid(col)) + omega(row,col)*city_mktsize;
            end

            % Aggregate market share weighted by market size
            hy_s(modelid(row),1) = hy_s(modelid(row),1)+ s_j_city(row)*city_mktsize;
        end
        
        clear temprows omega omega_sta
    end
    
    % ============================================================
    % Solve Bertrand equilibrium for optimal markup
    % ============================================================
    nanproduct = find(hy_s==0);
    hy_s(nanproduct)=[];
    hy_omega(nanproduct,:) =[];
    hy_omega(:,nanproduct)=[];
    
    mkup = -inv(hy_omega)*hy_s;
    
    % Get list of valid model IDs
    temp = [1:max(unique(mid))]';
    temp(nanproduct)=[];
    
    % ============================================================
    % Assign unified markup to all products of the same model
    % ============================================================
    for index = 1:length(temp)
        f(find(mid==temp(index)&hyd==hy))= mkup(index);
    end
    
    clear tempn omega omega_sta H
end
