function [mu, mval, pd_i, p_e, sj_e, mkup_e, subsidy_total, rg_w, sub_ind, cv] = optsubsidy2(x)
% OPTSUBSIDY2 Recover equilibrium and welfare outcomes for a subsidy rule.
%
% PURPOSE
%   After OPTSB2/FMINCON selects subsidy coefficients, resolve the Bertrand
%   equilibrium and construct the policy outcomes reported in the replication:
%   prices, shares, subsidy payments, range, and compensating variation.
%
% INPUT
%   x - [3 x 1] coefficients in max(x(1)+x(2)*income+x(3)*range,0). Units
%       inherit the income, range, and price scaling used by the main script.
%
% OUTPUTS
%   mu            - heterogeneous utility component.
%   mval          - utility component at equilibrium prices.
%   pd_i          - consumer prices after individual subsidies.
%   p_e           - counterfactual equilibrium producer prices.
%   sj_e          - counterfactual product shares.
%   mkup_e        - model-level markups assigned to products.
%   subsidy_total - scalar aggregate subsidy expenditure.
%   rg_w          - scalar share-weighted range for domestic EVs.
%   sub_ind       - flattened city-by-draw average subsidy vector.
%   cv            - scalar market-size-weighted compensating variation; its
%                   monetary units follow the price normalization.
%
% REQUIRED GLOBAL STATE
%   Uses demand estimates, baseline utilities (including eg0), simulation
%   draws, prices and costs, market sizes, eligibility indicators, ownership
%   matrices, and market IDs initialized by the main script.
%
% FILE INPUT
%   Loads opt_subsidy.mat. It must contain tol, temp, theta2w, theta1, nj,
%   subsidy_0, sub_tt0, p_s_0, hy, cd_range, and model_variant.

global range_tt  vi di tts qcost x1e x2e ehat mktsizets ev_tt import_tt eg0
global  omega_sta cost_sim p_s_o local_tt sc_tt  k ns alpha_i subdec subextra_id cdid

load opt_subsidy
rw_ev = find(ev_tt == 1 & import_tt == 0);
norm = 1;
avgnorm = 1;
i = 0;
p_e = p_s_o;

hy_temp = hy(temp);
cdid_hy = cdid(temp);
mid = model_variant(temp);
mktsizets_hy = mktsizets;


while norm > 1e-1*tol*10^(floor(i/50)) & avgnorm > 1e-3*tol*10^(floor(i/50))
    p_d = (p_e).*tts + qcost;
    if isreal(log(p_d)) == 1
        mval = zeros(size(temp,1),ns);
        for ri = 1:ns
            v_i = vi(:,ri:ns:k*ns);
            d_i = di(:,ri:ns:nj*ns);

            pd_i(:,ri) = p_d - subextra_id(:,ri).*max((x(1) + x(2)*d_i + x(3)*range_tt), 0);
            mval(:,ri) = [log(pd_i(:,ri)) x1e]*theta1 + ehat;
            mu(:,ri) = ([log(pd_i(:,ri)),x2e(:,2:end)].*v_i*theta2w(:,1))...
                +[log(pd_i(:,ri)),x2e(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
        end
        eg = exp(mu + mval);
        
        % Calculate market share for each city separately
        s_i = zeros(size(eg));
        
        hy_given = 8;
        hy_omega = zeros(max(unique(mid)));  
        hy_s = zeros(max(unique(mid)),1);
        
        for t = min(cd_range):max(cd_range)
            trows = find(cdid_hy == t); 
            time_trows = trows(hy_temp(trows) == hy_given);
            
            % Calculate market share for this city
            s_i(time_trows,:) = eg(time_trows,:)./(ones(length(time_trows),1)*(1 + sum(eg(time_trows,:))));
            
            % Calculate omega matrix and aggregate
            temp0 = (s_i(time_trows,:).*alpha_i(time_trows,:)./pd_i(time_trows,:).*tts(time_trows))*s_i(time_trows,:)'/ns;
            temp1 = sum((s_i(time_trows,:).*alpha_i(time_trows,:)./pd_i(time_trows,:).*tts(time_trows))')'/ns;
            H = (diag(temp1) - temp0);
            
            omega = H.*omega_sta{t};
            
            modelid = mid(time_trows);
            city_mktsize = unique(mktsizets_hy(time_trows));
            if length(city_mktsize) ~= 1, error('Market size not constant within city'); end
            s_j_city = mean(s_i(time_trows,:), 2);
            
            for row = 1:length(time_trows)
                for col = 1:length(time_trows)
                    hy_omega(modelid(row),modelid(col)) = hy_omega(modelid(row),modelid(col)) + omega(row,col)*city_mktsize;
                end
                hy_s(modelid(row),1) = hy_s(modelid(row),1)+ s_j_city(row)*city_mktsize;
            end
        end
        sj_e = mean(s_i,2);
        tot_ded = mean(subextra_id.*max((x(1) + x(2)*di + x(3)*range_tt), 0).*s_i,2);
        
        nanproduct = find(hy_s==0);
        hy_s(nanproduct)=[];
        hy_omega(nanproduct,:) =[];
        hy_omega(:,nanproduct)=[];
        
        mkup_model = -inv(hy_omega)*hy_s;
        
        temp_models = [1:max(unique(mid))]';
        temp_models(nanproduct)=[];
        
        mkup_e = zeros(length(p_e),1);
        for index = 1:length(temp_models)
            mkup_e(find(mid==temp_models(index)))= mkup_model(index);
        end
        
        pnew = mkup_e + cost_sim;
        norm = max(abs(pnew - p_e));
        avgnorm = mean(abs(pnew - p_e));
        p_e = pnew;
    end
    i = i + 1;
end
disp(['# of iterations for price mapping    ' num2str(i)])
mval = zeros(size(temp,1),ns);
p_d = (p_e).*tts + qcost;

for ri = 1:ns
    v_i = vi(:,ri:ns:k*ns);
    d_i = di(:,ri:ns:nj*ns);
    pd_i(:,ri) = p_d - subextra_id(:,ri).*max((x(1) + x(2)*d_i + x(3)*range_tt),0);
    mval(:,ri) = [log(pd_i(:,ri)) x1e]*theta1 + ehat;
    mu(:,ri) = ([log(pd_i(:,ri)),x2e(:,2:end)].*v_i*theta2w(:,1))...
        +[log(pd_i(:,ri)),x2e(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
end
subsidy_total = sum((tot_ded).*mktsizets);
rg_w = sum(range_tt(rw_ev, :).*sj_e(rw_ev, :))/sum(sj_e(rw_ev, :));

% Calculate CV for each city separately and weight by market size
cvs_city = zeros(length(unique(cdid_hy)),ns);
city_idx = 0;
sub_ind = zeros(length(unique(cdid_hy)), ns);
for t = min(cd_range):max(cd_range)
    trows = find(cdid_hy == t);
    time_trows = trows(hy_temp(trows) == hy_given);
    
    city_idx = city_idx + 1;
    city_mktsize = unique(mktsizets_hy(time_trows));
    
    subextra_id_city = subextra_id(time_trows,:);
    di_city = di(time_trows,:);
    range_tt_city = range_tt(time_trows,:);
    s_i_city = s_i(time_trows,:);
    sub_tt0_city = sub_tt0(time_trows);
    sub_ind(city_idx,:) = sum(subextra_id_city(sub_tt0_city>0,:).*max((x(1) + x(2)*di_city(sub_tt0_city>0,:) + x(3)*range_tt_city(sub_tt0_city>0,:)), 0)...
            .*s_i_city(sub_tt0_city>0, :))./sum(s_i_city(sub_tt0_city>0, :));
    % Calculate CV for this city
    cvs_city(city_idx,:) = -city_mktsize*mean((log(1 + sum(eg(time_trows,:))) - log(1 + sum(eg0(time_trows,:))))./...
        (alpha_i(time_trows,:)./pd_i(time_trows,:))) ;
end

sub_ind = reshape(sub_ind', 1, []);
% Total CV (sum of all cities' CV)
cv = sum(mean(cvs_city,2));
