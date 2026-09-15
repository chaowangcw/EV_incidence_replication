function f = optsb2(x)
% OPTSB2 Evaluate the objective for the parameterized subsidy-rule search.
%
% PURPOSE
%   Objective called by FMINCON in the policy counterfactual. For proposed
%   subsidy coefficients, solve the Bertrand price fixed point, compute shares
%   and subsidy payments, and return the scalar criterion implemented below.
%
% SYNTAX
%   objective = optsb2(x)
%
% INPUT
%   x - [3 x 1] subsidy-rule coefficients. Individual subsidies equal
%       max(x(1)+x(2)*income+x(3)*range,0) for eligible observations. Units
%       inherit the scaling of DI, RANGE_TT, and prices in the main script.
%
% OUTPUT
%   f - scalar criterion implemented as
%       abs(subsidy_total-subsidy_0)-sum(sj_e.*mktsizets). Invalid or complex
%       prices return 1e10.
%
% REQUIRED GLOBAL STATE
%   Uses demand estimates, simulation draws, prices and costs, market sizes,
%   eligibility indicators, ownership matrices, and market IDs initialized by
%   the main script; see the GLOBAL declarations below.
%
% FILE INPUT
%   Loads opt_subsidy.mat. It must contain tol, temp, theta2w, theta1, nj,
%   subsidy_0, sub_tt0, p_s_0, hy, cd_range, and model_variant.

global range_tt  vi di tts qcost x1e x2e ehat mktsizets ev_tt import_tt
global  omega_sta cost_sim p_s_o local_tt sc_tt  k  ns alpha_i subdec subextra_id cdid

load opt_subsidy

norm = 1;
avgnorm = 1;
i = 0;
p_e = p_s_0;

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
            pd_i(:,ri) = p_d - subextra_id(:,ri).*max((x(1) + x(2)*d_i + x(3)*range_tt),0);
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
        tot_ded = mean(subextra_id.*max((x(1) + x(2)*di + x(3)*range_tt),0).*s_i,2);
        
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
if max(isnan(p_d)) == 1 | isreal(log(p_d)) == 0
    f =       1e10;
else
    subsidy_total = sum((tot_ded).*mktsizets);
    rg_w = sum(range_tt.*sj_e)/sum(sj_e);
    f = abs(subsidy_total - subsidy_0) - sum(sj_e.*mktsizets);
end
