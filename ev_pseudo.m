function f = ev_pseudo(nobs_c, nseudo)
% EV_PSEUDO Generate Extreme-Value taste-shock draws for policy simulations.
%
% PURPOSE
%   Construct the pseudo draws used in the counterfactual demand simulations.
%   The first column is a transformed uniform Extreme-Value draw; later
%   columns select between the current and lagged candidates using the density
%   ratio implemented below.
%
% INPUTS
%   nobs_c - positive integer number of simulated observations/consumers.
%   nseudo - positive integer number of pseudo draws per observation.
%
% OUTPUT
%   f / epsilon - [nobs_c x nseudo] simulated Extreme-Value disturbances.
%
uni_seed = rand(nobs_c, nseudo);
n_seed = rand(nobs_c, nseudo); 
epsilon_tilda = -log(-log(uni_seed)); 
f = zeros(nobs_c, nseudo); 
f(:,1) = epsilon_tilda(:,1); 
for i = 2: nseudo
    pdf_ev = evpdf(epsilon_tilda(:,i), 0, 1); 
    pdf_gev = gevpdf(epsilon_tilda(:,i), 0, 1, 0);
    pdf_ev_lag = evpdf(epsilon_tilda(:,i-1), 0, 1); 
    pdf_gev_lag = gevpdf(epsilon_tilda(:,i-1), 0, 1, 0);
    I_temp = (n_seed(:,i) < (pdf_gev./pdf_ev)./(pdf_gev_lag./pdf_ev_lag));
    f(:,i) = epsilon_tilda(:,i).*I_temp + epsilon_tilda(:,i-1).*(1-I_temp);
end
    
