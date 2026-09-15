% Replication driver for the alternative BLP-instrument demand specification.
%
% REPLICATION WORKFLOW
%   - Initialize paths, software dependencies, random seeds, and simulation
%     settings.
%   - Load and prepare the vehicle-market data.
%   - Generate random-coefficient and income simulation draws.
%   - Construct traditional BLP instruments and baseline IV estimates.
%   - Estimate demand by two-step GMM with income micro moments.
%   - Recover the covariance matrix and report the alternative demand results.
%
% EXECUTION NOTES
%   Run the complete script; it identifies the replication directory from its
%   own file location. The program uses global variables and intermediate MAT
%   files, so sections should be run in order unless all required checkpoints
%   have already been generated. The BLP-IV specification is computationally
%   intensive and starts a parallel pool. Random seeds and simulation settings
%   match the preferred specification so that the excluded instruments are the
%   intended source of variation.
%   This file should be opened in UTF-8 mode.
%
% PRINCIPAL OUTPUTS
%   BLP-IV demand estimates, standard errors, a formatted Table A3 diary
%   output (only BLP-IV demand estimates, the differentiation-IV specification
%   is already done in the main file), and a specification-specific MAT results file.

%% Initialize replication environment
clear all
clc
global invA ns x1 x2 s_jt IV IVs vfull dfull theta1 theti thetj cdid cdindex nobs nmkt sale_tax
global thet2 op xs  nd  NewW  niter nprod  tax fnmb nseudo seudos alpha_i nobs_c
global range_tt ev_tt  import_tt  vi di tts qcost x1e x2e ehat mktsizets
global  omega_sta cost_sim p_s_o local_tt sc_tt k NewincW inc_r s_inc
global brand_matrix city_matrix inc_c_W inc_b_W subdec subextra_id eg0  inc_IV
global GLOBAL_SCALE_FACTOR PROJECT_ROOT_DIR nobs_c cdid

PROJECT_ROOT_DIR = fileparts(mfilename('fullpath'));
cd(PROJECT_ROOT_DIR)
TOOLBOX_DIR = fullfile(PROJECT_ROOT_DIR, 'econometric toolbox');
% Add project functions and the bundled econometric toolbox.
addpath(genpath(TOOLBOX_DIR));
RUN_TAG = datestr(now, 'yyyymmdd_HHMMSS');
format short


GLOBAL_SCALE_FACTOR = 12/100000; % scaling factor for family annual income (hundred thousand)

if isempty(gcp('nocreate')) == 1
    parpool(24) % change according to the number of cores in your machine
end

rand('state', 1234);%4321);
randn('state', 1234);% 1234);  %fixed random draws to reproduce the estimation results

flag_gen = 1; % 0: load pre-generated random draws; 1: generate draws
flag_mail = 1; %1: set email notice
ns = 500;       % number of simulated "indviduals" per market %

%% Load and prepare vehicle-market data
DATA_FILE = fullfile(PROJECT_ROOT_DIR, 'data_incidence.csv');
data = dataset('file',DATA_FILE, 'Delimiter','comma','ReadVarNames',true);

% define variables
data = sortrows(data, [1,2,3,4,7,5]);
firm_name = data.firm_name_unique;
tbl = tabulate(firm_name);
t = cell2table(tbl,'VariableNames', ...
    {'Value','Count','Percent'});
t.Value = categorical(t.Value);
for i = 1:length(t.Value)
    fnmb(find(firm_name == t.Value(i)),:) = i;
end
firm_dum=dummify(fnmb);
firm_dum(:,1)=[];
p_id = data.p_id;
p_dum = dummify(p_id);
p_dum(:,1) = [];
city_id = data.city_d;
city_dum = dummify(city_id);
city_dum(:,1) = [];
hy = data.hy;
hy_dum = dummify(hy);
hy_dum(:,1) = [];
%id = data.id;
nobs = size(hy, 1);
sale = data.sales;
displace = data.displacement;
fc_fv = data.fuel_consumption;
fc_ev = data.kw_100km;
fuel_cost = data.cost_unit; % gas expenditure per 100 km
weight = data.weight;
power = data.power;
import = data.import;
brand = data.brand_d; % amer japa kore euro ///
brand_dum = dummify(brand);
IMb1 = brand(find(import == 1, 1, 'first'));
IMb = find(brand == IMb1);
[tempr tempc] = find(brand_dum(IMb,:) == 1);
brand_dum(:, tempc(1)) = [];
DOb1 = brand(find(import == 0, 1, 'first'));
DOb = find(brand == DOb1);
[tempr tempc] = find(brand_dum(DOb,:) == 1);
brand_dum(:, tempc(1)) = [];
% Drop brand fixed effects supported by only one subsample observation.
brand_dum = brand_dum(:, sum(brand_dum, 1) >= 2);
clear tempr tempc DOb1 DOb IMb1 IMb
lnweight = data.ln_weight;
lnpower = data.ln_power;
lnuni_cost = data.ln_cost_unit;
lnuni_cost_weight = lnuni_cost - lnweight;

lnpower_w = lnpower -lnweight;
lnwidth = data.ln_width;
lnsize = data.ln_size;
AT = data.AT;
ev = data.ev;
mpv = data.mpv;
suv = data.suv;

boc = data.boc;
boc_dum = dummify(boc);
boc_dum(:,1) = [];
euro = data.euro;
japan = data.japan;
usa = data.usa;
korea = data.korea;
price = data.msrp;
subsidy_local_net = data.subsidy_local_net;
quota_cost = data.quota_cost/10^4;
evrange = data.range;
op = price - subsidy_local_net;
op(find(ev == 0 ),1) = op(find(ev == 0 ),1)...
    *(1+1/1.13*.1) + quota_cost(find(ev == 0 ),1);
sale_tax = 1+ (ev == 0 )*(1/1.13*.1) ;
s_c = data.s_c;
p_s = price + s_c;
mktn = data.mktn;
iv_mkt_lnuni_cost = data.blp_mkt_cost_unit;
iv_firm_lnuni_cost = data.blp_firm_cost_unit;
iv_mkt_lnweight =data.blp_mkt_weight;
iv_firm_lnweight = data.blp_firm_weight;
iv_mkt_lnpower = data.blp_mkt_power;
iv_firm_lnpower = data.blp_firm_power;
iv_mkt_lnwidth = data.blp_mkt_width;
iv_firm_lnwidth = data.blp_firm_width;
iv_firm_lnsize = data.blp_firm_size;
iv_mkt_lnsize = data.blp_mkt_size;
etc = data.etc;
iv_local_msrp_hat = data.iv_local_msrp_hat;
iv_mkt_msrp_hat = data.iv_mkt_msrp_hat;
iv_local_power = data.iv_local_power;
iv_local_weight = data.iv_local_weight;
iv_local_size = data.iv_local_size;
iv_local_cost_unit = data.iv_local_cost_unit;
iv_firm_num = data.blp_firm_num;
iv_mkt_num = data.blp_mkt_num;
iv_cost = data.iv_cost;
% iv_firm_cw = data.blp_firm_c_w;
% iv_mkt_cw = data.blp_mkt_c_w;
odds = data.odds;
quota_d = data.quota_d;
%category = data.category;
local_brand = data.local_brand;
steel_hy = data.steel_hy;
mktsize = data.population*10^4/3; % city and town household number in 2000
s_jt = data.mkt_share*3;

year= data.year;
year_dum=dummify(year);
year_dum(:,1)=[];
ln_c_w = data.ln_c_w;
survey_inc = data.survey_income;
inc_r = find(survey_inc>100);
s_inc = survey_inc(inc_r,1);
tempbrd = brand*10 + ceil(hy/2);
tempbrdr = tempbrd(inc_r,:);
brand_matrix = dummify(tempbrdr);
tempcd = city_id*10 + ceil(hy/2);
tempcdr = tempcd(inc_r,:);
city_matrix = dummify(tempcdr);
clear tempbrdr tempbrd tempcd tempcdr
b2 = size(brand_matrix, 2);
tempsparse = [];
for i = 1: size(city_matrix,2)
    rk2 = rank([city_matrix(:,1:i), brand_matrix]);
    if rk2 < b2 + i
        tempsparse = [tempsparse, i];
        b2 = b2 + 1;
    end
end
city_matrix(:,tempsparse) = [];
clear rk2 b2 tempsparse



tax = zeros(nobs,1);
ttt=find(displace<=1 & ev == 0 );
tax(ttt)=(1.01*1.17);
ttt=find(displace>1 & displace<=1.5 & ev == 0 );
tax(ttt)=(1.03*1.17);
ttt=find(displace>1.5 & displace<=2.0 & ev == 0 );
tax(ttt)=(1.05*1.17);
ttt=find(displace>2.0 & displace<=2.5 & ev == 0 );
tax(ttt)=(1.09*1.17);
ttt=find(displace>2.5 & displace<=3.0 & ev == 0 );
tax(ttt)=(1.12*1.17);
ttt=find(displace>3.0 & displace<=4.0 & ev == 0);
tax(ttt)=(1.25*1.17);
ttt=find(displace>4 & ev == 0);
tax(ttt)=(1.40*1.17);
clear ttt

ttt = find(ev == 0);
tax_amount = zeros(size(ev,1),1);
tax_amount(ttt) = price(ttt)/1.17.*(tax(ttt) -1 );

x1 = [ones(nobs,1) lnpower lnuni_cost lnweight lnsize ev import AT...
    suv city_dum brand_dum hy_dum];
x2 = [log(op) ones(nobs,1) lnpower lnuni_cost lnweight lnsize];

%% Configure nonlinear demand parameters

theta2w=[1.2  1.85281
    -0.75   0
    0.16   0
    -0.8, 0
    -0.3, 0
    .1, 0];

% create a vector of the non-zero elements in the above matrix, and the %
% corresponding row and column indices. this facilitates passing values %
% to the functions below. %
[theti, thetj, theta2_0]=find(theta2w);

thet2 = length(theta2_0);

horz=['    mean       sigma      income '];
vert=['price     ';
    'const     ';
    'power     ';
    'fuel_cost ';
    'weight    ';
    'size      '];

%% Construct market indices and simulation draws
cdid = data.cdid;
hhinc_mu = data.hhinc_mu;
hhinc_sigma = data.hhinc_sigma_prov;
% this vector provides for each index of the last observation %
% in the data used here all brands appear in all markets. if this %
% is not the case the two vectors, cdid and cdindex, have to be   %
% created in a different fashion but the rest of the program works fine.%
nmkt = max(cdid);
nd = size(x1,2);
cdindex=zeros(nmkt,1);
nprod=zeros(nmkt,1);
for i=1:nmkt;
    cdindex(i) = find(cdid==i,1,'last');
    nprod(i)=length(find(cdid==i));
end

if flag_gen==1
    % do not comment out this part while generating randoms.
    temp = haltonset(nmkt,'Skip',1e3);
    halt = net(temp, ns*(size(x2,2)));
    draws=norminv(halt,0,1);

    clear temp halt
    halt_inc=haltonseq(ns +10,1)+ones(ns +10,1)*rand(1,1);
    temp=find(halt_inc>1);
    halt_inc(temp)=halt_inc(temp)-1;
    draws_inc=norminv(halt_inc,0,1);
    draws_inc(1:10,:)=[];
    
    seednum = draws_inc';
    demogr = zeros(nmkt, ns);
    for i = 1:nmkt
        v(i,:) =(draws(1:ns*size(x2,2),i))';
        tempdfull =  hhinc_sigma(cdindex(i))*seednum + hhinc_mu(cdindex(i))*ones(1, ns);
        demogr(i,:) = exp(winsor(tempdfull',[.1,99.9])');
        clear tempdfull
    end
    dfull = log(demogr(cdid,:)/100*2.7); % as family income
    vfull = v(cdid, :);

    filename=strcat('ps2_', num2str(ns));
    save (filename, 'x1', 'x2', 'demogr', 'v', 's_jt', 'op');
    clear halt halt_inc draws draws_inc temp

else
    filename=strcat('ps2_', num2str(ns));
    load(filename, 'demogr', 'v')
    dfull = log(demogr(cdid,:)/100*2.7); % as family income
    vfull = v(cdid, :);
end

%% Construct traditional BLP instruments and weighting inputs
% The remaining demand-estimation logic matches the preferred specification.
% Only this excluded-instrument set is changed.
iv = [iv_mkt_lnuni_cost iv_firm_lnuni_cost ...
    iv_mkt_lnweight iv_firm_lnweight ...
    iv_mkt_lnpower iv_firm_lnpower ...
    iv_mkt_lnsize iv_firm_lnsize ...
    mktn iv_cost iv_local_msrp_hat];

IV = [x1, iv];

inc_IV = IV(inc_r,:);
% Drop zero and singleton instrument columns in the micro-moment sample.
inc_IV = inc_IV(:, sum(inc_IV ~= 0, 1) >= 2);
nobs_c = size(inc_IV,1);
invA = inv([IV'*IV]);

%% Estimate baseline Logit and linear IV specifications
% compute the outside good market share by market
for i = 1:nmkt
    sum1(i,1) = sum(s_jt(find(cdid == i)));
end
outshr = 1.0 - sum1(cdid,:);

y = log(s_jt) - log(outshr);
mid = [x1, log(op)]'*IV*invA*IV';
bt = inv(mid*[x1, log(op)])*mid*y;
mvalold = [x1, log(op)]*bt;
oldt2 = zeros(size(theta2_0));
mvalold = exp(mvalold);

ols_results=ols(y, [log(op), x1]);
tsls_results=tsls(y, [log(op)], x1, IV);

theta1=tsls_results.beta;
[theta1(1:15), tsls_results.tstat(1:15)]
save theta1 theta1;
save mvalold_norm mvalold oldt2;
clear mid y outshr bt oldt2 mvalold temp sum1;


%% Estimate demand using the initial GMM weight
tic
niter=1;

fprintf('======GMM using initial weight=========== \n');

opts = optimoptions(@fminunc,'Algorithm','quasi-newton','StepTolerance',...
    1e-9,'TolFun',1e-7,'MaxFunctionEvaluations',1e9);

%% test 5 initial values
starts = 5;
theta20 = theta2_0;
%% original startvalues setup
startvalues = repmat(theta20',starts,1).* cat(2,ones(size(theta20)),abs(rand([starts-1,size(theta20',2)] )' ))';

GMPEC = 1.0e20;
CPUtMPEC = 0;
FuncEvalMPEC = 0;
LB = [zeros(size(theta20,1)-1, 1); -inf];
for reps=1:starts,
    %%% NOTE: START VALUES SET TO GIVE SAME INITIAL GMM OBJECTIVE FUNCTIONAL
    %%% VALUE AS WITH NFP
    theta20 = startvalues(reps,:)';  % starting values for standard deviations
    t1 = cputime;
    [theta2t,fval_rep,exitflag,output] = fminunc('gmmobj_micro_2', theta20, opts) ;

    CPUtMPEC_rep = cputime - t1;

    theta1MPEC_rep = theta1;
    theta2MPEC_rep = theta2t;
    GMPEC_rep = fval_rep;
    INFOMPEC_rep = exitflag;

    CPUtMPEC = CPUtMPEC + CPUtMPEC_rep;
    FuncEvalMPEC = FuncEvalMPEC + output.funcCount;
    if (GMPEC_rep < GMPEC && INFOMPEC_rep>0),
        thetaMPEC1 = theta1MPEC_rep;
        thetaMPEC2 = theta2MPEC_rep;
        GMPEC = GMPEC_rep;
        INFOMPEC = INFOMPEC_rep;
        load gmmresid
        GMMRESID = gmmresid;
        INCresid = inc_resid;
    end
end

theta2t = thetaMPEC2;

save theta2t0 theta2t

iv2=IV;
gmmresid = GMMRESID;

shat = iv2.*(gmmresid*ones(1,cols(iv2)));
NewW = shat'*shat; 
shatc = inc_IV.*(INCresid*ones(1,cols(inc_IV)));
inc_c_W = shatc'*shatc;

clear shat shatc

fprintf('======GMM using optimal weight=========== \n');

%% Estimate demand using the optimal GMM weight

niter=1; %reset the count of iterations

opts_OS = optimset('Display','iter', 'TolX',...
    1e-4,'TolFun',1e-4,'MaxFunEvals', 1e6, 'MaxIter', 1e6);

[theta2ds,fval, flag2]=fminsearch('gmmobj_microw_2', theta2t, opts_OS);

save theta2t1 theta2ds
vcov = var_cov_micro_2(theta2ds);
se = sqrt(diag(vcov));

theta2w = full(sparse(theti,thetj,theta2ds));
t1n = size(se,1) - size(theta2ds,1) ;
se2w = full(sparse(theti,thetj,se(t1n+1:size(se,1))));

[delta, expmu] = meanval_micro_PP_norm(theta2ds);

%% Report and save BLP-IV demand estimates
time = toc;

diary_path = fullfile(PROJECT_ROOT_DIR, 'diary');
if ~exist(diary_path, 'dir'), mkdir(diary_path); end
diary_file = fullfile(diary_path, ['demand_BLPIV_' RUN_TAG '.txt']);
diary off
diary(diary_file);
diary on

rep_blp_names = {'log(Price)', 'Constant', 'log(Horsepower)', ...
    'log(Operation costs)', 'log(Weight)', 'log(Size)', ...
    'EV', 'Import', 'AT', 'SUV'};
rep_blp_coef = nan(10,3);
rep_blp_se = nan(10,3);
rep_blp_coef(:,1) = ols_results.beta(1:10);
rep_blp_se(:,1) = abs(ols_results.beta(1:10)./ols_results.tstat(1:10));
rep_blp_coef(:,2) = tsls_results.beta(1:10);
rep_blp_se(:,2) = abs(tsls_results.beta(1:10)./tsls_results.tstat(1:10));
rep_blp_coef(:,3) = theta1(1:10);
rep_blp_se(:,3) = se(1:10);
rep_blp_tstat = rep_blp_coef./rep_blp_se;

rep_blp_random_names = {'log(Price)', 'Constant', 'log(Horsepower)', ...
    'log(Operation costs)', 'log(Weight)', 'log(Size)'};
rep_blp_random_coef = theta2w(1:6,1);
rep_blp_random_se = se2w(1:6,1);
rep_blp_random_tstat = rep_blp_random_coef./rep_blp_random_se;
rep_blp_income_coef = theta2w(1,2);
rep_blp_income_se = se2w(1,2);
rep_blp_income_tstat = rep_blp_income_coef/rep_blp_income_se;

fprintf('\n=============================================================================================\n');
fprintf('TABLE A3: Demand Estimation Results: BLP Random Coefficients Models\n');
fprintf('=============================================================================================\n');
fprintf('%-28s %34s %24s\n', ...
    'Variables', 'Homogeneous Preference', 'Heterogeneous Preference');
fprintf('%-28s %16s %16s %16s\n', '', 'OLS', '2SLS', 'GMM');
fprintf('%-28s %16s %16s %16s\n', '', '(1)', '(2)', '(3)');
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('%-28s %50s\n', '', 'Panel A: Mean Utility Parameters');
for rep_blp_i = 1:10
    rep_blp_coef_text = repmat({''},1,3);
    rep_blp_se_text = repmat({''},1,3);
    for rep_blp_j = 1:3
        if ~isnan(rep_blp_coef(rep_blp_i,rep_blp_j))
            rep_blp_stars = '';
            rep_blp_abs_t = abs(rep_blp_tstat(rep_blp_i,rep_blp_j));
            if rep_blp_abs_t >= 2.5758
                rep_blp_stars = '***';
            elseif rep_blp_abs_t >= 1.9600
                rep_blp_stars = '**';
            elseif rep_blp_abs_t >= 1.6449
                rep_blp_stars = '*';
            end
            rep_blp_coef_text{rep_blp_j} = sprintf('%.3f%s', ...
                rep_blp_coef(rep_blp_i,rep_blp_j), rep_blp_stars);
            rep_blp_se_text{rep_blp_j} = sprintf('(%.3f)', ...
                rep_blp_se(rep_blp_i,rep_blp_j));
        end
    end
    fprintf('%-28s %16s %16s %16s\n', ...
        rep_blp_names{rep_blp_i}, rep_blp_coef_text{:});
    fprintf('%-28s %16s %16s %16s\n', ...
        '', rep_blp_se_text{:});
end
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('%-28s %50s\n', '', 'Panel B: BLP Random Coefficients');
fprintf('%s\n', 'Variance parameters (sigmas)');
for rep_blp_i = 1:6
    rep_blp_stars = '';
    rep_blp_abs_t = abs(rep_blp_random_tstat(rep_blp_i));
    if rep_blp_abs_t >= 2.5758
        rep_blp_stars = '***';
    elseif rep_blp_abs_t >= 1.9600
        rep_blp_stars = '**';
    elseif rep_blp_abs_t >= 1.6449
        rep_blp_stars = '*';
    end
    fprintf('%-28s %16s %16s %16s\n', rep_blp_random_names{rep_blp_i}, ...
        '', '', sprintf('%.3f%s', rep_blp_random_coef(rep_blp_i), rep_blp_stars));
    fprintf('%-28s %16s %16s %16s\n', '', '', '', ...
        sprintf('(%.3f)', rep_blp_random_se(rep_blp_i)));
end
fprintf('%s\n', 'Demographic interactions');
rep_blp_stars = '';
rep_blp_abs_t = abs(rep_blp_income_tstat);
if rep_blp_abs_t >= 2.5758
    rep_blp_stars = '***';
elseif rep_blp_abs_t >= 1.9600
    rep_blp_stars = '**';
elseif rep_blp_abs_t >= 1.6449
    rep_blp_stars = '*';
end
fprintf('%-28s %16s %16s %16s\n', 'Income x log(Price)', '', '', ...
    sprintf('%.3f%s', rep_blp_income_coef, rep_blp_stars));
fprintf('%-28s %16s %16s %16s\n', '', '', '', ...
    sprintf('(%.3f)', rep_blp_income_se));
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('%-28s %16s %16s %16s\n', 'Instrumental Variables', '', ...
    'BLP IV', 'BLP IV');
fprintf('=============================================================================================\n');
fprintf('Observations: %d\n', nobs);
fprintf('Final GMM objective: %.6f\n', fval);
fprintf('Final optimizer exit flag: %d\n', flag2);
fprintf('Running time: %.2f minutes\n\n', time/60);

diary off

results_file = fullfile(PROJECT_ROOT_DIR, ['demand_BLPIV_' RUN_TAG '.mat']);
save(results_file, 'ols_results', 'tsls_results', 'theta1', 'theta2ds', ...
    'theta2w', 'vcov', 'se', 'se2w', 'fval', 'flag2', 'GMPEC', ...
    'INFOMPEC', 'IV', 'inc_IV', 'iv', '-v7.3');

clear rep_blp_names rep_blp_coef rep_blp_se rep_blp_tstat
clear rep_blp_random_names rep_blp_random_coef rep_blp_random_se
clear rep_blp_random_tstat rep_blp_income_coef rep_blp_income_se
clear rep_blp_income_tstat
clear rep_blp_i rep_blp_j rep_blp_coef_text rep_blp_se_text
clear rep_blp_stars rep_blp_abs_t
