% Replication driver for the random-coefficients vehicle-demand model,
% supply-side markup recovery, and EV-subsidy counterfactual analysis.
%
% REPLICATION WORKFLOW
%   - Initialize paths, software dependencies, random seeds, and simulation
%     settings.
%   - Load and prepare the vehicle-market data.
%   - Generate random-coefficient and income simulation draws.
%   - Construct differentiation instruments and baseline IV estimates.
%   - Estimate demand by two-step GMM with income micro moments.
%   - Recover standard errors, elasticities, marginal costs, and markups.
%   - Run subsidy counterfactuals and optimize an alternative subsidy rule.
%   - Produce policy summaries, distributional results, and figures.
%
% EXECUTION NOTES
%   Run the complete script; it identifies the replication directory from its
%   own file location. The program uses global variables and intermediate MAT
%   files, so sections should be run in order unless all required checkpoints
%   have already been generated. The default specification is computationally
%   intensive and starts a parallel pool. Random seeds are set below for
%   reproducibility.
%   This file should be opened in UTF-8 mode.
%
% PRINCIPAL OUTPUTS
%   Demand estimates and covariance matrices, model-level markups, policy
%   simulation results, diary logs, MAT checkpoints, and EPS figures.

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
    parpool(24)
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

%% Construct differentiation instruments and weighting inputs
x_diffiv = [lnweight,lnpower,lnuni_cost,lnsize];
% for local and quadratic diffIVs, the first half is other products within the firm, and the second
% half is rival products
[diffIV_local_interact, diffIV_local_interact_colnames]= build_differentiation_iv(x_diffiv,cdid,fnmb,'local',true);
[diffIV_local, diffIV_local_colnames]= build_differentiation_iv(x_diffiv,cdid,fnmb,'local',false);
[diffIV_quadratic, diffIV_quadratic_colnames]= build_differentiation_iv(x_diffiv,cdid,fnmb,'quadratic',true);
diffIV_local = [diffIV_local, diffIV_local_interact];

% using differentiation IV from Gandhi and Houde
iv = [
    mktn iv_cost ...
     iv_local_msrp_hat ...
     diffIV_quadratic(:,[2,12]) diffIV_local(:,[6,8,38]) ...
];

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

%% Recover model-level markups and report demand estimates
% Call unified markup calculation function (compatible with demand estimation model)
mkup = markup_ds_uniform(theta2ds, delta, expmu, hy, data.model_variant, mktsize);
% Calculate markup rate
mkup_rate = mkup./(price + s_c);

diary_path = fullfile(PROJECT_ROOT_DIR, 'diary');
if ~exist(diary_path, 'dir'), mkdir(diary_path); end 
disp(['diary path is set as: ' diary_path]);
diary(fullfile(diary_path, ['main_UniformMSRP_' RUN_TAG '.txt']));

diary on
datastat = [[mean([sale, op, weight, power, fuel_cost])]',...
    [std([sale, op, weight,  power, fuel_cost])]',...
    [min([sale, op, weight,  power, fuel_cost])]',...
    [max([sale, op, weight,  power, fuel_cost])]'];

display('Results of Berry logit');
num2str([ols_results.beta(1:13) ols_results.tstat(1:13) tsls_results.beta(1:13) tsls_results.tstat(1:13)])

disp(horz)
disp('  ')
for i=1:size(theta2w,1)
    disp(vert(i,:))
    disp([theta1(i) theta2w(i,:)])
    disp([se(i) se2w(i,:)])
end

disp(['markup rate: average   max   min '])
[mean(mkup_rate)   max(mkup_rate) min(mkup_rate)]

vert_mean = ['price  '
    'const  '
    'power  '
    'fuel_co'
    'weight '
    'size   '
    'ev     '
    'import '
    'AT     '
    'suv    '];

horz_mean = ['variable',',', 'OLS coeff' ,',',  'OLS std' ,',', 'TSLS coeff', ',', 'TSLS std', ',', 'Full model coeff',',',  'Full model std'];
R0= size(vert_mean, 1);
disp(horz_mean)
disp('  ')
for i=1:R0
    disp([vert_mean(i,:),',', num2str(ols_results.beta(i)) ,',', '(' num2str(ols_results.beta(i)/ols_results.tstat(i)) ')',...
        ',', num2str(tsls_results.beta(i)),',', '(' num2str(tsls_results.beta(i)/tsls_results.tstat(i)) ')',...
        ',', num2str(theta1(i)),',', '(' num2str(se(i)) ')'])
end

time = toc;

fprintf('running time: ', time/60, ' min\n')

diary off

save(fullfile(PROJECT_ROOT_DIR, ['main_UniformMSRP_' RUN_TAG '.mat']), '-v7.3')


%% Compute elasticities and descriptive supply-side results
load incresid_sep % residuals from demand estimation

% experiment using Shanghai, second half year of 2019
etc = (data.p_id == 2 & hy == 8); % change the experiment city to beijing
mu = [];
[n k] = size(x2);
nj = size(theta2w,2)-1;
sigmap = theta2w(1, 1);
betapi = theta2w(1, 2:nj+1)';
epsilon = zeros(nobs,1);
vfull = v(cdid, :);
for et = 1:max(cdid)
    temp = find(cdid == et);
    x1e = x1(temp,:); p_e = log(op(temp));
    x2e = x2(temp,:);
    vi = vfull(temp,:); 
    di = exp(dfull(temp, :));
    mval = x1e*theta1(2:end) + theta1(1)*p_e;
    for ri = 1:ns
        v_i = vi(:,ri:ns:k*ns);
        d_i = di(:,ri:ns:nj*ns);
        mu(:,ri) = (x2e.*v_i*theta2w(:,1))+x2e.*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
    end
    eg = exp(mu + kron(ones(1,ns),mval));
    s_i = eg./(ones(size(temp,1),1)*(1 + sum(eg)));
    demo_i = zeros(size(temp,1), ns);
    for n = 1:nj
        demo_i = di(:, (n-1)*ns + 1:n*ns)*betapi(n,:) + demo_i;
    end
    temp0 = (s_i.*(theta1(1) + vi(:,1:ns)*sigmap + demo_i))*s_i'/ns;
    temp1 = sum((s_i.*(theta1(1) + vi(:,1:ns)*sigmap + demo_i))')'/ns;

    H = diag(temp1) - temp0;
    epsilon(temp,:) = diag(H)./mean(s_i,2);
    if min(etc(temp,:)) == 1 % highest frequency of ev == 1 & import == 1
        cep = H./(mean(s_i,2)*ones(1,length(temp)));
    end
    clear temp temp0 temp1 mu
end
[esp_mean, esp_std, esp_num, esp_min, esp_max] = grpstats(epsilon, cdid, {'mean','std','numel','min', 'max'});
wepsilon = epsilon.*sale;
[weps_mean] = grpstats(wepsilon, cdid, {'mean'});
[sale_m] = grpstats(sale, cdid, {'mean'});
welasticity = weps_mean./sale_m;

% compute markup rate for ICEV and EV
wmark_rate = mkup_rate.*sale;
[wmark_rate_mean] = grpstats(wmark_rate, ev, {'mean'});
[sale_evICEV] = grpstats(sale, ev, {'mean'});
wmkup = wmark_rate_mean./sale_evICEV;

hy_given = 8;
temp = find(data.hy == hy_given); % p_id is province id, hy is half year

% marginal cost is assumed to be unchanged.
mc = (price + s_c) - mkup; % marginal cost estimates

% To generate all subsidy following table 1===================
% Initialize full-sample subsidy
temp_sc = zeros(size(price));
% Eligibility
base_ev = import == 0 & ev == 1;
% t = 1 or 2
idx = base_ev & ismember(hy, [1 2]);
temp_sc(idx & evrange < 100) = 0;
temp_sc(idx & evrange >= 100 & evrange < 150) = 2.5;
temp_sc(idx & evrange >= 150 & evrange < 250) = 4.5;
temp_sc(idx & evrange >= 250) = 5.5;
% t = 3, 4, or 5
idx = base_ev & ismember(hy, [3 4 5]);
temp_sc(idx & evrange < 100) = 0;
temp_sc(idx & evrange >= 100 & evrange < 150) = 2;
temp_sc(idx & evrange >= 150 & evrange < 250) = 3.6;
temp_sc(idx & evrange >= 250) = 4.4;
% t = 6 or 7
idx = base_ev & ismember(hy, [6 7]);
temp_sc(idx & evrange < 150) = 0;
temp_sc(idx & evrange >= 150 & evrange < 200) = 1.5;
temp_sc(idx & evrange >= 200 & evrange < 250) = 2.4;
temp_sc(idx & evrange >= 250 & evrange < 300) = 3.4;
temp_sc(idx & evrange >= 300 & evrange < 400) = 4.5;
temp_sc(idx & evrange >= 400) = 5;
% t = 8
idx = base_ev & hy == 8;
temp_sc(idx & evrange < 250) = 0;
temp_sc(idx & evrange >= 250 & evrange < 400) = 1.8;
temp_sc(idx & evrange >= 400) = 2.5;
% Optional check
if any(~ismember(hy, 1:8))
    error("unknown halfyear")
end
% Cap subsidy
sc = min([temp_sc, price * .6], [], 2);
sub_all0 = sc + subsidy_local_net;
clear temp_sc idx
sub_all0_rate = sub_all0./(price + s_c);
% compute average markup for ICEV and EV, and average subsidy for EV
mark_welfare = op - mc; % to compute welfare, we need the markup defined by price out of consumer pocket
wmark_welfare = mark_welfare.*sale;
wsub0 = sub_all0.*sale;
tbl = table(ev, wmark_welfare, wsub0, sale);
tblstats1 = grpstats(tbl,"ev","mean","DataVars",["wmark_welfare","wsub0","sale"])
wmark_welfare_mean = grpstats(wmark_welfare, ev,"mean");
sale_evICEV = grpstats(sale,ev,"mean");
wsub0_mean = grpstats(wsub0,ev,"mean");
wmarkup_welfare = wmark_welfare_mean./sale_evICEV;
wsubsidy0_welfare = wsub0_mean./sale_evICEV;
% just for computing hy == 8
wmark_welfare_mean8 = grpstats(wmark_welfare(temp), ev(temp), "mean");
sale_evICEV8 = grpstats(sale(temp), ev(temp), "mean");
wsub0_mean8 = grpstats(wsub0(temp), ev(temp), "mean");
wmarkup_welfare8 = wmark_welfare_mean8./sale_evICEV8;
wsubsidy0_welfare8 = wsub0_mean8./sale_evICEV8;
wprice = (price+s_c).*sale;
wprice_mean = grpstats(wprice(temp), ev(temp), "mean");
wprice_weighted = wprice_mean./sale_evICEV8;
wmc = mc.*sale;
wmc_weighted = grpstats(wmc(temp), ev(temp), "mean")./sale_evICEV8
wprice_weighted - wmc_weighted
wqcost_weighted = grpstats(quota_cost(temp).*sale(temp), ev(temp),"mean")./sale_evICEV8
wmsrp_weighted = grpstats(price(temp).*sale(temp), ev(temp),"mean")./sale_evICEV8
wtax = (price-subsidy_local_net).*1/1.13*0.1.*sale;
wtax_weighted = grpstats(wtax(temp), ev(temp),"mean")./sale_evICEV8
wop_weighted = grpstats(op(temp).*sale(temp), ev(temp),"mean")./sale_evICEV8
wop_weighted - wmc_weighted
wmkup_weighted = grpstats(mkup(temp).*sale(temp), ev(temp), "mean")./sale_evICEV8
original_sub = s_c + subsidy_local_net;
woriginal_sub_weighted = grpstats(original_sub(temp).*sale(temp), ev(temp), "mean")./sale_evICEV8
% all subsidy generate end================================================

temp = find(etc == 1);
% Cross-price elasticity at the sub-brand level for selected brands. Numeric
% brand codes preserve the original selection after vehicle names are removed.
brand_c = data.brand_d(temp);
figure_brand_codes = [6; 76; 5; 11; 19; 60; 3; 49; 57; 118; 54; 87];
pp = price(temp);
sale_pp = sale(temp);
brand_id = zeros(length(temp),1);
nsbr = numel(figure_brand_codes);
for i = 1:nsbr
    tb = find(brand_c == figure_brand_codes(i));
    if isempty(tb)
        error('Figure 5 brand code %d is absent from the selected market.', ...
            figure_brand_codes(i));
    end

    p_m(i,:) = mean(pp(tb)); b_m(i,:) = figure_brand_codes(i);
    sale_mpp(i,:) = sum(sale_pp(tb));
    brand_id(tb) = i;
end
[tempr, tempc, br_id] = find(brand_id);
reduced_cep = cep(tempr, tempr);
for i = 1:nsbr
    tbr_i = find(br_id == i);
    for j = 1:nsbr
        tbr_j = find(br_id == j);
        cep_brand(i,j) = mean(sum(reduced_cep(tbr_i, tbr_j), 1));
    end
end

Amt = [[zeros(7,1); ones(nsbr-7,1)], p_m];
[Bmt, index] = sortrows(Amt, [1,2]);
diary on
disp(['brand level cross elasticity'])
[p_m(index,:), cep_brand(index, index)]
b_m(index, :)
diary off

%% Paper Figure 5: Price elasticity of selected EVs
% Replication output:
%   elasticity.eps (Figure 5 in the paper).
% Figure content:
%   Scatter plot of brand-level quantity against price for domestic and imported
%   vehicles. Point labels report the corresponding own-price elasticity.
b_m_new = ["FAW-Audi","BYD","FAW-VW","SAIC-VW","Dongfeng-Honda","GAC-Toyota","FAW-Toyota","VW (imported)","BMW (imported)","Lexus (imported)","Audi (imported)","Tesla (imported)"]
elastable = table(b_m_new',[zeros(7,1); ones(nsbr-7,1)], p_m*10, sale_mpp, diag(cep_brand))
elastable(6:7,:) = [];
figure;
domestic = (elastable.Var2 == 0);
scatter(elastable.sale_mpp(domestic),elastable.Var3(domestic), 'go', 'filled');
hold on;
imported = (elastable.Var2 == 1);
scatter(elastable.sale_mpp(imported), elastable.Var3(imported), 'b^', 'filled');
for i = 1:10
    label = sprintf('%s\n%.2f', elastable.Var1{i}, elastable.Var5(i));
    text(elastable.sale_mpp(i), elastable.Var3(i), label, ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 8);
end

xlabel('Quantity');
ylabel('Price (RMB ,000)');
legend('domestic', 'imported');
hold off
saveas(gcf, 'elasticity','epsc');

ev_tt = ev(temp); import_tt = import(temp);
fv_tt = 1 - ev_tt ;


for ri = 1:1
    if ri == 1
        varname = 'ev_tt';
    end
    tempr = find(import_tt == 0 & eval(varname) == 1); % first dimension
    for im_i = 0:1
        for ci = 1:2
            if ci == 1
                cvarname = 'ev_tt';
            else
                cvarname = 'fv_tt';
            end
            tempc = find(import_tt == im_i & eval(cvarname) == 1);
            if isempty(tempc) == 0
                tempels = cep(tempr, tempc);
                if  ri == ci & im_i == 0
                    tempels_own = diag(diag(tempels));
                    elsmat = tempels - tempels_own;
                else
                    elsmat = tempels;
                end
                [tempi, tempj, epsilon_non_zero] = find(elsmat);
                epsilon_avg(ri,im_i*2 + ci) = mean(epsilon_non_zero);
                epsilon_sum(ri,im_i*2 + ci) = mean(sum(elsmat,2));
            else
                epsilon_avg(ri,im_i*2 + ci) = 0;
                epsilon_sum(ri,im_i*2 + ci) = 0;
            end
        end
    end
end
diary on
vert_name_els = ['average elasticity  ';
    'aggregate elasticity'];
disp(['elasticities', 'ev_d', 'fv_d', 'ev_i', 'fv_i'])
disp([vert_name_els(1,:), num2str(epsilon_avg(1,:))])
disp([vert_name_els(2,:), num2str(epsilon_sum(1,:))])
diary off


%% Prepare subsidy counterfactual inputs
% Purpose:
%   Select the policy period and sample, construct observed and alternative
%   subsidy schedules, and prepare price, cost, income, ownership, and demand
%   objects for the counterfactual equilibrium calculations.

hy_given = 8;
cd_range = unique(cdid(find(hy==hy_given))); % choosing 2019 the second half year
temp = find(data.hy == hy_given); % p_id is province id, hy is half year

nobs_c = length(temp);
p_e = price(temp);
range_tt = evrange(temp);
ev_tt = ev(temp); import_tt = import(temp);
 fv_tt = 1 - ev_tt ;
temp_sc = zeros(nobs_c,1);
temp_sc(import_tt == 0 & ev_tt == 1 & range_tt <250) = 0;
temp_sc(import_tt == 0 & ev_tt == 1 & range_tt >=250 & range_tt < 400) = 1.8;
temp_sc(import_tt == 0 & ev_tt == 1 & range_tt >=400) = 2.5;
sc_tt = min([temp_sc, p_e*.6],[],2);

local_tt = subsidy_local_net(temp);
sub_tt0 = sc_tt + local_tt;

subsidy_h = zeros(size(temp));
subsidy_l = zeros(size(temp));
temp_h = zeros(size(temp));
temp_h(range_tt>=80 & range_tt < 150 & ev_tt==1 & import_tt ==0) = 3.15;
temp_h(range_tt>=150 & range_tt < 250 & ev_tt==1 & import_tt ==0) = 4.5;
temp_h(range_tt>=250 & ev_tt==1 & import_tt ==0) = 5.4;
subsidy_h = min([temp_h, p_e*.6/1.6],[],2);
clear temp_h
temp_h = zeros(size(temp));
temp_h(range_tt>=80 & range_tt < 150 & ev_tt==1) = 3.15;
temp_h(range_tt>=150 & range_tt < 250 & ev_tt==1) = 4.5;
temp_h(range_tt>=250 & ev_tt==1) = 5.4;
subsidy_a = min([temp_h, p_e*.6/1.6],[],2);
clear temp_h

%% Run benchmark subsidy counterfactuals
%   Solve Bertrand price equilibria under the benchmark subsidy scenarios and
%   calculate prices, shares, sales, profits, fiscal costs, externalities, and
%   compensating variation.

mc_ex_ev = .165*7; % 7 is the exchange rate; $.165/kw
mc_ex_ev2 = .0352*7;
mc_ex_fv = .11*7; % $.11/liter is the external costs for emission only.
lifetime = 15;
vmt = 17988;

cost_sim = mc(temp);
[n k] = size(x2);
nj = size(theta2w,2)-1;
sigmap = theta2w(1, 1);

nseudo = 1000;
seudos = ev_pseudo(nobs_c + 1, nseudo);
tol=1e-13;
maxit = 500;
shr_e = []; price_e = [];
pi_e = []; pi_e0 = [];
fi_diff = [];
pi_category = [];
pr_sim = [];
cvs = [];
avgprice_d =  zeros(5,4); avgprice_s = zeros(5,4);
pd_mat = zeros(nobs_c,5); ps_mat = zeros(nobs_c, 5);
sales_sim = zeros(5,4); mktshr_e =zeros(nobs_c,4);
subsidy_c = zeros(5,1);
ex_ev = zeros(5,1);  ex_fv = zeros(5,1);
ex_evd = zeros(5,1); 
ehat = gmmresid(temp,:);
mktsizets = mktsize(temp); % in person
qcost = quota_cost(temp).*(ev(temp) == 0);
tts = (1 + 1/1.13*.1*(ev(temp) == 0));
x1e = x1(temp,:);
x2e = x2(temp,:);
f_i = fnmb(temp,:);
cdid_hy = cdid(temp);
hy_temp = hy(temp);
omega_sta = cell(length(unique(cdid_hy)),1);
for t = min(cd_range):max(cd_range)
    trows = find(cdid_hy == t); tempn = size(trows,1);
    time_trows = trows(hy_temp(trows) == hy_given);
    omega_sta_c  = zeros(tempn);
    f_i_city = f_i(time_trows,:);
    for j = 1:tempn
        for h = 1:tempn
            omega_sta_c(j,h) = (f_i_city(j) == f_i_city(h));
        end
    end
    omega_sta{t} = omega_sta_c;
end
vi = vfull(temp,:);
di = exp(dfull(temp, :)); % here is still family income

demo_i = zeros(size(temp,1), ns);
for n = 1:nj
    demo_i = di(:, (n-1)*ns + 1:n*ns)*betapi(n,:) + demo_i;
end
alpha_i = theta1(1) + vi(:,1:ns)*sigmap + demo_i;

%% income tax
tempsub =  (di*100/2.7 > 60 & di*100/2.7 <= 90)*.03 +...
            (di*100/2.7 > 90 & di*100/2.7 <= 144)*.03 +...
                    (di*100/2.7 > 144 & di*100/2.7 <= 300)*.05 + (di*100/2.7 > 300)*.05;
subdec = tempsub.*(sub_tt0*ones(1,ns));
subextra_id = (sub_tt0>0)*ones(1,ns);

p_s_o = price(temp) + sc_tt; % msrp + subsidy to manufacturer = actual price paid to manufacturer; it is not zero here.
norm = 1;
avgnorm = 1;
i = 0;
mu = [];
p_e = p_s_o;
mid = data.model_variant(temp);
population_hy = data.population(temp);
mktsizets_hy = mktsizets;
while norm > 1e-1*tol*10^(floor(i/50)) & avgnorm > 1e-3*tol*10^(floor(i/50))
    p_d = (p_e- sub_tt0).*tts + qcost;
    mval = [log(p_d), x1e]*theta1 + ehat;
    for ri = 1:ns
        v_i = vi(:,ri:ns:k*ns);
        d_i = di(:,ri:ns:nj*ns);
        mu(:,ri) = ([log(p_d),x2e(:,2:end)].*v_i*theta2w(:,1))+[log(p_d),x2e(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
    end
    MU0 = mu + kron(ones(1,ns),mval);
    eg0 = exp(MU0);
    
    % Calculate market share for each city separately
    s_i = zeros(size(eg0));
    
    hy_omega = zeros(max(unique(mid)));  
    hy_s = zeros(max(unique(mid)),1);
    
    for t = min(cd_range):max(cd_range)
        trows = find(cdid_hy == t);
        time_trows = trows(hy_temp(trows) == hy_given);
        
        % Calculate market share for this city
        s_i(time_trows,:) = eg0(time_trows,:)./(ones(length(time_trows),1)*(1 + sum(eg0(time_trows,:))));
        
        % Calculate omega matrix and aggregate
        tempn = size(time_trows,1);
        temp0 = (s_i(time_trows,:).*alpha_i(time_trows,:))*s_i(time_trows,:)'/ns;
        temp1 = sum((s_i(time_trows,:).*alpha_i(time_trows,:))')'/ns;
        H = (diag(temp1) - temp0)./(p_d(time_trows,:)*ones(1,tempn)).*tts(time_trows,:);
        
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
    sj_e0 = mean(s_i,2);
    
    nanproduct = find(hy_s==0);
    hy_s(nanproduct)=[];
    hy_omega(nanproduct,:) =[];
    hy_omega(:,nanproduct)=[];
    
    mkup_model = -inv(hy_omega)*hy_s;
    
    temp_m = [1:max(unique(mid))]';
    temp_m(nanproduct)=[];
    
    mkup_e0 = zeros(length(p_e),1);
    for index = 1:length(temp_m)
        mkup_e0(find(mid==temp_m(index)))= mkup_model(index);
    end
    
    pnew = mkup_e0 + cost_sim;

    norm = max(abs(pnew - p_e));
    avgnorm = mean(abs(pnew - p_e));
    p_e = pnew;

    i = i + 1;
end
p_d_0 = (p_e - sub_tt0).*tts + qcost;
pd_mat(:,1) = p_d_0;

mval = [log(p_d_0), x1e]*theta1 + ehat;
for ri = 1:ns
    v_i = vi(:,ri:ns:k*ns);
    d_i = di(:,ri:ns:nj*ns);
    mu(:,ri) = ([log(p_d_0),x2e(:,2:end)].*v_i*theta2w(:,1))+[log(p_d_0),x2e(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
end
MU0 = mu + kron(ones(1,ns),mval);
clear mu mval
p_s_0 = p_e; 
ps_mat(:,1) = p_s_0;
avgprice_d(1,:) = [mean(p_d_0(find(ev_tt == 1 & import_tt == 0))),...
    mean(p_d_0(find(fv_tt == 1 & import_tt == 0))),...
    mean(p_d_0(find(ev_tt == 1 & import_tt == 1))),...
    mean(p_d_0(find(fv_tt == 1 & import_tt == 1)))];
avgprice_s(1,:) = [mean(p_s_0(find(ev_tt == 1 & import_tt == 0))),...
    mean(p_s_0(find(fv_tt == 1& import_tt == 0))),...
    mean(p_s_0(find(ev_tt == 1 & import_tt == 1))),...
    mean(p_s_0(find(fv_tt == 1& import_tt == 1)))];
s_i_0 = s_i; % save for later drawing income vs ev adoption rate

% Calculate sales and sub_inc for each city separately
sub_inc_city = zeros(length(unique(cdid_hy)), ns);
expected_sub_city = zeros(length(unique(cdid_hy)), ns); % expected subsidy (unconditional subsidy)
sales_ev0 = 0; sales_fv0 = 0; sales_ev1 = 0; sales_fv1 = 0;
city_idx = 0;
sub_inc = zeros(5,length(unique(cdid_hy))*ns);
for t = min(cd_range):max(cd_range)
    trows = find(cdid_hy == t);
    time_trows = trows(hy_temp(trows) == hy_given);
    city_mktsize = unique(mktsizets_hy(time_trows));
    
    city_idx = city_idx + 1;
    
    % Calculate sub_inc for this city
    s_i_city = s_i(time_trows,:);
    sub_tt0_city = sub_tt0(time_trows);
    ev_tt_city = ev_tt(time_trows);
    import_tt_city = import_tt(time_trows);
    ev_dom_idx = (ev_tt_city == 1 & import_tt_city == 0);
    sub_inc_city(city_idx,:) = sum(s_i_city(ev_dom_idx, :).*...
        (sub_tt0_city(ev_dom_idx)*ones(1, ns)))./sum(s_i_city(ev_dom_idx, :));
    expected_sub_city(city_idx,:) = sum(s_i_city(ev_dom_idx, :) .* ...
        (sub_tt0_city(ev_dom_idx) * ones(1, ns)));

    % Calculate sales for this city
    sales_ev0 = sales_ev0 + sum(sj_e0(time_trows(find(ev_tt(time_trows) == 1 & import_tt(time_trows) == 0))))*city_mktsize;
    sales_fv0 = sales_fv0 + sum(sj_e0(time_trows(find(fv_tt(time_trows) == 1 & import_tt(time_trows) == 0))))*city_mktsize;
    sales_ev1 = sales_ev1 + sum(sj_e0(time_trows(find(ev_tt(time_trows) == 1 & import_tt(time_trows) == 1))))*city_mktsize;
    sales_fv1 = sales_fv1 + sum(sj_e0(time_trows(find(fv_tt(time_trows) == 1 & import_tt(time_trows) == 1))))*city_mktsize;
end
sub_inc(1,:) = reshape(sub_inc_city(:, :)', 1, []);
sales_sim(1,:) = [sales_ev0, sales_fv0, sales_ev1, sales_fv1];
pi_e0(1,1) = sum(sj_e0.*mkup_e0.*mktsizets);
subsidy_c(1)= sum((local_tt + sc_tt).*sj_e0.*mktsizets);
ex_fv(1) = sum(sj_e0.*mktsizets.*fc_fv(temp).*(1-ev_tt)*mc_ex_fv*lifetime*vmt);
ex_ev(1) = sum(sj_e0.*mktsizets.*fc_ev(temp).*ev_tt*mc_ex_ev*lifetime*vmt);
ex_evd(1) = sum(sj_e0.*mktsizets.*fc_ev(temp).*ev_tt*mc_ex_ev2*lifetime*vmt);

% check whether the sales correspond to the ev sales
compare = [ev_tt, import_tt, s_jt(temp), sj_e0, mktsizets, sale(temp), s_jt(temp).*mktsizets, sj_e0.*mktsizets];
% compute the average subsidy over income
% sub_iMat = [reshape(demogr(cd_range,:)',1,[]);sub_inc(1,:)]';
sub_iMat = [reshape(demogr(cd_range,:)',1,[]);10000*reshape(expected_sub_city',1,[])]';
expected_sub_sort = sort(sub_iMat,1);
expected_sub_prct = prctile(expected_sub_sort(:,1), [25,50,75] );
expected_sub_quantile = [ mean( sub_iMat(sub_iMat(:,1) <= expected_sub_prct(1) ,2)), ...
  mean( sub_iMat(expected_sub_prct(2) >= sub_iMat(:,1) &  expected_sub_prct(1) < sub_iMat(:,1) ,2)), ...
  mean( sub_iMat(expected_sub_prct(3) >= sub_iMat(:,1) &  expected_sub_prct(2) < sub_iMat(:,1) ,2)), ... 
  mean( sub_iMat(sub_iMat(:,1) > expected_sub_prct(3) ,2)) ];


% Identify different firms for profit aggregation by firm
% f_i: firm IDs for all products in temp (hy=8 across all cities)
% nobs_c: number of observations in temp
pi_category = [];
[fi_sorted, I_fi] = sortrows(f_i);
fi_diff(2:nobs_c,:) = diff(fi_sorted);
fi_diff(1,:) = fi_sorted(1);
[r_c, c_c, ftemp] = find(fi_diff);

for sci = 1:3
    norm = 1;
    avgnorm = 1;
    i = 0;
    mu = zeros(nobs_c, ns);
    p_e = p_s_0;
    if sci == 1
        subsidy_e = subsidy_h;
    elseif sci == 2
        subsidy_e = subsidy_l;
    elseif sci == 3
        subsidy_e = sub_tt0;
    end
    pd_i = zeros(size(temp,1),ns); 
    while norm > 1e-1*tol*10^(floor(i/50)) & avgnorm > 1e-3*tol*10^(floor(i/50))
        p_d = (p_e-subsidy_e).*tts + qcost;
        if sci == 3
            mval = zeros(size(temp,1),ns);
            for ri = 1:ns
                v_i = vi(:,ri:ns:k*ns);
                d_i = di(:,ri:ns:nj*ns);

                pd_i(:,ri) = p_d + subdec(:,ri).*tts;
                mval(:,ri) = [log(pd_i(:,ri)) x1e]*theta1 + ehat;
                mu(:,ri) = ([log(pd_i(:,ri)),x2e(:,2:end)].*v_i*theta2w(:,1))...
                    +[log(pd_i(:,ri)),x2e(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
            end
            eg = exp(mu + mval);
        else
            mval = [log(p_d) x1e]*theta1 + ehat;
            for ri = 1:ns
                v_i = vi(:,ri:ns:k*ns);
                d_i = di(:,ri:ns:nj*ns);
                mu(:,ri) = ([log(p_d),x2e(:,2:end)].*v_i*theta2w(:,1))...
                    +[log(p_d),x2e(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
            end
            eg = exp(mu + kron(ones(1,ns),mval));
        end
        
        % Calculate market share for each city separately
        s_i = zeros(size(eg));
        
        hy_omega = zeros(max(unique(mid)));  
        hy_s = zeros(max(unique(mid)),1);
        
        for t = min(cd_range):max(cd_range)
            trows = find(cdid_hy == t);
            time_trows = trows(hy_temp(trows) == hy_given);
            
            % Calculate market share for this city
            s_i(time_trows,:) = eg(time_trows,:)./(ones(length(time_trows),1)*(1 + sum(eg(time_trows,:))));
            
            % Calculate omega matrix and aggregate
            tempn = size(time_trows,1);
            
            if sci == 3
                temp0 = (s_i(time_trows,:).*alpha_i(time_trows,:)./pd_i(time_trows,:).*tts(time_trows))*s_i(time_trows,:)'/ns;
                temp1 = sum((s_i(time_trows,:).*alpha_i(time_trows,:)./pd_i(time_trows,:).*tts(time_trows))')'/ns;
                H = (diag(temp1) - temp0);
            else
                temp0 = (s_i(time_trows,:).*alpha_i(time_trows,:))*s_i(time_trows,:)'/ns;
                temp1 = sum((s_i(time_trows,:).*alpha_i(time_trows,:))')'/ns;
                H = (diag(temp1) - temp0)./(p_d(time_trows,:)*ones(1,tempn)).*tts(time_trows);
            end
                        
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
        
        if sci == 2
            s_i_l = s_i; % save s_i for later drawing income vs ev adoption rate
        end

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
%         pnew = (p_e-subsidy_e).*tts + qcost;
        norm = max(abs(pnew - p_e));
        avgnorm = mean(abs(pnew - p_e));
        p_e = pnew;
        i = i + 1;
    end
    
    if sci == 3
        p_d00 = p_d;
        p_d = mean(pd_i,2);
    end
    pd_mat(:,sci+1) = p_d;
    ps_mat(:,sci+1) = p_e;
    avgprice_d( sci+1,:) = [mean(p_d(find(ev_tt == 1 & import_tt == 0))),...
        mean(p_d(find(fv_tt == 1& import_tt == 0))),...
        mean(p_d(find(ev_tt == 1 & import_tt == 1))),...
        mean(p_d(find(fv_tt == 1& import_tt == 1)))];
    avgprice_s( sci+1,:) = [mean(p_e(find(ev_tt == 1 & import_tt == 0))),...
        mean(p_e(find(fv_tt == 1 & import_tt == 0))),...
        mean(p_e(find(ev_tt == 1 & import_tt == 1))),...
        mean(p_e(find(fv_tt == 1 & import_tt == 1)))];
    
    % Calculate sales and sub_inc for each city separately
    sub_inc_city = zeros(length(unique(cdid_hy)), ns);
    sales_ev0 = 0; sales_fv0 = 0; sales_ev1 = 0; sales_fv1 = 0;
    city_idx = 0;
    for t = min(cd_range):max(cd_range)
        trows = find(cdid_hy == t);
        time_trows = trows(hy_temp(trows) == hy_given);
        city_mktsize = unique(mktsizets_hy(time_trows));
        
        city_idx = city_idx + 1;
        
        % Calculate sub_inc for this city
        s_i_city = s_i(time_trows,:);
        sub_tt0_city = sub_tt0(time_trows);
        subsidy_e_city = subsidy_e(time_trows);
        sub_idx = (sub_tt0_city > 0);
        if sci ~= 3
            sub_inc_city(city_idx,:) = sum(s_i_city(sub_idx, :).*...
                (subsidy_e_city(sub_idx)*ones(1, ns)))./sum(s_i_city(sub_idx, :));
        else
            subdec_city = subdec(time_trows,:);
            sub_inc_city(city_idx,:) = sum(s_i_city(sub_idx, :).*...
                (subsidy_e_city(sub_idx)*ones(1, ns) - subdec_city(sub_idx,:)))./sum(s_i_city(sub_idx, :));
        end
        
        % Calculate sales for this city
        sales_ev0 = sales_ev0 + sum(sj_e(time_trows(find(ev_tt(time_trows) == 1 & import_tt(time_trows) == 0))))*city_mktsize;
        sales_fv0 = sales_fv0 + sum(sj_e(time_trows(find(fv_tt(time_trows) == 1 & import_tt(time_trows) == 0))))*city_mktsize;
        sales_ev1 = sales_ev1 + sum(sj_e(time_trows(find(ev_tt(time_trows) == 1 & import_tt(time_trows) == 1))))*city_mktsize;
        sales_fv1 = sales_fv1 + sum(sj_e(time_trows(find(fv_tt(time_trows) == 1 & import_tt(time_trows) == 1))))*city_mktsize;
    end
    sub_inc(sci + 1, :) = reshape(sub_inc_city(1:city_idx, :)', 1, []);
    sales_sim(sci+1,:) = [sales_ev0, sales_fv0, sales_ev1, sales_fv1];
    %     shr_e = [shr_e; sj_e]; price_e = [price_e; p_e];
    mktshr_e(:,sci) = sj_e;
    
    % Calculate profit for each firm across all cities
    % Get unique firm IDs
    unique_firms = unique(f_i);
    [unique_firms,ia,ic] = unique(f_i);
    n_firms = length(unique_firms);
    pi_firm = zeros(n_firms, 5);
    
    % For each firm, aggregate profit across all cities
    for firm_idx = 1:n_firms
        firm_id = unique_firms(firm_idx);
        firm_products = find(f_i == firm_id);
        n_products = sum(f_i == firm_id);

        % Initialize profit for this firm
        pi_ev = 0; pi_import = 0; pi_orig = 0; pi_counter = 0;
        
        % Aggregate profit from all cities
        for t = min(cd_range):max(cd_range)
            trows = find(cdid_hy == t);
            time_trows = trows(hy_temp(trows) == hy_given);
            
            % Get products of this firm in this city
            firm_products_city = intersect(firm_products, time_trows);
            
            if ~isempty(firm_products_city)
                % Calculate profit for this firm in this city
                % omega_sta{t} aggregates profit within a firm in this city
                % All products of the same firm get the same profit value
                % So we just need to take the first product's value
                first_product = firm_products_city(1);
                idx_in_city = find(time_trows == first_product);
                
                % Calculate aggregated values for this city
                ev_city = omega_sta{t}(idx_in_city, :)*ev_tt(time_trows, :);
                import_city = omega_sta{t}(idx_in_city, :)*import_tt(time_trows, :);
                pi_orig_city = omega_sta{t}(idx_in_city, :)*(sj_e0(time_trows, :).*mkup_e0(time_trows, :).*mktsizets(time_trows, :));
                pi_counter_city = omega_sta{t}(idx_in_city, :)*(sj_e(time_trows, :).*mkup_e(time_trows, :).*mktsizets(time_trows, :));
                
                % Add to firm's total profit
                pi_ev = pi_ev + ev_city;
                pi_import = pi_import + import_city;
                pi_orig = pi_orig + pi_orig_city;
                pi_counter = pi_counter + pi_counter_city;
            end
        end
        
        pi_firm(firm_idx, :) = [pi_ev, pi_import, pi_orig, pi_counter, n_products];
    end
    t_evd = find((pi_firm(:,1)>0 & pi_firm(:,5) == pi_firm(:,1)) & pi_firm(:,2) == 0);
    % the first two columns are profits of domestic ev manufacturers
    % in original and counterfactual scenario
    pi_category(sci, 1:2) = sum(pi_firm(t_evd,3:4), 1);
    % the second two columns are profits of domestic fv manufacturers
    % in original and counterfactual scenario
    t_fvd = find(pi_firm(:,1)== 0  & pi_firm(:,2) == 0);
    pi_category(sci, 3:4) = sum(pi_firm(t_fvd,3:4), 1);
    % the third two columns are profits of domestic hybrid manufacturers
    % in original and counterfactual scenario
    t_hybrd = find((pi_firm(:,1)>0 & pi_firm(:,5) > pi_firm(:,1)) & pi_firm(:,2) == 0);
    pi_category(sci, 5:6) = sum(pi_firm(t_hybrd,3:4), 1);

    t_evi = find((pi_firm(:,1)>0 & pi_firm(:,5) == pi_firm(:,1)) & pi_firm(:,2) > 0);
    % the fourth two columns are profits of imported ev manufacturers
    % in original and counterfactual scenario
    pi_category(sci, 7:8) = sum(pi_firm(t_evi,3:4), 1);
    % the fifth two columns are profits of imported fv manufacturers
    % in original and counterfactual scenario
    t_fvi = find(pi_firm(:,1)== 0  & pi_firm(:,2) > 0);
    pi_category(sci, 9:10) = sum(pi_firm(t_fvi, 3:4), 1);
    % the sixth two columns are profits of imported hybrid manufacturers
    % in original and counterfactual scenario
    t_hybri = find((pi_firm(:,1)>0 & pi_firm(:,5) > pi_firm(:,1)) & pi_firm(:,2) > 0);
    pi_category(sci, 11:12) = sum(pi_firm(t_hybri,3:4), 1);

    cvs_city = zeros(length(unique(cdid_hy)),ns);
    cvs_city_mean = zeros(length(unique(cdid_hy)),1);
    city_idx = 0;
    for t = min(cd_range):max(cd_range)
        trows = find(cdid_hy == t);
        time_trows = trows(hy_temp(trows) == hy_given);
        
        city_idx = city_idx + 1;
        city_mktsize = unique(mktsizets_hy(time_trows));
        
        % Calculate CV for this city
        if sci~=3
            cvs_city(city_idx,:) = -city_mktsize*mean((log(1 + sum(eg(time_trows,:))) - log(1 + sum(eg0(time_trows,:))))./...
                         (alpha_i(time_trows,:)./p_d(time_trows)) );
        else
            cvs_city(city_idx,:) = -city_mktsize*mean((log(1 + sum(eg(time_trows,:))) - log(1 + sum(eg0(time_trows,:))))./...
                         (alpha_i(time_trows,:)./pd_i(time_trows)) );
        end
        cvs_city_mean(city_idx) = mean(cvs_city(city_idx,:),2);
    end
    
    % Total CV (sum of all cities' CV) 
    cvs(sci,1) = sum(cvs_city_mean);

    % subsidy calculation
    if sci == 3
        subsidy_c(sci+1,1) = sum((subsidy_e - mean(subdec, 2)).*sj_e.*mktsizets);
    else
        subsidy_c(sci+1,1) = sum(subsidy_e.*sj_e.*mktsizets);
    end
    pr_sim(:,sci) = p_e;
    ex_fv(sci+1,1) = sum(sj_e.*mktsizets.*fc_fv(temp).*(1-ev_tt)*mc_ex_fv*lifetime*vmt);
    ex_ev(sci+1,1) = sum(sj_e.*mktsizets.*fc_ev(temp).*ev_tt*mc_ex_ev*lifetime*vmt);

    ex_evd(sci+1,1) = sum(sj_e.*mktsizets.*fc_ev(temp).*ev_tt*mc_ex_ev2*lifetime*vmt);
    disp(['# of iterations for price mapping    ' num2str(i)])
    clear temp0 temp1 mu
end

%% Optimize and evaluate the alternative subsidy rule
%   Estimate the coefficients of the alternative income-and-range subsidy rule,
%   resolve the resulting equilibrium, and append its policy outcomes to the
%   benchmark scenario results.

subsidy_0 = subsidy_c(1);
population = data.population;
model_variant = data.model_variant;
save opt_subsidy tol temp theta2w theta1 nj subsidy_0 sub_tt0 p_s_0 hy cd_range population model_variant
LB = [0, -Inf, -Inf]';
UB = [Inf, 0, Inf]';
opt_con = optimset( 'TolFun', 1e-7, 'TolX', 1e-6);
[optsb_coeff,FVAL,EXITFLAG] = fmincon('optsb2', [0.2, -.4, .005]', [],[],[],[],LB, UB,[],opt_con);

[mu, mval, pd_i, p_e, sj_e, mkup_e, subsidy_total_al, rg_w, sub_ind, cv_opt] = optsubsidy2(optsb_coeff);
sci = 4;
cvs(sci,1) = mean(cv_opt); % in 10 billion
sub_inc(sci + 1, :) = sub_ind;
p_d = mean(pd_i,2);
    pd_mat(:,sci+1) = p_d;
    ps_mat(:,sci+1) = p_e;

avgprice_d( sci+1,:) = [mean(p_d(find(ev_tt == 1 & import_tt == 0))),...
    mean(p_d(find(fv_tt == 1& import_tt == 0))),...
    mean(p_d(find(ev_tt == 1 & import_tt == 1))),...
    mean(p_d(find(fv_tt == 1& import_tt == 1)))];
avgprice_s( sci+1,:) = [mean(p_e(find(ev_tt == 1 & import_tt == 0))),...
    mean(p_e(find(fv_tt == 1& import_tt == 0))),...
    mean(p_e(find(ev_tt == 1 & import_tt == 1))),...
    mean(p_e(find(fv_tt == 1& import_tt == 1)))];

% Calculate sales for each city separately and sum
sales_ev0 = 0; sales_fv0 = 0; sales_ev1 = 0; sales_fv1 = 0;
for t = min(cd_range):max(cd_range)
    trows = find(cdid_hy == t);
    time_trows = trows(hy_temp(trows) == hy_given);
    city_mktsize = unique(mktsizets_hy(time_trows));
    
    sales_ev0 = sales_ev0 + sum(sj_e(time_trows(find(ev_tt(time_trows) == 1 & import_tt(time_trows) == 0))))*city_mktsize;
    sales_fv0 = sales_fv0 + sum(sj_e(time_trows(find(fv_tt(time_trows) == 1 & import_tt(time_trows) == 0))))*city_mktsize;
    sales_ev1 = sales_ev1 + sum(sj_e(time_trows(find(ev_tt(time_trows) == 1 & import_tt(time_trows) == 1))))*city_mktsize;
    sales_fv1 = sales_fv1 + sum(sj_e(time_trows(find(fv_tt(time_trows) == 1 & import_tt(time_trows) == 1))))*city_mktsize;
end
sales_sim(sci+1,:) = [sales_ev0, sales_fv0, sales_ev1, sales_fv1];
mktshr_e(:,sci) = sj_e;

% Calculate profit for each firm across all cities
% Get unique firm IDs
unique_firms = unique(f_i);
n_firms = length(unique_firms);
pi_firm = zeros(n_firms, 5);

% For each firm, aggregate profit across all cities
for firm_idx = 1:n_firms
    firm_id = unique_firms(firm_idx);
    firm_products = find(f_i == firm_id);
    n_products = sum(f_i == firm_id);
    
    % Initialize profit for this firm
    pi_ev = 0; pi_import = 0; pi_orig = 0; pi_counter = 0;
    
    % Aggregate profit from all cities
    for t = min(cd_range):max(cd_range)
        trows = find(cdid_hy == t);
        time_trows = trows(hy_temp(trows) == hy_given);
        
        % Get products of this firm in this city
        firm_products_city = intersect(firm_products, time_trows);
        
        if ~isempty(firm_products_city)
            % Calculate profit for this firm in this city
            % omega_sta{t} aggregates profit within a firm in this city
            % All products of the same firm get the same profit value
            % So we just need to take the first product's value
            first_product = firm_products_city(1);
            idx_in_city = find(time_trows == first_product);
            
            % Calculate aggregated values for this city
            ev_city = omega_sta{t}(idx_in_city, :)*ev_tt(time_trows, :);
            import_city = omega_sta{t}(idx_in_city, :)*import_tt(time_trows, :);
            pi_orig_city = omega_sta{t}(idx_in_city, :)*(sj_e0(time_trows, :).*mkup_e0(time_trows, :).*mktsizets(time_trows, :));
            pi_counter_city = omega_sta{t}(idx_in_city, :)*(sj_e(time_trows, :).*mkup_e(time_trows, :).*mktsizets(time_trows, :));
            
            % Add to firm's total profit
            pi_ev = pi_ev + ev_city;
            pi_import = pi_import + import_city;
            pi_orig = pi_orig + pi_orig_city;
            pi_counter = pi_counter + pi_counter_city;
        end
    end
    
    pi_firm(firm_idx, :) = [pi_ev, pi_import, pi_orig, pi_counter, n_products];
end
t_evd = find((pi_firm(:,1)>0 & pi_firm(:,5) == pi_firm(:,1)) & pi_firm(:,2) == 0);
% the first two columns are profits of domestic ev manufacturers
% in original and counterfactual scenario
pi_category( sci, 1:2) = sum(pi_firm(t_evd,3:4), 1);
% the second two columns are profits of domestic fv manufacturers
% in original and counterfactual scenario
t_fvd = find(pi_firm(:,1)== 0 &  pi_firm(:,2) == 0);
pi_category(sci, 3:4) = sum(pi_firm(t_fvd,3:4), 1);
% the third two columns are profits of domestic hybrid manufacturers
% in original and counterfactual scenario
t_hybrd = find((pi_firm(:,1)>0 & pi_firm(:,5) > pi_firm(:,1)) & pi_firm(:,2) == 0);
pi_category(sci, 5:6) = sum(pi_firm(t_hybrd,3:4), 1);

t_evi = find((pi_firm(:,1)>0  & pi_firm(:,5) == pi_firm(:,1)) & pi_firm(:,2) > 0);
% the fourth two columns are profits of imported ev manufacturers
% in original and counterfactual scenario
pi_category(sci, 7:8) = sum(pi_firm(t_evi,3:4), 1);
% the fifth two columns are profits of imported fv manufacturers
% in original and counterfactual scenario
t_fvi = find(pi_firm(:,1)== 0  & pi_firm(:,2) > 0);
pi_category(sci, 9:10) = sum(pi_firm(t_fvi, 3:4), 1);
% the sixth two columns are profits of imported hybrid manufacturers
% in original and counterfactual scenario
t_hybri = find((pi_firm(:,1)>0 & pi_firm(:,5) > pi_firm(:,1)) & pi_firm(:,2) > 0);
pi_category(sci, 11:12) = sum(pi_firm(t_hybri,3:4), 1);
subsidy_c(sci+1,1) = subsidy_total_al;
ex_fv(sci+1,1) = sum(sj_e.*mktsizets.*fc_fv(temp).*(1-ev_tt)*mc_ex_fv*lifetime*vmt);
ex_ev(sci+1,1) = sum(sj_e.*mktsizets.*fc_ev(temp).*ev_tt*mc_ex_ev*lifetime*vmt);
ex_evd(sci+1,1) = sum(sj_e.*mktsizets.*fc_ev(temp).*ev_tt*mc_ex_ev2*lifetime*vmt);

%% Report counterfactual simulation results
%   Display and log the principal price, sales, welfare, profit, subsidy, and
%   externality outcomes across the benchmark and optimized policy scenarios.

diary on
disp(['average price'])
avgprice_d
avgprice_s

disp(['total sales'])
sales_sim/1000
disp(['compensating variation in billion'])
cvs/10^5
disp(['profits in billion'])
pi_category/10^5
disp(['subsidy in billion'])
subsidy_c/10^5
disp(['externality'])
[ex_ev,  ex_fv]/10^9
[ex_evd,  ex_fv]/10^9
diary off

%% Paper Figure 9: Subsidy Distribution over Vehicle Ranges
% Replication output:
%   range_dist.eps (Figure 9 in the paper).
% Figure content:
%   Share-weighted densities of EV range under the status quo and the
%   four counterfactual policy scenarios reported in the paper.
range_e = evrange(temp);
range_eev = range_e(find(ev_tt == 1));
sj_e0ev = sj_e0(find(ev_tt == 1));
sj_e1ev = mktshr_e(find(ev_tt == 1),1);
sj_e2ev = mktshr_e(find(ev_tt == 1), 2);
sj_e3ev = mktshr_e(find(ev_tt == 1), 3);

pts = (min(range_eev):10:max(range_eev));
[F0, xi] = ksdensity(range_eev,pts, 'Weights', sj_e0ev, 'Support', 'positive');
[F1, xi] = ksdensity(range_eev,pts, 'Weights', sj_e1ev, 'Support', 'positive');
[F2, xi] = ksdensity(range_eev,pts, 'Weights', sj_e2ev, 'Support', 'positive');
[F3, xi] = ksdensity(range_eev,pts, 'Weights', sj_e3ev, 'Support', 'positive');
[F4, xi] = ksdensity(range_eev,pts, 'Weights', sj_e(ev_tt  ==1), 'Support', 'positive');
plot(xi,F0, 'LineStyle', '-')
hold on
plot(xi, F1, 'LineStyle', ':')
plot(xi, F2, 'LineStyle', '--')
plot(xi, F3, 'LineStyle', '-.')
plot(xi, F4, 'LineStyle', '-')
legend('Scenario Null','Scenario (1)','Scenario (2)','Scenario (3)', 'Scenario (4)')
title('')
xlabel('Range')
ylabel('Density')
saveas(gcf,'range_dist','epsc')
hold off
close(gcf)

%% Analyze subsidy pass-through
%   Compare supply and consumer prices across policy scenarios and summarize
%   subsidy pass-through by producer type and price quantile.

price_supply=[p_s_0, pr_sim];
subsidy_summarize = [sub_tt0, subsidy_h, subsidy_l, sub_tt0 ];
price_demand = (price_supply-subsidy_summarize).*tts + qcost;
pass_through_demand =(price_demand < price_supply).* abs(  price_demand - price_demand(:,3))./subsidy_summarize;
pass_through_demand(isnan(pass_through_demand)) = 0;
subsidy_sci4 = p_e - (pd_i - qcost)./tts ;
subsidy_sci4 = subsidy_sci4.* (sub_tt0>0);
pass_through_demandsci4 =abs(  pd_i - price_demand(:,3))./subsidy_sci4;
pass_through_demandsci4(isinf(pass_through_demandsci4)) = 0;
pass_through_demand = [ pass_through_demand, mean(pass_through_demandsci4,2)];
evfirm = find(ismember(fnmb(temp), unique_firms(t_evd)));
evfvfirm = find(ismember(fnmb(temp),unique_firms([t_evi;t_fvd;t_hybrd;t_fvi;t_hybri])));
pass_through_rate_demand=[sum(pass_through_demand(evfirm,:))./sum(pass_through_demand(evfirm,:)~=0);sum(pass_through_demand(evfvfirm,:))./sum(pass_through_demand(evfvfirm,:)~=0);sum(pass_through_demand)./sum(pass_through_demand~=0)];
pass_through_rate_demand(:,all(isnan(pass_through_rate_demand),1))=[];
diary on 
disp(['average price'])
disp('Pass-through rate to Consumers')
pass_through_rate_demand
diary off

pt_target = [price_demand(:,1) pass_through_demand(:,1)];
pt_target(pt_target(:,2)==0,:) = [];
pt_sort = sort(pt_target,1);
pt_prct = prctile(pt_sort(:,1), [25,50,75] );
pt_quantile = [ mean( pt_target(pt_target(:,1) <= pt_prct(1) ,2)), ...
  mean( pt_target(pt_prct(2) >= pt_target(:,1) &  pt_prct(1) < pt_target(:,1) ,2)), ...
  mean( pt_target(pt_prct(3) >= pt_target(:,1) &  pt_prct(2) < pt_target(:,1) ,2)), ... 
  mean( pt_target(pt_target(:,1) > pt_prct(3) ,2)) ];
diary on
disp('Pass-through quantile')
pt_quantile
diary off






%% Analyze subsidy incidence across income
%   Match simulated incomes to scenario-specific individual subsidies and
%   summarize the distribution of subsidy benefits across income bins.


[unique_rows, ia,ic] = unique(cdid_hy, 'rows', 'first');
SMat_tem = [reshape(di(ia,:)'*100/2.7,1,[])',sub_inc'];
SMat = sortrows(SMat_tem, 1);
clear SMat_tem

% since the income distribution is log normal, there are some extremely
% high income, we just get rid of points with individual income higher than 350 of data points
%% Paper Figure 8: Subsidy Distribution over Income
% Replication output:
%   subsidy_dist.eps (Figure 8 in the paper).
% Figure content:
%   Mean expected subsidy across income bins for the status quo, high-subsidy,
%   and progressive-subsidy scenarios.

inc_threshold = 350;
SMat = SMat(SMat(:,1)<inc_threshold,:);
xi = SMat(:,1); 
xi_min = min(xi);
inc_gap = max(xi) - xi_min;
inc_cat_n = floor(inc_gap/10) + 1;
inc_cat = zeros(inc_cat_n,1); 
Sb_m = zeros(inc_cat_n, 4);
for i = 1:inc_cat_n - 1
    tempxi = find(xi >= xi_min + (i-1)*10 & xi < xi_min + i*10 );
    if isempty(tempxi) == 0
    inc_cat(i,1) = xi_min + (i-1)*10;
    Sb_m(i,:) = mean(SMat(tempxi,[2,3,5,6]),1);
    end
end
tempxi = find(xi>= xi_min + i*10);
Sb_m(i+1,:) = mean(SMat(tempxi,[2,3,5,6]),1);
inc_cat(i+1,1) = xi_min + i*10;

% remove zeros
inc_cat(~any(inc_cat,2),:) = [];
Sb_m( ~any(Sb_m, 2), :) = [];
plot(inc_cat,Sb_m(:,1), 'LineStyle', '-')
hold on
% plot(inc_cat, Sb_m(:,3), 'LineStyle', ':')
plot(inc_cat, Sb_m(:,2), 'LineStyle', '--')
plot(inc_cat, Sb_m(:,3), 'LineStyle', '-.')
plot(inc_cat, Sb_m(:,4), 'LineStyle', '-')
legend('Scenario Null','High Subsidy','Progressive Subsidy (1)', 'Progressive Subsidy (2)')
legend('Location','northwest')
title('')
xlabel('Income')
ylabel('Expected Subsidy')
saveas(gcf,'subsidy_dist','epsc')
hold off
close(gcf)

tempincr1 = find(SMat(:,1)<= 60); 
tempincr2 = find(SMat(:,1)> 60 & SMat(:,1)<= 90); 
tempincr3 = find(SMat(:,1)> 90 & SMat(:,1) <= 144);
tempincr4 = find(SMat(:,1)> 144);
avg_sub = [mean(SMat(tempincr1, [2:6])); 
mean(SMat(tempincr2, [2:6]));
mean(SMat(tempincr3, [2:6]));
mean(SMat(tempincr4, [2:6]))];
clear -regexp ^tempincr
inc_int = ['Less than  36,000';
           '36,000  ~  96,000';
           '96,000 ~  144,000';
           'more than 144,000'];
diary on
disp(['Income Intervals', 'Null', '(i)', '(ii)', '(iii)', '(iv)'])
for i = 1:4 
    disp([inc_int(i,:), num2str(avg_sub(i,:))])
end
diary off
save(fullfile(diary_path, ['results_' RUN_TAG '.mat']), '-v7.3')


%% Analyze EV adoption across income
%   Compare simulated EV-adoption rates across income quantiles under the
%   observed-policy and no-subsidy cases.

% status quo
s_i_0ev = zeros(length(cd_range),ns);
% no subsidy case
s_i_lev = zeros(length(cd_range),ns);
% we have to guarantee cd_range contains only consecutive rows
for i = 1:length(cd_range)
    t = cd_range(i);
    trows = find(cdid_hy == t);
    time_trows = trows(hy_temp(trows) == hy_given);
    s_i_0ev(i,:) = sum(s_i_0(time_trows,:).*repmat(ev_tt(time_trows),1,ns),1);
    s_i_lev(i,:) = sum(s_i_l(time_trows,:).*repmat(ev_tt(time_trows),1,ns),1);
end
s_i0ev_row = reshape(s_i_0ev,[],1);
s_ilev_row = reshape(s_i_lev,[],1);
income_row = reshape(demogr(cd_range,:),[],1);
% sort data: status quo
sorted_data = sortrows([income_row, s_i0ev_row]);
income_sorted = sorted_data(:, 1);
s_i0ev_sorted = sorted_data(:, 2);
n_bins = 10;
[Income_Mean0, EV_Adoption_Mean0,EV_Adoption_Std0,...
    GroupCount0] = drawincome_ev(income_sorted, s_i0ev_sorted, n_bins);
% sort data: no subsidy case
sorted_datal = sortrows([income_row, s_ilev_row]);
income_sorted = sorted_datal(:, 1);
s_ilev_sorted = sorted_datal(:, 2);
n_bins = 10;
[Income_Meanl, EV_Adoption_Meanl,EV_Adoption_Stdl,...
    GroupCountl] = drawincome_ev(income_sorted, s_ilev_sorted, n_bins);

%% Paper Figure 6: Income Distribution (Simulated) and EV Adoption Rate
% Replication output:
%   EVadop_income.eps (Figure 6 in the paper).
% Figure content:
%   Mean EV adoption rates under the status quo and no-subsidy scenario across
%   ten household-income quantile bins.

figure(1);
hold on
h1 = plot(Income_Mean0, EV_Adoption_Mean0, 'o-', ...
    'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'b',...
    'Color','b','DisplayName', 'Observed');

h2 = plot(Income_Meanl, EV_Adoption_Meanl, 's-', ...
    'LineWidth', 1, 'MarkerSize', 4, 'MarkerFaceColor', 'r',...
    'Color','r','DisplayName', 'No Subsidy Case');

legend([h1, h2],{'Observed', 'No Subsidy Case'}, ...
    'Location', 'northwest', 'FontSize', 10);
xlabel('Mean Income in Bin');
ylabel('EV Adoption Rate');
%title('EV Adoption Rate by Income Quantiles');
saveas(gcf,'EVadop_income','epsc')
hold off
close(gcf)

sum(sj_e0(cdid_hy == cd_range(2)))
sum(s_jt(etc))



%% Extend equilibrium and pass-through calculations across time
% Purpose:
%   Repeat the relevant equilibrium calculations for every half-year and store
%   period-specific prices, subsidies, and consumer pass-through measures.
% Produces:
%   Half-year cell arrays for equilibrium prices, subsidies, tax adjustments,
%   quota costs, and pass-through measures.
% notice that after running this part, you should not rerun any code before
% this part!


% save current p_s_0 for checking
p_s_0_hy8 = p_s_0;
pass_through_demand_hy8 = pass_through_demand;
% Initialize storage for equilibrium prices across all half-years
p_s_0_all = cell(8, 1);
pr_sim_all = cell(8, 1);
subsidy_l_all = cell(8, 1);
sub_tt0_all = cell(8, 1);
tts_all = cell(8, 1);
qcost_all = cell(8, 1);
sc_tt_all = cell(8, 1);

% Loop through all half-years (hy_given = 1 to 8)
for hy_iter = 1:8
    t = hy_iter; % Current half-year

    disp(["current hy is: ", t]);
    
    % Get data for this half-year
    temp_iter = find(data.hy == t);
    nobs_c_iter = length(temp_iter);
    
    if nobs_c_iter == 0
        continue;
    end
    
    % Extract variables for this half-year
    p_e_iter = price(temp_iter);
    range_tt_iter = evrange(temp_iter);
    ev_tt_iter = ev(temp_iter);
    import_tt_iter = import(temp_iter);
    fv_tt_iter = 1 - ev_tt_iter;
    
    % Calculate temp_sc using the existing switch-case logic
    temp_sc = zeros(nobs_c_iter, 1);
    switch true
        case t == 1 || t == 2
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter <100) = 0;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=100 & range_tt_iter < 150) = 2.5;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=150 & range_tt_iter < 250) = 4.5;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=250) = 5.5;
        case t == 3 || t == 4 || t == 5
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter <100) = 0;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=100 & range_tt_iter < 150) = 2;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=150 & range_tt_iter < 250) = 3.6;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=250) = 4.4;
        case t == 6 || t == 7
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter <150) = 0;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=150 & range_tt_iter < 200) = 1.5;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=200 & range_tt_iter < 250) = 2.4;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=250 & range_tt_iter < 300) = 3.4;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=300 & range_tt_iter < 400) = 4.5;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=400) = 5;
        case t == 8
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter <250) = 0;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=250 & range_tt_iter < 400) = 1.8;
            temp_sc(import_tt_iter == 0 & ev_tt_iter == 1 & range_tt_iter >=400) = 2.5;
        otherwise
            error("unknown halfyear")
    end
    sc_tt = min([temp_sc, p_e_iter*.6],[],2);
    
    % Calculate subsidy variables
    local_tt_iter = subsidy_local_net(temp_iter);
    sub_tt0 = sc_tt + local_tt_iter;
    
    subsidy_l = zeros(size(temp_iter));
    
    % Calculate tts and qcost
    qcost = quota_cost(temp_iter).*(ev(temp_iter) == 0);
    tts = (1 + 1/1.13*.1*(ev(temp_iter) == 0));
    
    % Get other necessary variables
    cost_sim_iter = mc(temp_iter);
    ehat_iter = gmmresid(temp_iter,:);
    mktsizets_iter = mktsize(temp_iter);
    x1e_iter = x1(temp_iter,:);
    x2e_iter = x2(temp_iter,:);
    f_i_iter = fnmb(temp_iter,:);
    cdid_hy_iter = cdid(temp_iter);
    hy_temp_iter = hy(temp_iter);
    mid_iter = data.model_variant(temp_iter);
    
    % Build omega_sta for this half-year
    cd_range_iter = unique(cdid(find(hy==t)));
    omega_sta_iter = cell(length(unique(cdid_hy_iter)),1);
    for city_t = min(cd_range_iter):max(cd_range_iter)
        trows = find(cdid_hy_iter == city_t);
        tempn = size(trows,1);
        time_trows = trows(hy_temp_iter(trows) == t);
        omega_sta_c = zeros(tempn);
        f_i_city = f_i_iter(time_trows,:);
        for j = 1:tempn
            for h = 1:tempn
                omega_sta_c(j,h) = (f_i_city(j) == f_i_city(h));
            end
        end
        omega_sta_iter{city_t} = omega_sta_c;
    end
    
    vi_iter = vfull(temp_iter,:);
    di_iter = exp(dfull(temp_iter, :));
    
    demo_i_iter = zeros(size(temp_iter,1), ns);
    for n = 1:nj
        demo_i_iter = di_iter(:, (n-1)*ns + 1:n*ns)*betapi(n,:) + demo_i_iter;
    end
    alpha_i_iter = theta1(1) + vi_iter(:,1:ns)*sigmap + demo_i_iter;
    
    % Calculate observed equilibrium (status quo)
    p_s_o = price(temp_iter) + sc_tt;
    norm = 1;
    avgnorm = 1;
    i_iter = 0;
    mu_iter = [];
    p_e = p_s_o;
    
    while norm > 1e-1*tol*10^(floor(i_iter/50)) & avgnorm > 1e-3*tol*10^(floor(i_iter/50))
        p_d = (p_e - sub_tt0).*tts + qcost;
        mval = [log(p_d), x1e_iter]*theta1 + ehat_iter;
        for ri = 1:ns
            v_i = vi_iter(:,ri:ns:k*ns);
            d_i = di_iter(:,ri:ns:nj*ns);
            mu_iter(:,ri) = ([log(p_d),x2e_iter(:,2:end)].*v_i*theta2w(:,1))+[log(p_d),x2e_iter(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
        end
        MU0 = mu_iter + kron(ones(1,ns),mval);
        eg0 = exp(MU0);
        
        s_i = zeros(size(eg0));
        hy_omega = zeros(max(unique(mid_iter)));
        hy_s = zeros(max(unique(mid_iter)),1);
        
        for city_t = min(cd_range_iter):max(cd_range_iter)
            trows = find(cdid_hy_iter == city_t);
            time_trows = trows(hy_temp_iter(trows) == t);
            
            s_i(time_trows,:) = eg0(time_trows,:)./(ones(length(time_trows),1)*(1 + sum(eg0(time_trows,:))));
            
            tempn = size(time_trows,1);
            temp0 = (s_i(time_trows,:).*alpha_i_iter(time_trows,:))*s_i(time_trows,:)'/ns;
            temp1 = sum((s_i(time_trows,:).*alpha_i_iter(time_trows,:))')'/ns;
            H = (diag(temp1) - temp0)./(p_d(time_trows,:)*ones(1,tempn)).*tts(time_trows,:);
            
            omega = H.*omega_sta_iter{city_t};
            
            modelid = mid_iter(time_trows);
            city_mktsize = unique(mktsizets_iter(time_trows));
            if length(city_mktsize) ~= 1, error('Market size not constant within city'); end
            s_j_city = mean(s_i(time_trows,:), 2);
            
            for row = 1:length(time_trows)
                for col = 1:length(time_trows)
                    hy_omega(modelid(row),modelid(col)) = hy_omega(modelid(row),modelid(col)) + omega(row,col)*city_mktsize;
                end
                hy_s(modelid(row),1) = hy_s(modelid(row),1) + s_j_city(row)*city_mktsize;
            end
        end
        
        nanproduct = find(hy_s==0);
        hy_s(nanproduct)=[];
        hy_omega(nanproduct,:)=[];
        hy_omega(:,nanproduct)=[];
        
        mkup_model = -inv(hy_omega)*hy_s;
        
        temp_m = [1:max(unique(mid_iter))]';
        temp_m(nanproduct)=[];
        
        mkup_e0 = zeros(length(p_e),1);
        for index = 1:length(temp_m)
            mkup_e0(find(mid_iter==temp_m(index))) = mkup_model(index);
        end
        
        pnew = mkup_e0 + cost_sim_iter;
        norm = max(abs(pnew - p_e));
        avgnorm = mean(abs(pnew - p_e));
        p_e = pnew;
        i_iter = i_iter + 1;
    end
    
    p_s_0 = p_e; % Observed equilibrium price
    
    % Calculate sci == 2 (no subsidy) equilibrium
    norm = 1;
    avgnorm = 1;
    i_iter = 0;
    mu_iter = zeros(nobs_c_iter, ns);
    p_e = p_s_0;
    subsidy_e = subsidy_l; % No subsidy
    
    while norm > 1e-1*tol*10^(floor(i_iter/50)) & avgnorm > 1e-3*tol*10^(floor(i_iter/50))
        p_d = (p_e - subsidy_e).*tts + qcost;
        mval = [log(p_d) x1e_iter]*theta1 + ehat_iter;
        for ri = 1:ns
            v_i = vi_iter(:,ri:ns:k*ns);
            d_i = di_iter(:,ri:ns:nj*ns);
            mu_iter(:,ri) = ([log(p_d),x2e_iter(:,2:end)].*v_i*theta2w(:,1))...
                +[log(p_d),x2e_iter(:,2:end)].*(d_i*theta2w(:,2:nj+1)')*ones(k,1);
        end
        eg = exp(mu_iter + kron(ones(1,ns),mval));
        
        s_i = zeros(size(eg));
        hy_omega = zeros(max(unique(mid_iter)));
        hy_s = zeros(max(unique(mid_iter)),1);
        
        for city_t = min(cd_range_iter):max(cd_range_iter)
            trows = find(cdid_hy_iter == city_t);
            time_trows = trows(hy_temp_iter(trows) == t);
            
            s_i(time_trows,:) = eg(time_trows,:)./(ones(length(time_trows),1)*(1 + sum(eg(time_trows,:))));
            
            tempn = size(time_trows,1);
            temp0 = (s_i(time_trows,:).*alpha_i_iter(time_trows,:))*s_i(time_trows,:)'/ns;
            temp1 = sum((s_i(time_trows,:).*alpha_i_iter(time_trows,:))')'/ns;
            H = (diag(temp1) - temp0)./(p_d(time_trows,:)*ones(1,tempn)).*tts(time_trows);
            
            omega = H.*omega_sta_iter{city_t};
            
            modelid = mid_iter(time_trows);
            city_mktsize = unique(mktsizets_iter(time_trows));
            if length(city_mktsize) ~= 1, error('Market size not constant within city'); end
            s_j_city = mean(s_i(time_trows,:), 2);
            
            for row = 1:length(time_trows)
                for col = 1:length(time_trows)
                    hy_omega(modelid(row),modelid(col)) = hy_omega(modelid(row),modelid(col)) + omega(row,col)*city_mktsize;
                end
                hy_s(modelid(row),1) = hy_s(modelid(row),1) + s_j_city(row)*city_mktsize;
            end
        end
        
        nanproduct = find(hy_s==0);
        hy_s(nanproduct)=[];
        hy_omega(nanproduct,:)=[];
        hy_omega(:,nanproduct)=[];
        
        mkup_model = -inv(hy_omega)*hy_s;
        
        temp_models = [1:max(unique(mid_iter))]';
        temp_models(nanproduct)=[];
        
        mkup_e = zeros(length(p_e),1);
        for index = 1:length(temp_models)
            mkup_e(find(mid_iter==temp_models(index))) = mkup_model(index);
        end
        
        pnew = mkup_e + cost_sim_iter;
        norm = max(abs(pnew - p_e));
        avgnorm = mean(abs(pnew - p_e));
        p_e = pnew;
        i_iter = i_iter + 1;
    end
    
    pr_sim = p_e; % No subsidy equilibrium price
    
    % Store results for this half-year
    p_s_0_all{t} = p_s_0;
    pr_sim_all{t} = pr_sim;
    subsidy_l_all{t} = subsidy_l;
    sub_tt0_all{t} = sub_tt0;
    tts_all{t} = tts;
    qcost_all{t} = qcost;
    sc_tt_all{t} = sc_tt;
end

pass_through_demand_hy = cell(8,1);
pass_through_rate_demand_hy = cell(8,1);
for hy_iter = 1:8
    f_i = fnmb(data.hy == hy_iter);
    price_supply_hy=[p_s_0_all{hy_iter}, pr_sim_all{hy_iter}];
    subsidy_summarize_hy = [sub_tt0_all{hy_iter}, subsidy_l_all{hy_iter}];
    price_demand_hy = (price_supply_hy-subsidy_summarize_hy).*tts_all{hy_iter} + qcost_all{hy_iter};
    pass_through_demand_hy{hy_iter} =(price_demand_hy < price_supply_hy).* abs(  price_demand_hy - price_demand_hy(:,2))./subsidy_summarize_hy;
    pass_through_demand_hy{hy_iter} = pass_through_demand_hy{hy_iter}(:,1);
    pass_through_demand_hy{hy_iter}(isnan(pass_through_demand_hy{hy_iter})) = 0;  
end



%% Compare range and pass-through patterns across models
% Purpose:
%   Assemble EV model histories across half-years, identify comparable model
%   variants, and visualize the joint evolution of range and pass-through.
% Requires:
%   Period-specific pass-through arrays and vehicle model identifiers from the
%   full data panel.

evData = [];
for hy_iter = 1:8
    etc = (data.hy == hy_iter); % change the experiment city to beijing
    testmarket_temp = data(etc == 1,:);
    ptr_t = pass_through_demand_hy{hy_iter}; % we need to modify here
    testmarket_temp.ptr = ptr_t(:,1);
    %% Step 1: Filter for EV data only (ev == 1)
    evData_temp = testmarket_temp(testmarket_temp.ev == 1, :);  % Only analyze EV models
    
    evData = [evData; evData_temp];
end

%% Step 2: Find groups with same firm_name_unique + sub_brand + nameplate
% Create unique identifier for grouping
evData.combined_id = strcat(evData.firm_name_unique, '_', ...
                            evData.sub_brand, '_', ...
                            evData.nameplate);

% Find unique combinations
[unique_ids, ~, idx] = unique(evData.combined_id);

%% step 1: sort data
% refine data
evData_short = evData(:, {'combined_id', 'model_variant', 'cdid', 'hy', ...
                          'city_d', 'p_id', 'range', 'ev', 'ptr'});

% find models with multi variants
unique_ids = unique(evData_short.combined_id);
results = {};

for i = 1:length(unique_ids)
    id = unique_ids{i};
    all_data = evData_short(strcmp(evData_short.combined_id, id), :);
    
    variants = unique(all_data.model_variant);

    if length(variants) > 1
        % create series for each time
        variant_trends = struct();

        for v = 1:length(variants)
            variant_name = variants(v);
            v_data = all_data(all_data.model_variant==variant_name, :);
            
            % sort by hy
            v_data = sortrows(v_data, 'hy');
            
            variant_trends(v).name = variant_name;
            variant_trends(v).hy = v_data.hy;
            variant_trends(v).range = v_data.range;
            variant_trends(v).ptr = v_data.ptr;
        end
        
        results{end+1} = struct(...
            'combined_id', id, ...
            'variants', variants, ...
            'trends', variant_trends ...
        );
        
    end
end


%% filter data with ptr >0
evData_filtered = evData_short(evData_short.ptr ~= 0, :);

%% Paper Figure 7: Range and Pass-through Rate over Time
% Replication output:
%   range_ptr_hy.eps (Figure 7 in the paper).
% Figure content:
%   Two-panel plot of EV range and estimated subsidy pass-through rates over
%   time, with vehicle models distinguished by color.

figure('Position', [100, 100, 1200, 500]);

% filter all models and set colors
unique_ids = unique(evData_filtered.combined_id);
colors = lines(length(unique_ids));

time_labels = {'2016H1', '2016H2', '2017H1', '2017H2', ...
               '2018H1', '2018H2', '2019H1', '2019H2'};
% range subplot
subplot(1,2,1);
hold on;
for i = 1:length(unique_ids)
    id = unique_ids{i};
    id_mask = strcmp(evData_filtered.combined_id, id);
    scatter(evData_filtered.hy(id_mask), evData_filtered.range(id_mask), ...
            20, colors(i,:), 'filled');
end
xlabel('Time');
ylabel('Range');
title(sprintf('Range vs Time'));
% set x label
set(gca, 'XTick', 1:8, 'XTickLabel', time_labels);
xtickangle(45);  % 45 degree 
grid on;

% ptr subplot
subplot(1,2,2);
hold on;
for i = 1:length(unique_ids)
    id = unique_ids{i};
    id_mask = strcmp(evData_filtered.combined_id, id);
    scatter(evData_filtered.hy(id_mask), evData_filtered.ptr(id_mask), ...
            20, colors(i,:), 'filled');
end
xlabel('Time');
ylabel('Passthrough Rate');
title(sprintf('Passthrough Rate vs Time'));
% set x label
set(gca, 'XTick', 1:8, 'XTickLabel', time_labels);
xtickangle(45);  % 45 degree
grid on;

filename_ptr_range = sprintf('range_ptr_hy.eps');
saveas(gcf,filename_ptr_range)
close(gcf)

%% Report paper tables
% Purpose:
%   Collect the formatted replication outputs for Paper Tables 5-9 and A5 in one
%   place after all estimation, elasticity, and counterfactual calculations.
%   This section only reorganizes and reports results computed above.

diary on
%% Paper Table 5: Estimation results for the demand side
% This reporting block only reformats estimates already computed above.
rep_t5_names = {'log(Price)', 'Constant', 'log(Horsepower)', ...
    'log(Operation costs)', 'log(Weight)', 'log(Size)', ...
    'EV', 'Import', 'AT', 'SUV'};
rep_t5_coef = nan(10,5);
rep_t5_se = nan(10,5);
rep_t5_coef(:,1) = ols_results.beta(1:10);
rep_t5_se(:,1) = abs(ols_results.beta(1:10)./ols_results.tstat(1:10));
rep_t5_coef(:,2) = tsls_results.beta(1:10);
rep_t5_se(:,2) = abs(tsls_results.beta(1:10)./tsls_results.tstat(1:10));
rep_t5_coef(:,3) = theta1(1:10);
rep_t5_se(:,3) = se(1:10);
rep_t5_coef(1:6,4) = theta2w(1:6,1);
rep_t5_se(1:6,4) = se2w(1:6,1);
rep_t5_coef(1,5) = theta2w(1,2);
rep_t5_se(1,5) = se2w(1,2);
rep_t5_tstat = rep_t5_coef./rep_t5_se;

fprintf('\n=================================================================================================================\n');
fprintf('TABLE 5: Estimation Results for the Demand Side\n');
fprintf('=================================================================================================================\n');
fprintf('%-24s %14s %14s %14s %14s %14s\n', ...
    'Variable', 'OLS', 'TSLS', 'GMM Mean', 'GMM Random', 'GMM Income');
fprintf('-----------------------------------------------------------------------------------------------------------------\n');
for rep_t5_i = 1:10
    rep_t5_coef_text = repmat({''},1,5);
    rep_t5_se_text = repmat({''},1,5);
    for rep_t5_j = 1:5
        if ~isnan(rep_t5_coef(rep_t5_i,rep_t5_j))
            rep_t5_stars = '';
            rep_t5_abs_t = abs(rep_t5_tstat(rep_t5_i,rep_t5_j));
            if rep_t5_abs_t >= 2.5758
                rep_t5_stars = '***';
            elseif rep_t5_abs_t >= 1.9600
                rep_t5_stars = '**';
            elseif rep_t5_abs_t >= 1.6449
                rep_t5_stars = '*';
            end
            rep_t5_coef_text{rep_t5_j} = sprintf('%.4f%s', ...
                rep_t5_coef(rep_t5_i,rep_t5_j), rep_t5_stars);
            rep_t5_se_text{rep_t5_j} = sprintf('(%.4f)', ...
                rep_t5_se(rep_t5_i,rep_t5_j));
        end
    end
    fprintf('%-24s %14s %14s %14s %14s %14s\n', ...
        rep_t5_names{rep_t5_i}, rep_t5_coef_text{:});
    fprintf('%-24s %14s %14s %14s %14s %14s\n', ...
        '', rep_t5_se_text{:});
end
fprintf('-----------------------------------------------------------------------------------------------------------------\n');
fprintf('%-24s %14d\n', 'Observations', nobs);
fprintf('%-24s %14s\n', 'City fixed effects', 'Yes');
fprintf('%-24s %14s\n', 'Brand fixed effects', 'Yes');
fprintf('%-24s %14s\n', 'Time fixed effects', 'Yes');
fprintf('=================================================================================================================\n\n');
clear rep_t5_names rep_t5_coef rep_t5_se rep_t5_tstat rep_t5_i rep_t5_j
clear rep_t5_coef_text rep_t5_se_text rep_t5_stars rep_t5_abs_t

%% Paper Table 6: Cross-price elasticity between product categories
% This reporting block reorganizes the category elasticities computed above.
rep_t6_values = [epsilon_avg(1,:); epsilon_sum(1,:)];
fprintf('\n=========================================================================================\n');
fprintf('TABLE 6: Cross-Price Elasticity between Product Categories\n');
fprintf('=========================================================================================\n');
fprintf('%-8s %-14s %12s %12s %12s %12s\n', ...
    '', '', 'Domestic EV', 'Domestic ICEV', 'Imported EV', 'Imported ICEV');
fprintf('-----------------------------------------------------------------------------------------\n');
fprintf('%-8s %-14s %12.5f %12.6f %12.6f %12.6f\n', ...
    'EV', 'Individual', rep_t6_values(1,:));
fprintf('%-8s %-14s %12.4f %12.4f %12.6f %12.5f\n', ...
    '', 'Aggregate', rep_t6_values(2,:));
fprintf('=========================================================================================\n\n');
clear rep_t6_values

%% Paper Table 7: Subsidy pass-through to consumers
% This reporting block converts the pass-through ratios above to percentages.
rep_t7_values = 100*pass_through_rate_demand;
rep_t7_rows = {'EV manufacturers', 'Hybrid manufacturers', 'All manufacturers'};
fprintf('\n================================================================================\n');
fprintf('TABLE 7: Subsidy Pass-through to Consumers\n');
fprintf('================================================================================\n');
fprintf('%-26s %12s %12s %12s %12s\n', ...
    'Scenarios', 'Observed', '(1)', '(3)', '(4)');
fprintf('--------------------------------------------------------------------------------\n');
for rep_t7_i = 1:3
    fprintf('%-26s %11.2f%% %11.2f%% %11.2f%% %11.2f%%\n', ...
        rep_t7_rows{rep_t7_i}, rep_t7_values(rep_t7_i,:));
end
fprintf('================================================================================\n\n');
clear rep_t7_values rep_t7_rows rep_t7_i

%% Paper Table 8: Pass-through by EV price quartile
% This reporting block converts the quartile pass-through ratios to percentages.
rep_t8_values = 100*pt_quantile;
fprintf('\n=======================================================================================================\n');
fprintf('TABLE 8: Subsidy Pass-through to Consumers by EV Price Quartiles\n');
fprintf('=======================================================================================================\n');
fprintf('%-30s %16s %16s %16s %16s\n', ...
    'EV price quartiles', '1st quartile', '2nd quartile', '3rd quartile', '4th quartile');
fprintf('-------------------------------------------------------------------------------------------------------\n');
fprintf('%-30s %15.2f%% %15.2f%% %15.2f%% %15.2f%%\n', ...
    'Pass-through to consumers', rep_t8_values);
fprintf('=======================================================================================================\n\n');
clear rep_t8_values

%% Paper Table 9: Cost-benefit analysis of EV subsidies
% This reporting block reorganizes welfare components computed above.
% Monetary values are reported in RMB billions.
rep_t9_cv = [0; cvs(:)]'/1e5;
rep_t9_profit = [pi_category(1,1:2:11)', pi_category(:,2:2:12)']/1e5;
rep_t9_subsidy = subsidy_c(:)'/1e5;
rep_t9_sales_subtotal = rep_t9_cv + sum(rep_t9_profit,1) - rep_t9_subsidy;
rep_t9_ev_coal = ex_ev(:)'/1e9;
rep_t9_ev_gas = ex_evd(:)'/1e9;
rep_t9_icev = ex_fv(:)'/1e9;
rep_t9_ext_coal = rep_t9_ev_coal + rep_t9_icev;
rep_t9_ext_gas = rep_t9_ev_gas + rep_t9_icev;
rep_t9_total_coal = rep_t9_sales_subtotal - rep_t9_ext_coal;
rep_t9_total_gas = rep_t9_sales_subtotal - rep_t9_ext_gas;

fprintf('\n====================================================================================================================================================\n');
fprintf('TABLE 9: Cost-Benefit Analysis of EV Subsidies\n');
fprintf('All monetary values are in RMB billions.\n');
fprintf('====================================================================================================================================================\n');
fprintf('%-34s %-34s %12s %12s %12s %12s %12s\n', ...
    'Component', 'Category', 'Observed', '(1)', '(2)', '(3)', '(4)');
fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    'Compensating variation', '', rep_t9_cv);
fprintf('%-34s\n', 'Profits');
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    'Domestic', 'EV manufacturers', rep_t9_profit(1,:));
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '', 'ICEV manufacturers', rep_t9_profit(2,:));
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '', 'Hybrid manufacturers', rep_t9_profit(3,:));
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    'Imported', 'EV manufacturers', rep_t9_profit(4,:));
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '', 'ICEV manufacturers', rep_t9_profit(5,:));
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    '', 'Hybrid manufacturers', rep_t9_profit(6,:));
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    'Subsidy', '', rep_t9_subsidy);
fprintf('%-34s %-34s %12.4f %12.4f %12.4f %12.4f %12.4f\n', ...
    'Subtotal for sales', '', rep_t9_sales_subtotal);
fprintf('%-34s\n', 'Externalities');
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'EVs', 'Coal-fired electricity', rep_t9_ev_coal);
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'EVs', 'Natural-gas-powered electricity', rep_t9_ev_gas);
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'ICEVs', '', rep_t9_icev);
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'Subtotal for externalities', 'Coal-fired electricity', rep_t9_ext_coal);
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'Subtotal for externalities', 'Natural-gas-powered electricity', rep_t9_ext_gas);
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'Total', 'Coal-fired electricity', rep_t9_total_coal);
fprintf('%-34s %-34s %12.1f %12.1f %12.1f %12.1f %12.1f\n', ...
    'Total', 'Natural-gas-powered electricity', rep_t9_total_gas);
fprintf('====================================================================================================================================================\n\n');
clear rep_t9_cv rep_t9_profit rep_t9_subsidy rep_t9_sales_subtotal
clear rep_t9_ev_coal rep_t9_ev_gas rep_t9_icev rep_t9_ext_coal rep_t9_ext_gas
clear rep_t9_total_coal rep_t9_total_gas

%% Paper Table A5: Expected subsidy by income quartile
% This reporting block formats the unconditional expected subsidies computed above.
% Values include zero subsidy for simulated non-purchasers and are reported in RMB.
rep_ta5_values = expected_sub_quantile;
fprintf('\n=======================================================================================================\n');
fprintf('TABLE A5: Expected Subsidy to Consumers by Income Quartiles\n');
fprintf('=======================================================================================================\n');
fprintf('%-38s %14s %14s %14s %14s\n', ...
    'Income quantiles', '1st quartile', '2nd quartile', '3rd quartile', '4th quartile');
fprintf('-------------------------------------------------------------------------------------------------------\n');
fprintf('%-38s %14.2f %14.2f %14.2f %14.2f\n', ...
    'Expected subsidy to consumers (RMB)', rep_ta5_values);
fprintf('=======================================================================================================\n\n');
clear rep_ta5_values
diary off
