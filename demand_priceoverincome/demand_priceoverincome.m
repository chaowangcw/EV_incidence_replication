% Replication driver for the price-over-income random-coefficients
% vehicle-demand specification reported in Table A4.
%
% REPLICATION WORKFLOW
%   - Initialize paths, software dependencies, random seeds, and simulation
%     settings.
%   - Load the vehicle-market data and trim the bottom one percent of MSRP.
%   - Express price and continuous product characteristics in levels.
%   - Generate random-coefficient and income simulation draws.
%   - Construct differentiation instruments and estimate OLS and 2SLS models.
%   - Estimate demand by two-step GMM with price/income heterogeneity and
%     income micro moments.
%   - Recover standard errors and elasticity diagnostics.
%
% EXECUTION NOTES
%   Run this script from its demand_priceoverincome subdirectory. Data, common
%   functions, the econometric toolbox, and diary output are located in the
%   replication root directory. Specification-specific functions remain in
%   this subdirectory and must take precedence over same-named root functions.
%   The program uses global variables and intermediate MAT files, so sections
%   should be run in order. Random seeds are set below for reproducibility.
%   This file should be opened in UTF-8 mode.
%
% PRINCIPAL OUTPUTS
%   Price-over-income demand estimates underlying Table A4, covariance
%   matrices, diary logs, and MAT checkpoints.

%% Initialize replication environment
clear all
clc
global invA ns x1 x2 s_jt IV IVs vfull dfull theta1 theti thetj cdid cdindex nobs nmkt sale_tax
global thet2 op xs  nd  NewW  niter nprod  tax fnmb nseudo seudos alpha_i nobs_c
global range_tt ev_tt  import_tt  vi di tts qcost x1e x2e ehat mktsizets
global  omega_sta cost_sim p_s_o local_tt sc_tt k NewincW inc_r s_inc
global brand_matrix city_matrix inc_c_W inc_b_W subdec subextra_id eg0  inc_IV
global GLOBAL_SCALE_FACTOR PROJECT_ROOT_DIR nobs_c

SPECIFICATION_DIR = fileparts(mfilename('fullpath'));
PROJECT_ROOT_DIR = fileparts(SPECIFICATION_DIR);
TOOLBOX_DIR = fullfile(PROJECT_ROOT_DIR, 'econometric toolbox');
% Keep specification-specific functions in the current directory and add
% common replication functions and the bundled econometric toolbox.
cd(SPECIFICATION_DIR)
addpath(PROJECT_ROOT_DIR, '-begin');
addpath(genpath(TOOLBOX_DIR));
RUN_TAG = datestr(now, 'yyyymmdd_HHMMSS');
format short

GLOBAL_SCALE_FACTOR = 12/100000; % scaling factor for family annual income (hundred thousand)
PRICE_SCALE_FACTOR = 10; % scaling factor for price (hundred thousand thousand)
if isempty(gcp('nocreate')) == 1
    parpool(24) % change according to the number of cores in your machine
end

rand('state', 1234);%4321);
randn('state', 1234);% 1234);  %fixed random draws to reproduce the estimation results
flag_gen = 1; % 0: load pre-generated random draws; 1: generate draws
flag_mail = 1; %1: set email notice
ns = 500;       % number of simulated "indviduals" per market %
DATA_FILE = fullfile(PROJECT_ROOT_DIR, 'data_incidence.csv');
data = dataset('file',DATA_FILE, 'Delimiter','comma','ReadVarNames',true);

msrpcutoff1 = prctile(data.msrp,1);
data = data(data.msrp>msrpcutoff1,:);
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
price = data.msrp/PRICE_SCALE_FACTOR;
subsidy_local_net = data.subsidy_local_net/PRICE_SCALE_FACTOR;
s_c = data.s_c/PRICE_SCALE_FACTOR;
quota_cost = data.quota_cost/10^4/PRICE_SCALE_FACTOR;
evrange = data.range;
op = price - subsidy_local_net;
op(find(ev == 0 ),1) = op(find(ev == 0 ),1)...
    *(1+1/1.13*.1) + quota_cost(find(ev == 0 ),1);
sale_tax = 1+ (ev == 0 )*(1/1.13*.1) ;
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
odds = data.odds;
quota_d = data.quota_d;
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

power = exp(lnpower)/100;
uni_cost = exp(lnuni_cost)/100;
weight = exp(lnweight);
x1 = [ones(nobs,1) power uni_cost weight exp(lnsize)/10 ev import AT...
    suv city_dum brand_dum hy_dum];
x2 = [op ones(nobs,1) power uni_cost weight exp(lnsize)/10];

%% random coefficient on price
% starting values. zero elements in the following matrix correspond to %
% coeff that will not be max over,i.e are fixed at zero. %
theta2w=[0  -0.1083
    0.3526   0
    0.0563   0
    -0.1719, 0
    0.3957, 0
    0.3203, 0];

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

%% generate vfull dfull

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
        demogr(i,:) = exp(winsor(tempdfull',[1,99])');
        clear tempdfull
    end
    dfull = log(demogr(cdid,:)/100*2.7); % see as family income
    vfull = v(cdid, :);

    filename=strcat('ps2_', num2str(ns));
    save (filename, 'x1', 'x2', 'demogr', 'v', 's_jt', 'op');
    clear halt halt_inc draws draws_inc temp
else
    filename=strcat('ps2_', num2str(ns));
    load(filename, 'demogr', 'v')
    dfull = log(demogr(cdid,:)/100*2.7); % see as family income
    vfull = v(cdid, :);
end



%% Construct differentiation instruments and weighting inputs
x_diffiv = [weight,power,uni_cost,exp(lnsize)/10];
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
     diffIV_quadratic(:,[1,4,7]) diffIV_local(:,[1,2,21]) ...
];
IV = [x1, iv];
inc_IV = IV(inc_r,:);
% Drop zero and singleton instrument columns in the micro-moment sample.
inc_IV = inc_IV(:, sum(inc_IV ~= 0, 1) >= 2);
nobs_c = size(inc_IV,1);
invA = inv([IV'*IV]);
% Logit results and save the mean utility as initial values for the search below

%% Estimate baseline Logit and linear IV specifications
% compute the outside good market share by market
for i = 1:nmkt
    sum1(i,1) = sum(s_jt(find(cdid == i)));
end
outshr = 1.0 - sum1(cdid,:);

y = log(s_jt) - log(outshr);
mid = [x1, op]'*IV*invA*IV';
bt = inv(mid*[x1, op])*mid*y;
mvalold = [x1, op]*bt;
oldt2 = zeros(size(theta2_0));
mvalold = exp(mvalold);

ols_results=ols(y, [op, x1]);
tsls_results=tsls(y, [op], x1, IV);
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
    1e-9,'TolFun',1e-6,'MaxFunctionEvaluations',1e9);

starts = 5;
theta20 = theta2_0;
startvalues = repmat(theta20',starts,1) + cat(2,zeros(size(theta20)),rand([starts-1,size(theta20',2)] )' *1 )';

GMPEC = 1.0e20;
CPUtMPEC = 0;
FuncEvalMPEC = 0;
LB = [zeros(size(theta20,1)-1, 1); -inf];
for reps=1:starts,
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

diary_path = fullfile(PROJECT_ROOT_DIR, 'diary');
if ~exist(diary_path, 'dir'), mkdir(diary_path); end 
disp(['diary path is set as: ' diary_path]);
diary_file = fullfile(diary_path, ...
    ['demand_priceoverincome_' RUN_TAG '.txt']);
diary(diary_file);
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
save(fullfile(PROJECT_ROOT_DIR, ...
    ['demand_priceoverincome_' RUN_TAG '.mat']), '-v7.3')

%% Paper Table A4: Demand estimation with price-over-income heterogeneity
rep_poi_names = {'Price', 'Constant', 'Power', 'Fuel Cost', 'Weight', ...
    'Size', 'EV', 'Import', 'AT', 'SUV'};
rep_poi_coef = nan(10,3);
rep_poi_se = nan(10,3);
rep_poi_coef(:,1) = ols_results.beta(1:10);
rep_poi_se(:,1) = abs(ols_results.beta(1:10)./ols_results.tstat(1:10));
rep_poi_coef(:,2) = tsls_results.beta(1:10);
rep_poi_se(:,2) = abs(tsls_results.beta(1:10)./tsls_results.tstat(1:10));
rep_poi_coef(:,3) = theta1(1:10);
rep_poi_se(:,3) = se(1:10);
rep_poi_tstat = rep_poi_coef./rep_poi_se;

rep_poi_random_names = {'Constant', 'Power', 'Fuel Cost', 'Weight', 'Size'};
rep_poi_random_coef = theta2w(2:6,1);
rep_poi_random_se = se2w(2:6,1);
rep_poi_random_tstat = rep_poi_random_coef./rep_poi_random_se;
rep_poi_income_coef = theta2w(1,2);
rep_poi_income_se = se2w(1,2);
rep_poi_income_tstat = rep_poi_income_coef/rep_poi_income_se;

fprintf('\n=============================================================================================\n');
fprintf('TABLE A4: Demand Estimation Results: BLP Random Coefficients Models\n');
fprintf('=============================================================================================\n');
fprintf('%-28s %34s %24s\n', ...
    'Variables', 'Homogeneous Preference', 'Heterogeneous Preference');
fprintf('%-28s %16s %16s %16s\n', '', 'OLS', '2SLS', 'BLP');
fprintf('%-28s %16s %16s %16s\n', '', '(1)', '(2)', '(3)');
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('%-28s %50s\n', '', 'Panel A: Mean Utility Parameters');
for rep_poi_i = 1:10
    rep_poi_coef_text = repmat({''},1,3);
    rep_poi_se_text = repmat({''},1,3);
    for rep_poi_j = 1:3
        rep_poi_stars = '';
        rep_poi_abs_t = abs(rep_poi_tstat(rep_poi_i,rep_poi_j));
        if rep_poi_abs_t >= 2.5758
            rep_poi_stars = '***';
        elseif rep_poi_abs_t >= 1.9600
            rep_poi_stars = '**';
        elseif rep_poi_abs_t >= 1.6449
            rep_poi_stars = '*';
        end
        rep_poi_coef_text{rep_poi_j} = sprintf('%.3f%s', ...
            rep_poi_coef(rep_poi_i,rep_poi_j), rep_poi_stars);
        rep_poi_se_text{rep_poi_j} = sprintf('(%.3f)', ...
            rep_poi_se(rep_poi_i,rep_poi_j));
    end
    fprintf('%-28s %16s %16s %16s\n', ...
        rep_poi_names{rep_poi_i}, rep_poi_coef_text{:});
    fprintf('%-28s %16s %16s %16s\n', '', rep_poi_se_text{:});
end
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('%-28s %50s\n', '', 'Panel B: BLP Random Coefficients');
fprintf('%s\n', 'Variance parameters (sigmas)');
for rep_poi_i = 1:5
    rep_poi_stars = '';
    rep_poi_abs_t = abs(rep_poi_random_tstat(rep_poi_i));
    if rep_poi_abs_t >= 2.5758
        rep_poi_stars = '***';
    elseif rep_poi_abs_t >= 1.9600
        rep_poi_stars = '**';
    elseif rep_poi_abs_t >= 1.6449
        rep_poi_stars = '*';
    end
    fprintf('%-28s %16s %16s %16s\n', rep_poi_random_names{rep_poi_i}, ...
        '', '', sprintf('%.3f%s', rep_poi_random_coef(rep_poi_i), rep_poi_stars));
    fprintf('%-28s %16s %16s %16s\n', '', '', '', ...
        sprintf('(%.3f)', rep_poi_random_se(rep_poi_i)));
end
fprintf('%s\n', 'Demographic interaction');
rep_poi_stars = '';
rep_poi_abs_t = abs(rep_poi_income_tstat);
if rep_poi_abs_t >= 2.5758
    rep_poi_stars = '***';
elseif rep_poi_abs_t >= 1.9600
    rep_poi_stars = '**';
elseif rep_poi_abs_t >= 1.6449
    rep_poi_stars = '*';
end
fprintf('%-28s %16s %16s %16s\n', 'Price/Income', '', '', ...
    sprintf('%.3f%s', rep_poi_income_coef, rep_poi_stars));
fprintf('%-28s %16s %16s %16s\n', '', '', '', ...
    sprintf('(%.3f)', rep_poi_income_se));
fprintf('---------------------------------------------------------------------------------------------\n');
fprintf('%-28s %16s %16s %16s\n', 'Instrumental Variables', '', ...
    'Differentiation IVs', 'Differentiation IVs');
fprintf('=============================================================================================\n');
fprintf('Observations: %d\n', nobs);
fprintf('Final GMM objective: %.6f\n', fval);
fprintf('Final optimizer exit flag: %d\n', flag2);
fprintf('Running time: %.2f minutes\n\n', time/60);

diary off

clear rep_poi_names rep_poi_coef rep_poi_se rep_poi_tstat
clear rep_poi_random_names rep_poi_random_coef rep_poi_random_se
clear rep_poi_random_tstat rep_poi_income_coef rep_poi_income_se
clear rep_poi_income_tstat rep_poi_i rep_poi_j
clear rep_poi_coef_text rep_poi_se_text rep_poi_stars rep_poi_abs_t
