function [f1, f3] = meanval_micro_PP_norm(theta2)
% MEANVAL_MICRO_PP_NORM Recover mean utilities by BLP contraction mapping using parallel processing.
%
% INPUT
%   theta2 - free nonlinear parameters.
%
% OUTPUTS
%   f1 / delta  - [N x 1] recovered mean utilities, equal to log(mval).
%   f3 / expmu  - [N x ns] exponentiated heterogeneous utility components.
%
% REQUIRED GLOBAL STATE
%   theti, thetj - locations of the free nonlinear parameters.
%   x2           - [N x K] random-coefficient product characteristics.
%   s_jt         - [N x 1] observed product shares.
%   ns           - number of simulation draws.
%   cdindex      - [M x 1] cumulative last-row index for each market.
%   nmkt, nprod  - number of markets and product counts by market.
%
% FILE INPUTS / OUTPUTS
%   Loads mvalold from mvalold_norm.mat as the starting value. After a valid
%   solution, overwrites that file with updated mvalold and oldt2.
%
% SOURCE
%   Modified from Aviv Nevo's random-coefficients demand code.

global theti thetj x2 s_jt ns niter cdindex nmkt nprod 

load mvalold_norm

index_begin=[0;cdindex(1:end-1)]+1;

flag=-1;
tol=1e-15;
maxiter=500;

theta2w = full(sparse(theti,thetj,theta2));

expmu = exp(mufunc_norm(x2,theta2w));

mval = zeros(length(x2), 1);

mval_mat=zeros(max(nprod), nmkt); %make the mvalold into a matrix. otherwise, parfor won't work
nprod_max=max(nprod);

parfor m=1:nmkt
    norm = 1;
    avgnorm = 1;
    i = 0;
    
    aa=index_begin(m);
    bb=cdindex(m);
    
    mvalold_m=mvalold(aa:bb,:);
    expmu_m=expmu(aa:bb,:);
    s_jt_m=s_jt(aa:bb,:);  
    
       while i<10
        indshr_m = ind_sh_m(mvalold_m,expmu_m, nprod(m), ns);
         mktshr_m = sum(indshr_m')'/ns;
           mval_m = mvalold_m.*(s_jt_m)./(mktshr_m); 
         

   		mvalold_m = mval_m;
         i = i + 1;
      end
    while i<maxiter & norm > tol & avgnorm > 1e-3*tol

        indshr_m = ind_sh_m(mvalold_m,expmu_m, nprod(m), ns);
        mktshr_m = sum(indshr_m')'/ns;
        
        log_mval =log(mvalold_m) + log(s_jt_m) - log(mktshr_m);

        mval_m=exp(log_mval);
  
        t = abs(mval_m-mvalold_m);
	  	norm = max(t);
        avgnorm = mean(t);
  		mvalold_m = mval_m;
        i = i + 1;
    end
    if i>=maxiter
        display(['contraction mapping over limit in: ' num2str(m)]);
        num2str([norm*1000000 avgnorm*1000000])
    end
    mval_m_temp=[mval_m; zeros(nprod_max-nprod(m),1)];
    mval_mat(:,m)=mval_m_temp;   
    
end

for m=1:nmkt
    mval(index_begin(m):cdindex(m),:)=mval_mat(1:nprod(m),m);
end

if max(isnan(mval)) < 1;
   mvalold = mval;
   oldt2 = theta2;
   save mvalold_norm mvalold oldt2;
else
   display(['mvalold has problem']);
end   
f1 = log(mval); f3 = expmu;
