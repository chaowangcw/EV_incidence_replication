function [Z, col_names] = build_differentiation_iv(X, market_ids, firm_ids, version, interact)
% BUILD_DIFFERENTIATION_IV  Excluded differentiation instruments a la Gandhi & Houde (2017).
%
%   [Z, names] = build_differentiation_iv(X, market_ids, firm_ids, version, interact)
%
% Inputs
%   X           : [N x L] matrix of product characteristics (columns = characteristics l=1..L)
%   market_ids  : [N x 1] vector (string/cell/categorical/numeric) of market identifiers
%   firm_ids    : [N x 1] vector (string/cell/categorical/numeric) of firm identifiers
%   version     : 'local' (default) or 'quadratic'
%                 - local:   indicators 1(|d_{jktl}| < SD_l) as weights
%                 - quadratic: continuous d_{jktl} * d_{jktl'} (or d^2 when l'=l)
%   interact    : logical (default = false)
%                 - local:    indicator(|d_l|<SD_l) * d_{l'} for all l'
%                 - quadratic: d_{l} * d_{l'} for all l' >= l (unique pairs)
%
% Outputs
%   Z           : [N x K] instrument matrix (columns ordered as documented below)
%   col_names   : {1 x K} cell array of column labels
%
% Column layout:
%   Without interactions (L characteristics):
%     [ Other: l=1..L , Rival: l=1..L ]  ->  2*L columns
%
%   With interactions:
%     version='local'    : for each base l, combine with all l' (1..L):
%                          [ Other: (l,l'=1..L), Rival: (l,l'=1..L) ] -> 2*L*L columns
%     version='quadratic': use unique pairs (l,l'), l' >= l:
%                          [ Other: pairs, Rival: pairs ] -> 2*L*(L+1)/2 columns
%
% Notes:
%   - SD_l is the std of pairwise differences d_{jktl} pooled across all markets.
%   - Self-pairs are excluded; "Other" sums only over same-firm (k ~= j), "Rival" over different firms.
%   - Markets with only one product produce zeros for that market.
%
% Reference: Gandhi, A. & Houde, J.-F. (2017). "Measuring Substitution Patterns in Differentiated Products
%            Industries using Market-Level Data."

    if nargin < 4 || isempty(version),  version  = 'local';     end
    if nargin < 5 || isempty(interact), interact = false;        end

    validateattributes(X, {'numeric'}, {'2d','nonempty','finite'});
    N = size(X,1);
    L = size(X,2);

    % Normalize/encode ids
    market_ids = normalize_ids(market_ids);
    firm_ids   = normalize_ids(firm_ids);

    % Precompute SD_l over all pairwise diffs across markets (for the local version)
    SD = compute_pairwise_sd(X, market_ids);

    % Prepare column plan (which interaction indices to include)
    if ~interact
        % no interactions -> use only (l) itself
        pairs = [(1:L)', (1:L)'];           % l with itself
        is_unique_pairs = false;
    else
        switch lower(version)
            case 'local'
                % local + interactions: for each base l, include all l' (1..L)
                [L2, L1] = ndgrid(1:L, 1:L);
                pairs = [L1(:), L2(:)];     % (l, l') all combinations
                is_unique_pairs = false;
            case 'quadratic'
                % quadratic + interactions: include unique pairs (l, l'), l' >= l
                pairs = [];
                for ell = 1:L
                    pairs = [pairs; [ell*ones(L-ell+1,1), (ell:L)']]; %
                end
                is_unique_pairs = true;
            otherwise
                error('Unknown version: %s', version);
        end
    end
    P = size(pairs,1);  % number of interaction "channels" per block (Other or Rival)

    % Total columns: 2 blocks (Other, Rival)
    K = 2 * P;
    Z = zeros(N, K);
    col_names = cell(1, K);

    % Build per market to avoid spurious cross-market interactions
    markets = unique(market_ids);
    col = 0;

    % We'll fill blocks in the order:
    %   Block 1: Other (within-firm), all pairs in 'pairs' in order
    %   Block 2: Rival (different-firm), all pairs in 'pairs' in order
    %
    % To do this cleanly, we'll accumulate within temporary arrays and then place them.

    Other_block = zeros(N, P);
    Rival_block = zeros(N, P);

    for m = markets.'
        idx = (market_ids == m);
        J = sum(idx);  % number of products in this market
        if J <= 1
            continue;  % leave zeros for singletons
        end

        Xm = X(idx, :);
        firms_m = firm_ids(idx);

        % 3D difference tensor D: [J x J x L], D(j,k,l) = x_{ktl} - x_{jtl}
        % (k along columns, j along rows so we can sum across k)
        D = compute_diff_tensor(Xm);  % JxJxL
        % Ignore self-pairs (vectorized diagonal masking across all l)
        diag_idx = repmat((1:J+1:J*J).', 1, L) + repmat((0:L-1)*(J*J), J, 1);
        D(diag_idx) = NaN;          % linear indices into the 3-D array


        % Masks for "Other" (same firm, k ~= j) and "Rival" (different firm)
        same_firm = bsxfun(@eq, firms_m, firms_m.');  % JxJ logical
        same_firm(1:J+1:end) = false;                 % exclude diagonal
        rival_firm = ~same_firm & ~eye(J);

        % Compute market contributions for each pair (l,l')
        switch lower(version)
            case 'local'
                % Indicators: 1(|d_l| < SD_l)
                for p = 1:P
                    ell  = pairs(p,1);
                    ellp = pairs(p,2);

                    W = abs(D(:,:,ell)) < SD(ell);   % JxJ logical weights

                    if interact
                        % summand = 1(|d_ell|<SD_ell) * d_{ell'}
                        summand = W .* D(:,:,ellp);
                    else
                        % no interactions -> just the indicator itself
                        summand = W;
                    end

                    % Sum over k for each j, split by firm relation
                    Other_block(idx, p) = nansum(summand .* same_firm, 2);
                    Rival_block(idx, p) = nansum(summand .* rival_firm, 2);
                end

            case 'quadratic'
                % Quadratic distances:
                %   no-interact : d_{ell}.^2
                %   interact    : d_{ell} .* d_{ellp}  (ell'>=ell in pairs)
                for p = 1:P
                    ell  = pairs(p,1);
                    ellp = pairs(p,2);

                    if interact
                        summand = D(:,:,ell) .* D(:,:,ellp);
                    else
                        summand = D(:,:,ell).^2;
                    end

                    Other_block(idx, p) = nansum(summand .* same_firm, 2);
                    Rival_block(idx, p) = nansum(summand .* rival_firm, 2);
                end

            otherwise
                error('Unknown version: %s', version);
        end
    end

    % Assemble Z and names
    Z(:, 1:P)       = Other_block;
    Z(:, P+1:2*P)   = Rival_block;

    % Column names
    % Other block:
    for p = 1:P
        ell  = pairs(p,1);
        ellp = pairs(p,2);
        if ~interact
            base = sprintf('Other_%s', base_tag(version));
        else
            if strcmpi(version,'local')
                base = sprintf('Other_local_l%u_by_l%u', ellp, ell);
            else
                if is_unique_pairs
                    base = sprintf('Other_quad_l%u_x_l%u', ell, ellp);
                else
                    base = sprintf('Other_quad_l%u_by_l%u', ell, ellp);
                end
            end
        end
        col_names{p} = pretty_name(base, ell, ellp, version, interact);
    end
    % Rival block:
    for p = 1:P
        ell  = pairs(p,1);
        ellp = pairs(p,2);
        if ~interact
            base = sprintf('Rival_%s', base_tag(version));
        else
            if strcmpi(version,'local')
                base = sprintf('Rival_local_l%u_by_l%u', ellp, ell);
            else
                if is_unique_pairs
                    base = sprintf('Rival_quad_l%u_x_l%u', ell, ellp);
                else
                    base = sprintf('Rival_quad_l%u_by_l%u', ell, ellp);
                end
            end
        end
        col_names{P+p} = pretty_name(base, ell, ellp, version, interact);
    end
end

% === helpers ===

function ids = normalize_ids(x)
% Map arbitrary id vector to compact integer codes 1..K (stable) as a column.
% Robust to numeric, logical, categorical, string, cellstr, or char.

    % force column
    x = x(:);

    if iscategorical(x)
        % categorical -> stable integer codes via grp2idx
        [ids, ~] = grp2idx(x);

    elseif isstring(x)
        % string array -> grp2idx
        [ids, ~] = grp2idx(x);

    elseif iscellstr(x)
        % cellstr -> grp2idx
        [ids, ~] = grp2idx(x);

    elseif ischar(x)
        % char matrix (NxM) -> cellstr -> grp2idx
        [ids, ~] = grp2idx(cellstr(x));

    else
        % numeric / logical -> unique on the raw numbers (no text conversion)
        [~,~,ids] = unique(x, 'stable');
    end

    ids = ids(:);  % ensure column
end


function SD = compute_pairwise_sd(X, market_ids)
    % Compute SD_l over pooled pairwise diffs across markets.
    N = size(X,1); L = size(X,2);
    assert(numel(market_ids) == N, ...
        'market_ids must have exactly size(X,1) elements.');
    market_ids = market_ids(:);
    SD = zeros(1,L);

    markets = unique(market_ids,'stable');
    for ell = 1:L
        diffs = [];
        for mi = 1:numel(markets)
            m = markets(mi);
            idx = (market_ids == m);   % logical Nx1
            Xm  = X(idx, ell);         % Jx1
            J   = numel(Xm);
            if J <= 1, continue; end
            D = Xm.' - Xm;             % JxJ pairwise diffs
            D(1:J+1:end) = [];         % drop diagonal to vectorize
            diffs = [diffs; D(:)];     
        end
        if isempty(diffs)
            SD(ell) = eps;             % avoid zero threshold
        else
            SD(ell) = std(diffs,0,'omitnan');
            if SD(ell) == 0, SD(ell) = eps; end
        end
    end
end


function D = compute_diff_tensor(Xm)
    % For Xm = [J x L], return D = [J x J x L] with D(j,k,l) = Xm(k,l) - Xm(j,l)
    [J, L] = size(Xm);
    D = zeros(J, J, L);
    for ell = 1:L
        v = Xm(:, ell);
        % D(:,:,ell) = v.' - v;  % k across columns, j across rows
        D(:,:,ell) = v - v.';  % Now D(j,k,l) = Xm(j,l) - Xm(k,l)
    end
end

function s = base_tag(version)
    if strcmpi(version,'local')
        s = 'local';
    else
        s = 'quad';
    end
end

function nm = pretty_name(base, ell, ellp, version, interact)
    if ~interact
        % Other_local or Rival_quad, plus the characteristic index
        nm = sprintf('%s_l%u', base, ell);
    else
        nm = base; % already includes l and l'
    end
end
