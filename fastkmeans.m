function [label,optMean,err]=fastkmeans(X, k, options)
    par = [];
    par.repeat = 1;
    par.blockLen = 1;
    par.start = 1;

    if isfield(options,'repeat'), par.repeat = options.repeat; end
    if isfield(options,'blockLen'), par.blockLen = options.blockLen; end
    if isfield(options,'start'), par.start = options.start; end
     
    blockNum = ceil(par.repeat/par.blockLen);
    %print(['block number is',num2str(blockNum)]);
    candSeedPool = randperm(max(10000, blockNum));
    
    labelCell = cell(1,blockNum);
    optMeanCell = cell(1,blockNum);
    optErrCell = cell(1,blockNum);
    
    parfor i=1: blockNum
        [label_tmp,optMean_tmp,err_tmp] = myKmeansL2(X,k,par,candSeedPool(i));
        optErrCell{i} = err_tmp;
        labelCell{i} = label_tmp;
        optMeanCell{i} = optMean_tmp;
    end
    
    [~,I]=min(cell2mat(optErrCell));
        
    label = labelCell{I};
    optMean = optMeanCell{I};
    err = optErrCell{I};
end


function [label,optMean,optErr] = myKmeansL2(X,k,par,randSeed)

if par.start==1, rng(randSeed); end 
label = []; optMean = []; optErr = -1; optK = -1; n = size(X,2);

for cnt = 1: par.blockLen
    last = 0; m = [];iter=1;
    switch par.start
        case 1
            m = initSeedsRand(X, k);
            n_m=size(m,2);
            %disp(['n_m ',num2str(n_m)]);
            %disp(['m size',num2str(size(m,2))]);
            
        case 2
            m = initSeedsKmeansPP(X, k);
        case 3
            m = X(:,randsample(n,k));
    end
    [~,currLabel] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1);
    n_cu=1;
    while any(currLabel ~= last') && iter<=100
        [u,~,currLabel] = unique(currLabel);
        k_init = length(u);
        E = sparse(1:n,currLabel,1,n,k_init,n);  % transform currLabel into indicator matrix
        m = X*(E*spdiags(1./sum(E,1)',0,k_init,k_init));    % compute m of each cluster,
        % creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d
        last = currLabel;
        [~,currLabel] = max(bsxfun(@minus,m'*X,dot(m,m,1)'/2),[],1); % assign samples to the nearest centers
        iter = iter+1;
        n_cu=n_cu+1;
    end
    
    D = bsxfun(@plus,dot(m,m,1)',dot(X,X,1))-2*(m'*X); 
  
    err = sum(min(D));
    %disp(['optErr: ', num2str(optErr),'  currErr: ', num2str(err),' optK: ', num2str(optK),' k: ', num2str(k_init)]);
    if (cnt == 1) || (err < optErr),label = currLabel'; optErr = err; optMean = m; optK = k_init;end
    %m is the centroids
end
end




% Random initialization
function m = initSeedsRand(X, k)
[d,n] = size(X);
%disp(['k is ',num2str(k)]);
label = ceil(k*rand(1,n));
[u,~,currLabel] = unique(label);   % remove empty clusters
k_init = length(u);
%disp(['k_init ',num2str(k_init)]);
E = sparse(1:n,currLabel,1,n,k_init,n);  % transform currLabel into indicator matrix
m = X*(E*spdiags(1./sum(E,1)',0,k_init,k_init)); 
%n_m=size(m,2);
%disp(['n_m ',num2str(n_m)]);
end

% Kmeans++ initialization
function m = initSeedsKmeansPP(X, k)
[d,n] = size(X);
m = zeros(d,k);
v = inf(1,n);
m(:,1) = X(:,ceil(n*rand));
for i = 2:k
    Y = abs(bsxfun(@minus,X,m(:,i-1)));
    dd = sum(Y);
    v = cumsum(min(v,dd));
    m(:,i) = X(:,find(rand < v/v(end),1));
end
end