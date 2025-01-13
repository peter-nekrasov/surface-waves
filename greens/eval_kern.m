function M = kernmat(src,targ,func,inds,corrs)
% Evaluates kernels and adds corrections to the source part based on inds

    kerns = func(src,targ);

    try 
        h = norm(targ(:,1) - targ(:,2)); % probably a better way to 
    catch
        h = norm(src(:,1) - src(:,2)); % probably a better way to 
    end

    if (numel(kerns) ~= numel(inds)) | (numel(inds) ~= numel(corrs))
        error('output of func must be the same length as inds and corrs')
    end

    M = cell(1,numel(kerns));

    for ii = 1:numel(kerns)
        kern = kerns{ii};
        kern = kern*h*h;
        ind = inds{ii};
        corr = corrs{ii};
        kern(ind) = kern(ind) + corr;
        M{ii} = kern;
    end

end