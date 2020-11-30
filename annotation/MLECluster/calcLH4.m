function lkh = calcLH4(data, kk, idx, chr)

train = [chr(kk).distG, chr(kk).distP];
train = sortrows(train,1);

%Region for normalization
LG = 0-0.01;
UG = chr(kk).length;
LP = 0;
UP = max(chr(kk).distP)+1;

[distP_construct,distG_construct] = genomicAndPhysicalDist(data, [], [], []);

for ii = 1:max(idx) %loop through clusters
    
    %Get distances for each segment along path
    data_ii = data(idx==ii,:);
    
    [dP,dG] = getDistancesForEachSegment(data_ii);
    [f,xi] = ksdensity([chr(kk).distG,chr(kk).distP],[distG_construct,distP_construct],'Support',[LG,LP;UG,UP],'BoundaryCorrection','reflection');
    [c,ia,ib] = intersect([dG,dP],xi,'rows');
    
    lkh_tmp(ii) = sum(log(f(ib)));
    
end

lkh = sum(lkh_tmp);


end
