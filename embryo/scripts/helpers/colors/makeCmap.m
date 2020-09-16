function cmap = makeCmap(name,color)
% make colormap
cmap = color.(name);
cmap_m = zeros(length(cmap),3);
for j = 1:length(cmap)
    cmap_m(j,:) = rgbconv(cmap{j});
end

% expand color map
n = 50;  %%% was 5
cLen = length(cmap)-1;
cmap = zeros(cLen*n,3);
idxE = 1:n:length(cmap);

% start colors
for j = 1:cLen
    r1 = cmap_m(j,1); r2 = cmap_m(j+1,1); rInt = linspace(r1,r2,n);cmap(idxE(j):idxE(j)+n-1,1) = rInt;
    g1 = cmap_m(j,2); g2 = cmap_m(j+1,2); gInt = linspace(g1,g2,n);cmap(idxE(j):idxE(j)+n-1,2) = gInt;
    b1 = cmap_m(j,3); b2 = cmap_m(j+1,3); bInt = linspace(b1,b2,n);cmap(idxE(j):idxE(j)+n-1,3) = bInt;
end
end