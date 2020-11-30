function [distP,distG] = genomicAndPhysicalDist(data, jj, kk, ll)

if  ~isempty(jj) & ~isempty(kk) & ~isempty(ll)
    
    celljj_chrkk_clusterll =  (data.cell == jj) & (data.hg38_chr == kk) & (data.cluster == ll);
    data_jj_kk_ll = data(celljj_chrkk_clusterll,:);

else
    data_jj_kk_ll = data;
end


    

distP = doParwiseDistP(data_jj_kk_ll.x_um,data_jj_kk_ll.y_um,data_jj_kk_ll.z_um);
distG = doParwiseDistG(data_jj_kk_ll.hg38_pos);

upperTri = logical(triu(distP,1));
distP= distP(upperTri);
distG = distG(upperTri);

end
