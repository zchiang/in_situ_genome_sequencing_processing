%{
Identify chromosome territory inliers based on spatial density.

Please note that the inlier labels in the published datatable may vary
slightly compared to the output of this script, as those labels were
assigned using this script on an earlier and slighly less curated 
version of the datatable.
%}

%Read the data
data = readtable('Table_S2_embryo_data_table.csv');

%%

%Generate list of embryos
embryoLabels = unique(data.embryo_id);
nEmbryo = numel(embryoLabels);
inlier = ones(size(data,1),1)*-1;

for ii = 1:nEmbryo %loop through embryos
    
    tmp = (data.embryo_id == ii);
    cells = max(data.cell_id(tmp));
    
    for jj = 1:cells %loop through cells
        
        tmp = (data.embryo_id == ii) & (data.cell_id == jj);
        chrs = max(data.chr(tmp));
        
        %Get spatial coordinates of amplicons
        x = data.x_um_abs(tmp);
        y = data.y_um_abs(tmp);
        z = data.z_um_abs(tmp);
        
        %Use data to estimate distance cutoff
        dist = doParwiseDistP(x,y,z);
        dist(dist == 0) = 50000;
        dist = mean(min(dist));
        
        
        for kk = 1:chrs %loop through chromsomes
            
            data_ind = (data.embryo_id == ii) & (data.cell_id == jj) & (data.chr == kk);
            cond = (sum(data_ind)> 10); %Don't attempt outlier removal if less than 8 points
            
            if ~cond
                disp('Did not qualify')
            else
                
                x = data.x_um_abs(data_ind);
                y = data.y_um_abs(data_ind);
                z = data.z_um_abs(data_ind);
                
                coord = [x';y';z'];
                clust_lab = data.cluster(data_ind); %if cluster labels exist use genomic distance for filtering
                
                seed = 2;
                
                if max(clust_lab)>1
                    ptsC = zeros(size(clust_lab));
                    
                    for ll = 1:max(clust_lab)
                        clust_ind_ll = (data.cluster(data_ind) == ll);
                        coord_ll = coord(:,clust_ind_ll);
                        [idx,ptsC_ll,~] = dbscan(coord_ll,3.5*dist,seed);
                        
                        if sum(clust_ind_ll)>2
                            ptsC(clust_ind_ll) = ptsC_ll;
                        else
                            ptsC(clust_ind_ll) = -1;
                        end
                        
                    end
                else
                    [idx,ptsC,~] = dbscan(coord,3.5*dist,seed);
                end
               
                
                
                ptsC = logical(ptsC);
                inlier(data_ind) = ptsC;
                
                x_o = x(~ptsC);
                y_o = y(~ptsC);
                z_o = z(~ptsC);
                
            end
            
            
            
        end
    end
end

%% Write labled data
data = [data,inlier];
writetable(data,'clustered_data.csv');