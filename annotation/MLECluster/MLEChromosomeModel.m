%{
Cluster reads by spatial density. 

Please note that the cluster labels in the published datatable may vary
slightly compared to the output of this script, as those labels were
assigned using this script on an earlier and slighly less curated 
version of the datatable.
%}

%% Load the data

data = readtable('Table_S1_pgp1_data_table.csv');

%% Generate training data

%The training data is an set of clustered data where groupings can be
%established by eye. This is used to generate the emperical probability
%distribution for genomic vs. physical distance.

chr = getTrainingData('trainingData.txt');

%% Fit the model

nCell = max(data.cell_id);
nChr = max(data.hg38_chr);

jj = 1; %Cell ID
kk = 1; %Chromosome number


%Logical - marks rows of specified cell/chromosome
celljj = (data.cell_id == jj);
celljj_chrkk = (data.cell_id == jj) & (data.hg38_chr == kk);

%Data matrix for cell/row
data_jj = data(celljj,:);
data_jj_kk = data(celljj_chrkk,:);


if ~isempty(data_jj_kk)
    
    %Outlier removal using dbscan
    matDB = [data_jj_kk.x_um, data_jj_kk.y_um, data_jj_kk.z_um];
    [C, ptsC, centres] = dbscan(matDB',3,2);
    idx = zeros(size(ptsC));
    outliers = (ptsC == 0);
    
    %Initialize clustering
    data_to_cluster = data_jj_kk;
    data_to_cluster(ptsC==0,:) = [];
    
    if size(data_to_cluster,1) > 1
        idx_tmp = kmeans([data_to_cluster.x_um,data_to_cluster.y_um,data_to_cluster.z_um],2,'Replicates',10);
        idx(~outliers) = idx_tmp;
        numClusterPoints = numel(idx_tmp);
        
        %Show initial grouping
        figure(1)
        clf;
        subplot(1,2,1)
        plotChromoRandomWalkColorByDist(data,jj,kk,idx,chr)
        axis equal;
        title('Initialization')
        
        if (numClusterPoints > 9) && sum(idx == 1) > 1 && sum(idx == 2) > 1
            
            %Calculate liklihood associated with initial model
            lkh = calcLH4(data_jj_kk,kk,idx,chr);
            
            %Inatialize variables
            n = 0; %number of iterations
            zero_count = 0; %number of successive iterations where liklihood has not changed
            lkh_track = lkh; %trak liklhood after each iteration
            idx_init = idx;
            
            %Liklihood maximization
            %Successivly move each point to the other group
            %If liklihood increases, keep new grouping
            %Stop when moving any of the points in each group does not increase
            %liklihooid
            
            while (zero_count < numel(idx))
                
                idx_tmp = idx;
                
                %n_ defines the which point will be moved to other group
                n_ = mod(n,numel(idx))+1;
                
                
                % Move point to other group if group has at least one point
                if idx(n_) == 1 & sum(idx == 1) > 2
                    idx_tmp(n_) = 2;
                elseif idx(n_)==2 & sum(idx == 2) > 2
                    idx_tmp(n_)=1;
                end
                
                %Re-calculate liklihood
                lkh_tmp = calcLH4(data_jj_kk,kk,idx_tmp,chr);
                
                %If liklihood increased, keep new group assignment
                if lkh_tmp > lkh
                    
                    idx = idx_tmp;
                    lkh = lkh_tmp;
                    zero_count = 0;
                    
                else
                    %count how many cycles it's been since liklihood increased
                    zero_count = zero_count+1;
                end
                
                %Track liklihood across cycles
                lkh_track = [lkh_track,lkh];
                
                %count number of cycles
                n = n+1;
                
                %Plot grouping at end of cycle
                figure(3)
                clf;
                subplot(1,2,1)
                plotChromoRandomWalkColorByDist(data,jj,kk,idx,chr)
                colorbar
                
                subplot(1,2,2)
                plot(lkh_track);
                xlabel('Cycle number')
                ylabel('Log Liklihood')
                
            end
            
            figure(1)
            subplot(1,2,2)
            plotChromoRandomWalkColorByDist(data,jj,kk,idx,chr)
            axis equal;
            title('MLE Grouping')
            
            %Ask the user to check the grouping by eye
            accept = input('Accept this grouping? ')
            
            
            if accept
                data.mle_cluster(celljj_chrkk) = idx;
            end
            
        end
        
    end
    
end


%Write the data table with mle cluster assignments
writetable(data,'clustered_data.txt');









