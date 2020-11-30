function f = plotChromoRandomWalk2(data,jj,kk,idx,chr)

% Rows of data matrix corresponding to cell jj
celljj = (data.cell_id == jj);

% Rows of data matrix corresponding cell jj and chromosome kk
celljj_chrkk = (data.cell_id == jj) & (data.hg38_chr == kk);

%Relevant section of data matrix
data_jj_kk = data(celljj_chrkk,:);
data_jj = data(celljj,:);

%plot all reads in light gray
scatter1 = scatter3(data_jj.x_um,data_jj.y_um,data_jj.z_um,15, 'filled', 'MarkerFaceColor',[.9,.9,.9]);
scatter1.MarkerFaceAlpha = .4;
hold on;
view(45,70);
title([' Cell ID: ',num2str(jj),' Chr: ', num2str(kk)]);

if ~isempty(data_jj_kk)
    
    %List of genomic positions for chromosome pair
    cpos = data_jj_kk.hg38_pos; 
    
    %Which row of the colormap genomic position corresponds to
    scaled = ceil((cpos/chr(kk).length)*256);
    
    %Create colormap with 265 intervals;
    cc = parula(256);
    ccc = cc(scaled,:);

 
for ii = 1:max(idx)   
    
    data_ = data_jj_kk(idx == ii,:);
    data_ordered = [data_.x_um,data_.y_um,data_.z_um,data_.hg38_pos];
    data_ordered = sortrows(data_ordered,4);
    
   
    for jj = 1:size(data_ordered,1)-1
        
        x1 = data_ordered(jj,1);
        x2 = data_ordered(jj+1,1);
        
        y1 = data_ordered(jj,2);
        y2 = data_ordered(jj+1,2);
        
        z1 = data_ordered(jj,3);
        z2 = data_ordered(jj+1,3);
        
        %Calculate index in colormap
        cmapInd = ceil((data_ordered(jj,4)/chr(kk).length)*256);
        plot3([x1,x2],[y1,y2],[z1,z2], 'Color', cc(cmapInd,:));
        
    end
    
    
end

    scatter3(data_jj_kk.x_um,data_jj_kk.y_um,data_jj_kk.z_um,15, ccc,'filled');
    xlabel('x (um)');
    ylabel('y (um)');
    zlabel('z (um)');
    view(45,70)
  
end
end

