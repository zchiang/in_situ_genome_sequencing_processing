function[] = expt002_fov(exp,condition,fov,primers,ligations)

exp = "expt002";
condition = "well2";
fov = 3;

primers = ["N","N","N","N","N-1","N-1","N-1","N-1","N-2","N-2","N-2","N-2","N-3","N-3","N-3","N-3","N-4","N-4","N-4"];
ligations = ["L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","LL","LLL","LLLL"];

%%

tic
clearvars -except exp condition fov primers ligations
addpath(genpath('bin'));
colors = distinguishable_colors(101); colors(4,:) = [];

hyb_dir = sprintf('X:\\acp/191218_embryo_seq_full_run/%s/hyb',condition);
seq_dir =  sprintf('X:\\acp/191218_embryo_seq_full_run/%s',condition);
stain_dir = sprintf('X:\\acp/191218_embryo_seq_full_run/%s/stain',condition);

hyb_reader = bfGetReader(sprintf('%s/Experiment_*_F%02d.ims',hyb_dir,fov));
seq_reader = bfGetReader(sprintf('%s/%s/%s/Experiment_*_F%02d.ims',seq_dir,primers(1),ligations(1),fov));
stain_reader = bfGetReader(sprintf('%s/Experiment_*_F%02d.ims',stain_dir,fov));

hyb_xlen = hyb_reader.getSizeX;
hyb_ylen = hyb_reader.getSizeY;
hyb_zlen = hyb_reader.getSizeZ;

fov_xlen = seq_reader.getSizeX;
fov_ylen = seq_reader.getSizeY;
fov_zlen = seq_reader.getSizeZ;

num_channels = seq_reader.getSizeC;
num_cycles = size(primers,2);

hyb_channel = 1;
dapi_channel = 2;

disp(sprintf('%s: Processing %s FOV %d',sec2time(toc),condition,fov))

%% Create directories

con_dir = sprintf('%s/%s',exp,condition);
if ~exist(con_dir, 'dir') mkdir(con_dir), end

fov_dir = sprintf('%s/fov%03d',con_dir,fov);
if ~exist(fov_dir, 'dir') mkdir(fov_dir), end

%% Load full hyb and DAPI stacks

hyb_stack = uint16(zeros(hyb_xlen,hyb_ylen,hyb_zlen));
dapi_stack = uint16(zeros(hyb_xlen,hyb_ylen,hyb_zlen));

for z=1:hyb_zlen
    hyb_stack(:,:,z) = readPlane(hyb_reader,1,z,hyb_channel,1);
    dapi_stack(:,:,z) = readPlane(hyb_reader,1,z,dapi_channel,1);
end

disp(sprintf('%s: Loaded hyb and DAPI stacks',sec2time(toc)))

%% FOR FIXING SEGMENTATION
%{
%cell_bounds = dlmread(sprintf("%s/%s/fov%03d/cell_bounds.txt",exp,condition,fov));

figure;
for z=140:205
    imshow(dapi_stack(400:700,1350:1650,z),[])
end
%}
%% Segment nuclei using DAPI
%{
bw_stack = zeros(size(dapi_stack));


dapi_prc_thresh = 99.5;
thresh  = prctile(dapi_stack(:),dapi_prc_thresh);

%thresh = zeros(hyb_zlen,1);
%look_behind = 0;
%look_ahead = 0;

%figure;
for z=1:hyb_zlen
    
    %image = uint16(dapi_stack(:,:,z:min(z+look_ahead,hyb_zlen)));
    %image = uint16(dapi_stack(:,:,max(1,z-look_behind):min(z+look_ahead,hyb_zlen)));
    %thresh(z) = prctile(image(:),dapi_prc_thresh);

    bw_stack(:,:,z) = dapi_stack(:,:,z) > thresh;%thresh(z);
    bw_stack(:,:,z) = imopen(bw_stack(:,:,z),strel('disk',2));
    
    %imshowpair(bw_stack(:,:,z),capImage(dapi_stack(:,:,z),99,'prc'),'montage')
end

labels = bwlabeln(bw_stack);
regions = regionprops(labels, 'Area');
labels_filt = zeros(size(bw_stack));

vol_thresh = 100000;

for i=1:size(regions)
    if regions(i).Area > vol_thresh
        labels_filt(labels==i)=i;
    end
end
      
bw_stack = bwlabeln(labels_filt);
regions = regionprops(bw_stack, 'BoundingBox');

disp(sprintf('%s: Segmented nuclei',sec2time(toc)))

%% Crop based on nuclei segmentation

pad_xy = 200;
pad_z = 5;

nuclei_bounds = reshape([regions.BoundingBox],6,[])';

seg_min_x = floor(max(min(nuclei_bounds(:,1))-pad_xy,1));
seg_max_x = ceil(min(max(nuclei_bounds(:,1)+nuclei_bounds(:,4))+pad_xy,hyb_xlen));
seg_min_y = floor(max(min(nuclei_bounds(:,2))-pad_xy,1));
seg_max_y = ceil(min(max(nuclei_bounds(:,2)+nuclei_bounds(:,5))+pad_xy,hyb_ylen));
seg_min_z = floor(max(min(nuclei_bounds(:,3))-pad_z,1));
seg_max_z = ceil(min(max(nuclei_bounds(:,3)+nuclei_bounds(:,6))+pad_z,hyb_zlen));

bounds_3d = [seg_min_x seg_min_y seg_min_z seg_max_x-seg_min_x seg_max_y-seg_min_y seg_max_z-seg_min_z];

crop_dapi_stack = imcrop_xyz(dapi_stack,bounds_3d);
crop_hyb_stack = imcrop_xyz(hyb_stack,bounds_3d);
crop_bw_stack = imcrop_xyz(bw_stack,bounds_3d);

disp(sprintf('%s: Cropped hyb and DAPI based on segmentation',sec2time(toc)))

%% Segmentation video

crop_seg_stack = label2rgb3d(crop_bw_stack,colors);

video = VideoWriter(sprintf('figures/%s/dapi_vs_seg/%s_fov%03d',exp,condition,fov),'MPEG-4');
video.FrameRate = 15;

open(video);
for z=1:size(crop_seg_stack,3)
    f = figure('visible','off');
    imshowpair(capImage(crop_dapi_stack(:,:,z),99,'prc'),squeeze(crop_seg_stack(:,:,z,:)),'montage'); 
    text(25,25,sprintf('%d',z),'Color','Red')
    writeVideo(video,getframe(gcf))
end
close(video)
clear crop_dapi_stack crop_hyb_stack crop_bw_stack crop_seg_stack labels labels_filt

disp(sprintf('%s: Wrote segmentation video to file',sec2time(toc)))

%% Save cell bounds

pad_xy = 50;
pad_z = 3;
cell_bounds = [];

for cell=1:size(regions,1)
    
    cell_bounds(cell,:) = regions(cell).BoundingBox;
    cell_bounds(cell,:) = [max(0,cell_bounds(cell,1)-pad_xy) max(0,cell_bounds(cell,2)-pad_xy) max(0,cell_bounds(cell,3)-pad_z) ...
        min(cell_bounds(cell,4)+pad_xy*2,hyb_xlen-cell_bounds(cell,1)-pad_xy) min(cell_bounds(cell,5)+pad_xy*2,hyb_ylen-cell_bounds(cell,2)-pad_xy) ...
        min([cell_bounds(cell,6)+pad_z*2,hyb_zlen-cell_bounds(cell,3)-pad_z,fov_zlen-cell_bounds(cell,3)-pad_z])]; 
end

dlmwrite(sprintf("%s/%s/fov%03d/cell_bounds.txt",exp,condition,fov),cell_bounds);

return
%}
%% Load cell bounds

cell_bounds = dlmread(sprintf("%s/%s/fov%03d/cell_bounds.txt",exp,condition,fov));

%% Save DAPI and hyb for each cell
%{
warning('off')

for cell=1:size(cell_bounds,1)

    cell_dir = sprintf('%s/cell%03d',fov_dir,cell);
    if ~exist(cell_dir, 'dir') mkdir(cell_dir), end
    if ~exist(cell_dir+"/offset", 'dir') mkdir(cell_dir+"/offset"), end

    write_3d_tif(sprintf('%s/offset/hyb.tif',cell_dir),imcrop_xyz(hyb_stack,cell_bounds(cell,:)));
    write_3d_tif(sprintf('%s/offset/dapi.tif',cell_dir),imcrop_xyz(dapi_stack,cell_bounds(cell,:)));

    disp(sprintf('%s: Saved DAPI + hyb stacks for cell %d',sec2time(toc),cell))
    
end

%%

warning('off')

curr_stack = uint16(zeros(fov_xlen,fov_ylen,fov_zlen,num_channels));
offset_stack = uint16(zeros(size(curr_stack)));
proj_stack = uint16(zeros(fov_xlen,fov_ylen,num_cycles,2));

%cycle = 1;
for cycle=1:num_cycles
    
    if condition == "well1" & primers(cycle) == "N" & ligations(cycle) == "LLL"
        if fov >= 0 & fov <= 9
            seq_reader = bfGetReader(sprintf('%s/%s/%s/first_10_fovs/Experiment_*_F%02d.ims',seq_dir,primers(cycle),ligations(cycle),fov));
        elseif fov >= 9
            seq_reader = bfGetReader(sprintf('%s/%s/%s/Experiment_*_F%02d.ims',seq_dir,primers(cycle),ligations(cycle),fov-9));
        end
    elseif condition == "well2" & primers(cycle) == "N" & ligations(cycle) == "LL"
        seq_reader = bfGetReader(sprintf('%s/%s/%s/Experiment_*_F%02d.ims',seq_dir,primers(cycle),ligations(cycle),fov-2));
    else
        seq_reader = bfGetReader(sprintf('%s/%s/%s/Experiment_*_F%02d.ims',seq_dir,primers(cycle),ligations(cycle),fov));
    end
    
    for channel=1:num_channels
        for z=1:fov_zlen
            curr_stack(:,:,z,channel) = readPlane(seq_reader,1,z,channel,1);
        end
    end
    
    if condition == "scratch1" & primers(cycle) == "N" & ligations(cycle) == "L"
        for z=1:fov_zlen
            seq_reader = bfGetReader(sprintf('%s/%s/%s/reacquire_c4/Experiment_*_F%01d.ims',seq_dir,primers(cycle),ligations(cycle),fov));
            curr_stack(:,:,z,4) = readPlane(seq_reader,1,z,1,1);
        end
    end
    
    proj_stack(:,:,cycle,1) = max(max(curr_stack,[],3),[],4);
    
    disp(sprintf('%s: Loaded images for cycle %d',sec2time(toc),cycle))
    
    offset_stack = simple_offset_xyz(hyb_stack,curr_stack);
    proj_stack(:,:,cycle,2) = max(max(offset_stack,[],3),[],4);
    
    disp(sprintf('%s: Offset cycle %d',sec2time(toc),cycle))
    
    for cell=1:size(cell_bounds,1)
        cell_dir = sprintf('%s/cell%03d',fov_dir,cell);
        if ~exist(cell_dir+"/offset", 'dir') mkdir(cell_dir+"/offset"), end
        for channel=1:num_channels
            filename = sprintf('%s/offset/cy%02d_ch%02d.tif',cell_dir,cycle,channel);
            write_3d_tif(filename, imcrop_xyz(offset_stack(:,:,:,channel),cell_bounds(cell,:)));
        end
    end
     
    disp(sprintf('%s: Saved cycle %d cell stacks to file',sec2time(toc),cycle))
    
end

clearvars curr_stack offset_stack

%% Save projected stack

write_3d_tif(sprintf('%s/proj_before.tif',fov_dir),proj_stack(:,:,:,1))
write_3d_tif(sprintf('%s/proj_after.tif',fov_dir),proj_stack(:,:,:,2))
%}
%% Load stain images and offset using DAPI

num_stain_channels = stain_reader.getSizeC;
stain_xlen = stain_reader.getSizeX;
stain_ylen = stain_reader.getSizeY;
stain_zlen = stain_reader.getSizeZ;
stain_stacks = zeros(stain_xlen,stain_ylen,stain_zlen,num_stain_channels);

for channel=1:num_stain_channels
    for z=1:stain_zlen
        stain_stacks(:,:,z,channel) = readPlane(stain_reader,1,z,channel,1);
    end
end

disp(sprintf('%s: Loaded stain stacks',sec2time(toc)))

% Offset stain stacks

offset_stain_stacks = simple_offset_xyz_applied(dapi_stack,stain_stacks(:,:,:,3),stain_stacks);

disp(sprintf('%s: Offset stain stacks',sec2time(toc)))

%% Save stain stacks

for cell=1:size(cell_bounds,1)

    cell_dir = sprintf('%s/cell%03d',fov_dir,cell);
    if ~exist(cell_dir, 'dir') mkdir(cell_dir), end

    write_3d_tif(sprintf('%s/offset/cenpa.tif',cell_dir),imcrop_xyz(stain_stacks(:,:,:,1),cell_bounds(cell,:)));
    write_3d_tif(sprintf('%s/offset/lamin.tif',cell_dir),imcrop_xyz(stain_stacks(:,:,:,2),cell_bounds(cell,:)));

    disp(sprintf('%s: Saved cenpa and lamin stacks for cell %d',sec2time(toc),cell))
    
end

%% Offset visualization (with stain)
%{
figure('visible','off');
p = tight_subplot(2,num_cycles+1,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles+1
    
        if cycle == num_cycles+1
            
            axes(p(cycle)); imshowpair(capImage(max(dapi_stack,[],3),99,'prc'),capImage(max(stain_stacks(:,:,:,3),[],3),99,'prc'))
            title(sprintf('DAPI - Stain DAPI'));   
            axes(p(num_cycles+1+cycle)); imshowpair(capImage(max(dapi_stack,[],3),99,'prc'),capImage(max(offset_stain_stacks(:,:,:,3),[],3),99,'prc'))
            
        else
    
            axes(p(cycle)); imshowpair(capImage(max(hyb_stack,[],3),99,'prc'),capImage(proj_stack(:,:,cycle,1),99,'prc'))
            title(sprintf('Hyb - Cycle %d',cycle));  
            axes(p(num_cycles+1+cycle)); imshowpair(capImage(max(hyb_stack,[],3),99,'prc'),capImage(proj_stack(:,:,cycle,2),99,'prc'));
            
        end      
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2 5];

fig_dir = sprintf('figures/%s/offsets',exp);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d.png',fig_dir,condition,fov));


%% Offset visualization

figure('visible','off');
p = tight_subplot(2,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles

    axes(p(cycle)); imshowpair(capImage(max(hyb_stack,[],3),99,'prc'),capImage(proj_stack(:,:,cycle,1),99,'prc'))
	title(sprintf('Hyb - Cycle %d',cycle));  
	axes(p(num_cycles+cycle)); imshowpair(capImage(max(hyb_stack,[],3),99,'prc'),capImage(proj_stack(:,:,cycle,2),99,'prc'));
            

end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2 5];

fig_dir = sprintf('figures/%s/offsets',exp);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d.png',fig_dir,condition,fov));

return
%}