% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Processes a single field of view of in situ genome sequencing data

function[] = expt002_fov(exp,condition,fov,primers,ligations,img_dir)
%% Can manually set parameters here

%exp = "expt002";
%condition = "well2";
%fov = 3;

%primers = ["N","N","N","N","N-1","N-1","N-1","N-1","N-2","N-2","N-2","N-2","N-3","N-3","N-3","N-3","N-4","N-4","N-4"];%
%ligations = ["L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","LL","LLL","LLLL"];

%img_dir = 'X:\\acp/191218_embryo_seq_full_run/';

%% Clear environment, start time, add path

clearvars -except exp condition fov primers ligations img_dir
tic
addpath(genpath('bin'));
colors = distinguishable_colors(101); colors(4,:) = []; % remove black

% Create directories
con_dir = sprintf('%s/%s',exp,condition);
if ~exist(con_dir, 'dir') mkdir(con_dir), end

fov_dir = sprintf('%s/fov%03d',con_dir,fov);
if ~exist(fov_dir, 'dir') mkdir(fov_dir), end

%% Load image readers and get dimensions

hyb_dir = sprintf('%s/%s/hyb',img_dir,condition);
seq_dir =  sprintf('%s/%s',img_dir,condition);
stain_dir = sprintf('%s/%s/stain',img_dir,condition);

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

% manually set hyb reader channel ids
hyb_channel = 1;
dapi_channel = 2;

disp(sprintf('%s: Processing %s FOV %d',sec2time(toc),condition,fov))

%% Load full hyb and DAPI stacks

hyb_stack = uint16(zeros(hyb_xlen,hyb_ylen,hyb_zlen));
dapi_stack = uint16(zeros(hyb_xlen,hyb_ylen,hyb_zlen));

for z=1:hyb_zlen
    hyb_stack(:,:,z) = readPlane(hyb_reader,1,z,hyb_channel,1);
    dapi_stack(:,:,z) = readPlane(hyb_reader,1,z,dapi_channel,1);
end

disp(sprintf('%s: Loaded hyb and DAPI stacks',sec2time(toc)))

%% Segment nuclei using DAPI

dapi_prc_thresh = 99.5; % adjust threshold here
thresh  = prctile(dapi_stack(:),dapi_prc_thresh);
bw_stack = zeros(size(dapi_stack));

% segment each z-slice

%figure; % for debugging
for z=1:hyb_zlen
    
    bw_stack(:,:,z) = dapi_stack(:,:,z) > thresh;
    bw_stack(:,:,z) = imopen(bw_stack(:,:,z),strel('disk',2)); % morphological opening, reduces noise
    
    %imshowpair(bw_stack(:,:,z),capImage(dapi_stack(:,:,z),99,'prc'),'montage') % for debugging
    
end

% find area of segmented objects

labels = bwlabeln(bw_stack);
regions = regionprops(labels, 'Area');
labels_filt = zeros(size(bw_stack));

% filter objects less than volume threshold

vol_thresh = 100000;
for i=1:size(regions)
    if regions(i).Area > vol_thresh
        labels_filt(labels==i)=i;
    end
end

% re-label objects after filtering
      
bw_stack = bwlabeln(labels_filt);
regions = regionprops(bw_stack, 'BoundingBox');

disp(sprintf('%s: Segmented nuclei',sec2time(toc)))

%% Save cell bounds

% set padding for cell bounds here

pad_xy = 50;
pad_z = 3;
cell_bounds = [];

% correct cell bounds so that they do not exceed image bounds

for cell=1:size(regions,1)
    
    cell_bounds(cell,:) = regions(cell).BoundingBox;
    cell_bounds(cell,:) = [max(0,cell_bounds(cell,1)-pad_xy) max(0,cell_bounds(cell,2)-pad_xy) max(0,cell_bounds(cell,3)-pad_z) ...
        min(cell_bounds(cell,4)+pad_xy*2,hyb_xlen-cell_bounds(cell,1)-pad_xy) min(cell_bounds(cell,5)+pad_xy*2,hyb_ylen-cell_bounds(cell,2)-pad_xy) ...
        min([cell_bounds(cell,6)+pad_z*2,hyb_zlen-cell_bounds(cell,3)-pad_z,fov_zlen-cell_bounds(cell,3)-pad_z])]; 
end

dlmwrite(sprintf("%s/%s/fov%03d/cell_bounds.txt",exp,condition,fov),cell_bounds);

%% Load cell bounds (can skip to this step after performing segmentation once)

cell_bounds = dlmread(sprintf("%s/%s/fov%03d/cell_bounds.txt",exp,condition,fov));

%% Save DAPI and hyb for each cell

warning('off')

for cell=1:size(cell_bounds,1)

    cell_dir = sprintf('%s/cell%03d',fov_dir,cell);
    if ~exist(cell_dir, 'dir') mkdir(cell_dir), end
    if ~exist(cell_dir+"/offset", 'dir') mkdir(cell_dir+"/offset"), end

    write_3d_tif(sprintf('%s/offset/hyb.tif',cell_dir),imcrop_xyz(hyb_stack,cell_bounds(cell,:)));
    write_3d_tif(sprintf('%s/offset/dapi.tif',cell_dir),imcrop_xyz(dapi_stack,cell_bounds(cell,:)));

    disp(sprintf('%s: Saved DAPI + hyb stacks for cell %d',sec2time(toc),cell))
    
end

%% Save cropped images for all cells

warning('off')

curr_stack = uint16(zeros(fov_xlen,fov_ylen,fov_zlen,num_channels));
offset_stack = uint16(zeros(size(curr_stack)));

%cycle = 1;
for cycle=1:num_cycles
    
    % hard code to deal with image acquisition/naming discrepancies
    
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
    
    % load all images
    
    for channel=1:num_channels
        for z=1:fov_zlen
            curr_stack(:,:,z,channel) = readPlane(seq_reader,1,z,channel,1);
        end
    end
    
    % more hard code
    
    if condition == "scratch1" & primers(cycle) == "N" & ligations(cycle) == "L"
        for z=1:fov_zlen
            seq_reader = bfGetReader(sprintf('%s/%s/%s/reacquire_c4/Experiment_*_F%01d.ims',seq_dir,primers(cycle),ligations(cycle),fov));
            curr_stack(:,:,z,4) = readPlane(seq_reader,1,z,1,1);
        end
    end
    
    disp(sprintf('%s: Loaded images for cycle %d',sec2time(toc),cycle))
    
    % offset in 3D using cross-correlation
    
    offset_stack = simple_offset_xyz(hyb_stack,curr_stack);
    
    disp(sprintf('%s: Offset cycle %d',sec2time(toc),cycle))
    
    % save offset images
    
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

clearvars curr_stack offset_stack % saves memory

%% Load stain images and offset using DAPI

stain_dapi_channel = 3;

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

offset_stain_stacks = simple_offset_xyz_applied(dapi_stack,stain_stacks(:,:,:,stain_dapi_channel),stain_stacks);

disp(sprintf('%s: Offset stain stacks',sec2time(toc)))

%% Save stain stacks

for cell=1:size(cell_bounds,1)

    cell_dir = sprintf('%s/cell%03d',fov_dir,cell);
    if ~exist(cell_dir, 'dir') mkdir(cell_dir), end

    write_3d_tif(sprintf('%s/offset/cenpa.tif',cell_dir),imcrop_xyz(stain_stacks(:,:,:,1),cell_bounds(cell,:)));
    write_3d_tif(sprintf('%s/offset/lamin.tif',cell_dir),imcrop_xyz(stain_stacks(:,:,:,2),cell_bounds(cell,:)));

    disp(sprintf('%s: Saved cenpa and lamin stacks for cell %d',sec2time(toc),cell))
    
end

end