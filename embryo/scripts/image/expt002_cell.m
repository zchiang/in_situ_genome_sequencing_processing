function[] = expt002_cell(exp,condition,fov,cell,run,umi_file,primers,ligations)
%% Manually set parameters

exp = "expt002";
%condition = "scratch1";
condition = "well2";
fov = 30;
cell = 2;
run = "a1";

primers = ["N","N","N","N","N-1","N-1","N-1","N-1","N-2","N-2","N-2","N-2","N-3","N-3","N-3","N-3","N-4","N-4","N-4"];
ligations = ["L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","LL","LLL","LLLL"];

umi_file = "expt002_nova_all.mm10.final.txt";

%% Clear environment, start time, add path

clearvars -except exp condition fov embryo cell run umi_file primers ligations
tic
addpath(genpath('bin'));
colors = distinguishable_colors(101); colors(4,:) = [];

%primers = ["N","N","N","N","N-1","N-1","N-1","N-1","N-2","N-2","N-2","N-2","N-3","N-3","N-3","N-3","N-4","N-4","N-4"];
%ligations = ["L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","L","LL","LLL","LLLL","LL","LLL","LLLL"];

num_cycles = size(primers,2);
num_channels = 4;

disp(sprintf('%s: Processing cell %d from %s FOV %d',sec2time(toc),cell,condition,fov))

%% Create directories

fov_dir = sprintf('%s/%s/fov%03d',exp,condition,fov);
cell_dir = sprintf('%s//cell%03d',fov_dir,cell);

cell_reg_dir = sprintf('%s/registration',cell_dir);
if ~exist(cell_reg_dir, 'dir') mkdir(cell_reg_dir), end

%% Load cell stack

cell_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
cell_bounds = floor(cell_bounds(cell,:));

stack = zeros(cell_bounds(5),cell_bounds(4),cell_bounds(6),num_channels,num_cycles);

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/offset/cy%02d_ch%02d.tif',cell_dir,cycle,channel);
        stack(:,:,:,channel,cycle) = read_3d_tif(filename,cell_bounds(5),cell_bounds(4),cell_bounds(6));
    end
end

% Load hyb + DAPI stacks

hyb_stack = read_3d_tif(sprintf('%s/offset/hyb.tif',cell_dir),cell_bounds(5),cell_bounds(4),cell_bounds(6));
dapi_stack = read_3d_tif(sprintf('%s/offset/dapi.tif',cell_dir),cell_bounds(5),cell_bounds(4),cell_bounds(6));

disp(sprintf('%s: Loaded cell stack',sec2time(toc)))

%% Deconvolve with Gaussian filter

deconv_stack = zeros(size(stack));

for cycle=1:num_cycles
    for channel=1:num_channels
        for z=1:size(deconv_stack,3)
            high_pass_filter = imgaussfilt(stack(:,:,z,channel,cycle),2);
            high_pass_image = stack(:,:,z,channel,cycle) - high_pass_filter;
            deconv_stack(:,:,z,channel,cycle) = high_pass_image;
        end
    end
end

deconv_stack(deconv_stack<0) = 0;

% Deconvolve hyb stack

deconv_hyb_stack = zeros(size(hyb_stack));

for z=1:size(hyb_stack,3)
    high_pass_filter = imgaussfilt(hyb_stack(:,:,z),2);
	high_pass_image = hyb_stack(:,:,z) - high_pass_filter;
	deconv_hyb_stack(:,:,z) = high_pass_image;
end

disp(sprintf('%s: Deconvolved images with Gaussian filter',sec2time(toc)))

%% Cap high pixel values

cap = {};
cap_prctile = 99.95;
cap_stack = zeros(size(deconv_stack));

for cycle=1:num_cycles
    for channel=1:num_channels
        cap{channel,cycle} = prctile(reshape(deconv_stack(:,:,:,channel,cycle),[],1),cap_prctile);
        tmp_stack = deconv_stack(:,:,:,channel,cycle);
        
        pixels = find(tmp_stack>cap{channel,cycle});
        [sort_val sort_order] = sort(tmp_stack(pixels));       
                
        tmp_stack(tmp_stack>cap{channel,cycle}) = cap{channel,cycle};
        tmp_stack(pixels(sort_order)) = tmp_stack(pixels(sort_order))+(1:size(pixels,1))';
        cap_stack(:,:,:,channel,cycle) = tmp_stack;
        
    end
end

% Cap hyb stack

cap_hyb_stack = zeros(size(deconv_hyb_stack));

hyb_cap = prctile(deconv_hyb_stack(:),cap_prctile);
tmp_stack = deconv_hyb_stack(:,:,:,1);
        
pixels = find(tmp_stack>hyb_cap);
[sort_val sort_order] = sort(tmp_stack(pixels));       
                
tmp_stack(tmp_stack>hyb_cap) = hyb_cap;
tmp_stack(pixels(sort_order)) = tmp_stack(pixels(sort_order))+(1:size(pixels,1))';
cap_hyb_stack = tmp_stack;

disp(sprintf('%s: Capped high pixel values',sec2time(toc)));

%% Hyb + DAPI registration

min_overlap = 0.5;
reg_stack = zeros(size(cap_stack));
curr_cycle = zeros(size(cap_stack));
reg_stack(:,:,:,:,1) = cap_stack(:,:,:,:,1);
reg_tforms = {};

% Register hyb to cycle 1

[reg_dapi_stack reg_hyb_stack stats hyb_tform] = register_cycle(max(reg_stack(:,:,:,:,1),[],4),cap_hyb_stack,dapi_stack,min_overlap);
    
disp(sprintf('%s: Registered hyb and DAPI: %s',sec2time(toc),stats));

%% Nucleus segmentation (up here in case you want to use it to focus registration)
%{
bw_stack = zeros(size(reg_dapi_stack));

image = uint16(reg_dapi_stack);
thresh = multithresh(image,2);
thresh = thresh(2)-25;
    
%figure;
for z=1:size(bw_stack,3)

    bw_stack(:,:,z) = dapi_stack(:,:,z) > thresh;
    bw_stack(:,:,z) = imopen(bw_stack(:,:,z),strel('disk',2));
    %bw_stack(:,:,z) = imdilate(bw_stack(:,:,z),strel('disk',2));
    %bw_stack(:,:,z) = imfill(bw_stack(:,:,z),'holes');
    
    %imshowpair(bw_stack(:,:,z),capImage(dapi_stack(:,:,z),99,'prc'),'montage')
end

labels = bwlabeln(bw_stack);
regions = regionprops(labels, 'Area');
labels_filt = zeros(size(bw_stack));

for i=1:size(regions)
    if regions(i).Area == max([regions.Area])
        labels_filt(labels==i)=i;
    end
end
      
bw_stack = bwlabeln(labels_filt);
regions = regionprops(bw_stack, 'BoundingBox');

disp(sprintf('%s: Segmented nucleus',sec2time(toc)))
%}
%% Registration

for cycle=19:19%2:num_cycles
    
    if condition == "well2" & (fov  == 0 | fov == 1) & cycle == 2
        continue
    end
    
    curr_cycle = max(cap_stack(:,:,:,:,cycle),[],4);  
    [reg_stack(:,:,:,:,cycle) reg_cycle stats reg_tforms{cycle}] = register_cycle(max(reg_stack(:,:,:,:,1),[],4),curr_cycle,cap_stack(:,:,:,:,cycle),min_overlap);
    %[reg_stack(:,:,:,:,cycle) reg_cycle stats] = register_cycle(max(reg_stack(:,:,:,:,1),[],4).*bw_stack,curr_cycle.*bw_stack,cap_stack(:,:,:,:,cycle),min_overlap);
    
    disp(sprintf('%s: Cycle %d registration: %s',sec2time(toc),cycle,stats));
    
end

%% Color correct cycle 1 to hyb and save offsets

reg_stack2 = reg_stack;
correlations = zeros(num_channels+1,1);
tforms = {};

cycle = 1;
for channel=1:num_channels
    [reg_stack2(:,:,:,channel,cycle) correlations(channel) correlations(channel+1) tforms{channel}] = ...
        color_correct_initial(reg_hyb_stack,reg_stack2(:,:,:,:,cycle),channel);
    
end
disp(sprintf('%s: Cycle %d color correction: initial=%.03f, channel 1=%+.03f, channel 2=%+.03f, channel 3=%+.03f, channel 4=%+.03f, final=%.03f',sec2time(toc),cycle,...
    correlations(1),correlations(2)-correlations(1),correlations(3)-correlations(2),correlations(4)-correlations(3),correlations(5)-correlations(4),correlations(5)));

%% Secondary color correction

reg_stack2(:,:,:,:,2:end) = reg_stack(:,:,:,:,2:end);
correlations = zeros(num_channels+1,1);

for cycle=2:num_cycles
    
    if condition == "well2" & (fov  == 0 | fov == 1) & cycle == 2
        continue
    end
    
    for channel=1:num_channels
        [reg_stack2(:,:,:,channel,cycle) correlations(channel) correlations(channel+1)] = ...
            color_correct_heuristic(reg_stack2(:,:,:,:,1),reg_stack2(:,:,:,:,cycle),channel,tforms{channel});
        %[reg_stack2(:,:,:,channel,cycle) correlations(channel) correlations(channel+1)] = fine_cc_nucleus(reg_stack2(:,:,:,:,1).*bw_stack,reg_stack2(:,:,:,:,cycle).*bw_stack,reg_stack2(:,:,:,:,cycle),channel);

    end
    disp(sprintf('%s: Cycle %d color correction: initial=%.03f, channel 1=%+.03f, channel 2=%+.03f, channel 3=%+.03f, channel 4=%+.03f, final=%.03f',sec2time(toc),cycle,...
        correlations(1),correlations(2)-correlations(1),correlations(3)-correlations(2),correlations(4)-correlations(3),correlations(5)-correlations(4),correlations(5)));
end

%% Channel quantile normalization and zero matrix

norm_stack = zeros(size(reg_stack2));

for cycle=1:num_cycles
    
    tmp_stack = reshape(reg_stack2(:,:,:,:,cycle),cell_bounds(5)*cell_bounds(4)*cell_bounds(6),num_channels);
    quantile_norm = quantilenorm(tmp_stack);
    norm_stack(:,:,:,:,cycle) = reshape(quantile_norm,cell_bounds(5),cell_bounds(4),cell_bounds(6),num_channels);
    
    norm_stack(:,:,:,:,cycle) = norm_stack(:,:,:,:,cycle) - mean(reshape(norm_stack(:,:,:,:,cycle),[],1));
    
end

norm_stack(norm_stack<0) = 0;

disp(sprintf('%s: Performed quantile normaliation and zeroed matrix',sec2time(toc)))

%% Save visualizations

% All cycles

f = figure('visible','off');
p = tight_subplot(num_channels,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles
    for channel=1:num_channels   
        axes(p((channel-1)*num_cycles+cycle)); imshow(capImage(max(norm_stack(:,:,:,channel,cycle),[],3),95,'prc'),[]);
        title(sprintf('cycle %d, channel %d',cycle,channel))
    end
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2*size(norm_stack,2)/max(size(norm_stack,1),size(norm_stack,2)) 12*size(norm_stack,1)/max(size(norm_stack,1),size(norm_stack,2))];;

fig_dir = sprintf('figures/%s/all_stacks_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));

% Cycle overlap

f = figure('visible','off');
p = tight_subplot(4,num_cycles-1,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=2:num_cycles
    
    init_corr(cycle-1) = corr(reshape(max(stack(:,:,:,:,cycle-1),[],4),[],1),reshape(max(stack(:,:,:,:,cycle),[],4),[],1));
    deconv_corr(cycle-1) = corr(reshape(max(deconv_stack(:,:,:,:,cycle-1),[],4),[],1),reshape(max(deconv_stack(:,:,:,:,cycle),[],4),[],1));
    reg_corr(cycle-1) = corr(reshape(max(reg_stack2(:,:,:,:,cycle-1),[],4),[],1),reshape(max(reg_stack2(:,:,:,:,cycle),[],4),[],1));
    norm_corr(cycle-1) = corr(reshape(max(norm_stack(:,:,:,:,cycle-1),[],4),[],1),reshape(max(norm_stack(:,:,:,:,cycle),[],4),[],1));
    
    axes(p(cycle-1)); imshowpair(capImage(max(max(stack(:,:,:,:,cycle-1),[],4),[],3),95,'prc'),capImage(max(max(stack(:,:,:,:,cycle),[],4),[],3),95,'prc')); 
    title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,init_corr(cycle-1)));
    
    axes(p(num_cycles+cycle-2)); imshowpair(capImage(max(max(deconv_stack(:,:,:,:,cycle-1),[],4),[],3),95,'prc'),capImage(max(max(deconv_stack(:,:,:,:,cycle),[],4),[],3),95,'prc'));
    title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,deconv_corr(cycle-1)));
    
    axes(p(num_cycles*2+cycle-3)); imshowpair(capImage(max(max(reg_stack2(:,:,:,:,cycle-1),[],4),[],3),95,'prc'),capImage(max(max(reg_stack2(:,:,:,:,cycle),[],4),[],3),95,'prc'));
    title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,reg_corr(cycle-1)));
    
    axes(p(num_cycles*3+cycle-4)); imshowpair(capImage(max(max(norm_stack(:,:,:,:,cycle-1),[],4),[],3),95,'prc'),capImage(max(max(norm_stack(:,:,:,:,cycle),[],4),[],3),95,'prc'));
    title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,norm_corr(cycle-1)));
      
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles-1)*2*size(norm_stack,2)/max(size(norm_stack,1),size(norm_stack,2)) 12*size(norm_stack,1)/max(size(norm_stack,1),size(norm_stack,2))];

fig_dir = sprintf('figures/%s/cycle_overlap_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));

disp(sprintf('%s: Saved visualizations',sec2time(toc)))

%% Save stacks

% Save normalized stack

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/cy%02d_ch%02d.tif',cell_reg_dir,cycle,channel);
        write_3d_tif(filename, norm_stack(:,:,:,channel,cycle));
    end
end

% Save hyb + DAPI and segmentation

write_3d_tif(sprintf('%s/reg_hyb.tif',cell_dir),reg_hyb_stack);
write_3d_tif(sprintf('%s/reg_dapi.tif',cell_dir),reg_dapi_stack);
%write_3d_tif(sprintf('%s/bw.tif',cell_dir),bw_stack);

disp(sprintf('%s: Saved normalized stacks',sec2time(toc)))

%% Load stacks

cell_bounds = dlmread(sprintf('%s/cell_bounds.txt',fov_dir));
cell_bounds = floor(cell_bounds(cell,:));

norm_stack = zeros(cell_bounds(5),cell_bounds(4),cell_bounds(6),num_channels,num_cycles);

% Load normalized stack

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/cy%02d_ch%02d.tif',cell_reg_dir,cycle,channel);
        norm_stack(:,:,:,channel,cycle) = read_3d_tif(filename,cell_bounds(5),cell_bounds(4),cell_bounds(6));
    end
end

% Load hyb + DAPI and segmentation

reg_hyb_stack = read_3d_tif(sprintf('%s/reg_hyb.tif',cell_dir),cell_bounds(5),cell_bounds(4),cell_bounds(6));
reg_dapi_stack = read_3d_tif(sprintf('%s/reg_dapi.tif',cell_dir),cell_bounds(5),cell_bounds(4),cell_bounds(6));
%bw_stack = read_3d_tif(sprintf('%s/bw.tif',cell_dir),cell_bounds(5),cell_bounds(4),cell_bounds(6));

disp(sprintf('%s: Loaded normalized stacks',sec2time(toc)))

%% Fix data with missing cycles

if condition == "scratch1"
    new_stack = zeros(cell_bounds(5),cell_bounds(4),cell_bounds(6),num_channels,19);
    new_stack(:,:,:,:,1:3) = norm_stack(:,:,:,:,1:3);
    new_stack(:,:,:,:,5:7) = norm_stack(:,:,:,:,4:6);
    new_stack(:,:,:,:,9:11) = norm_stack(:,:,:,:,7:9);
    new_stack(:,:,:,:,13:19) = norm_stack(:,:,:,:,10:16);
    norm_stack = new_stack;
    num_cycles = 19;
elseif condition == "well2" & (fov  == 0 | fov == 1)
    norm_stack(:,:,:,:,2) = 0;
end

%% Nucleus segmentation

bw_stack = zeros(size(reg_dapi_stack));

image = uint16(reg_dapi_stack);
%thresh = multithresh(image,2);
%thresh = thresh(2)-25;
thresh = prctile(image(:),85);
    
%figure;
for z=1:size(bw_stack,3)

    bw_stack(:,:,z) = reg_dapi_stack(:,:,z) > thresh;
    bw_stack(:,:,z) = imopen(bw_stack(:,:,z),strel('disk',2));
    bw_stack(:,:,z) = imdilate(bw_stack(:,:,z),strel('disk',2));
    bw_stack(:,:,z) = imfill(bw_stack(:,:,z),'holes');
    
    %imshowpair(bw_stack(:,:,z),capImage(reg_dapi_stack(:,:,z),99,'prc'),'montage')
end

labels = bwlabeln(bw_stack);
regions = regionprops(labels, 'Area');
labels_filt = zeros(size(bw_stack));

for i=1:size(regions)
    if regions(i).Area == max([regions.Area])
        labels_filt(labels==i)=i;
    end
end
      
bw_stack = bwlabeln(labels_filt);
regions = regionprops(bw_stack, 'BoundingBox');

disp(sprintf('%s: Segmented nucleus',sec2time(toc)))


%% Peak calling

all_peaks = [];
channel_peaks = {};
cycle = 1;

for channel=1:num_channels
        
    [Maxima,MaxPos,Minima,MinPos]=MinimaMaxima3D(double(norm_stack(:,:,:,channel,cycle).*bw_stack),1,0);
    
    thresh = 100;
    
    channel_peaks{channel} = MaxPos(Maxima>thresh,:);
    all_peaks = cat(1,all_peaks,channel_peaks{channel});

    disp(sprintf('%s: Found %d 3D peaks in cycle %d, channel %d',sec2time(toc),length(channel_peaks{channel}),cycle,channel));
          
end

dlmwrite(sprintf('%s/all_peaks.txt',cell_dir),all_peaks)
disp(sprintf('%s: Saved %d total 3D peaks',sec2time(toc),length(all_peaks)));

%% Peak call video

fig_dir = sprintf('figures/%s/peak_calls_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

video = VideoWriter(sprintf('%s/%s_fov%03d_cell%03d',fig_dir,condition,fov,cell),'MPEG-4');
video.FrameRate = 2.5;
open(video);

z_buffer = 2;

for z=1:size(norm_stack,3)
    
    f = figure('visible','off');
    p = tight_subplot(1,5,[0.001 0.001],[0.001 0.001],[0.001 0.001]);
    
    axes(p(1));
    imshowpair(capImage(reg_dapi_stack(:,:,z),95,'prc'),capImage(reg_hyb_stack(:,:,z),95,'prc'))
    title('DAPI')
    
    for channel=1:num_channels
    
        peaks = channel_peaks{channel};
        z_channel_peaks = peaks(ismember(peaks(:,3),(z+[-z_buffer:+z_buffer])),:);

        axes(p(channel+1));
        imshow(capImage(norm_stack(:,:,z,channel,1),95,'prc'),[]); hold on;
        plot(z_channel_peaks(:,2),z_channel_peaks(:,1),'Marker','.','MarkerEdgeColor','r','LineStyle','none');
        title(sprintf('Peak calls (n = %d)',size(peaks,1)))
        
    end
   
    %fig.Units = 'pixels';
    f.Position = [0 0 size(norm_stack,2).*2.*(5) size(norm_stack,1).*2+50];    
    writeVideo(video,getframe(gcf))
end
close(video); %close the file
close all

disp(sprintf('%s: Saved peak call video',sec2time(toc)));

%% Add background

noise_stack = norm_stack + 10;
%noise_stack(:,:,:,:,11) = [];
%num_cycles = 14;

%% Load UMIs

umi_table = readtable(umi_file,'Format','%s%d%s%d%d%s%s%d%d','FileType','text');
umis = char(cellstr(umi_table.Var6));
umi_table.Var10 = zeros(size(umi_table,1),1);

disp(sprintf('%s: Read in UMIs from %s',sec2time(toc),umi_file));

%% Matching

corr_thresh = 0.95;
spot_dim = [5 5 5];
spot_vol = [5 500];

hamming_thresh = 4;

vis_size = [5 5];

region_dim = [2 2 1];
region = zeros(region_dim(1)*2+1,region_dim(2)*2+1,region_dim(3)*2+1,num_channels,num_cycles);

purity_thresh = 6;
max_iterations = 5;

intensity_mat = zeros(size(all_peaks,1),num_channels,num_cycles+1);
norm_intensity_mat = zeros(size(all_peaks,1),num_channels,num_cycles+1);
all_spots = zeros(size(all_peaks,1),spot_dim(1)*2+1,spot_dim(2)*2+1,spot_dim(3)*2+1);

match_table = table;

spatial_locations = zeros(size(all_peaks,1),4);
match_stats = ones(size(all_peaks,1),4);
max_inds = strings(size(all_peaks,1),1);

%umi_matches = table;
%match_scores(i,:) = zeros(size(all_peaks,1),3);
%closest_matches(i) = zeros(size(all_peaks,1),1);

warning('off')

for i=1:size(all_peaks,1)
    
    peak = all_peaks(i,:);

    intensity_mat(i,:,1) = [peak 1];
    [intensity_mat(i,:,2:num_cycles+1) all_spots(i,:,:,:)] = ...
        spot_caller([peak(1) peak(2) peak(3)], noise_stack, spot_vol, spot_dim, corr_thresh);
    
    norm_intensity_mat(i,:,1) = intensity_mat(i,:,1);
    norm_intensity_mat(i,:,2:num_cycles+1) = normc(squeeze(intensity_mat(i,:,2:num_cycles+1)));
    
    [max_val max_ind] = max(squeeze(norm_intensity_mat(i,:,2:num_cycles+1).^2),[],1);
    curr_umi = char(string(squeeze(max_ind)));
    curr_prob = -log10(prod(max_val,2));

    [match_ind match_val hamming close_umis] = match_spot(curr_umi,umis,squeeze(norm_intensity_mat(i,:,:)),hamming_thresh);
    %match_info = umi_table(match_ind,:);
    %match_umi = char(match_info.Var6);

    spatial_locations(i,:) = [squeeze(norm_intensity_mat(i,:,1))];
    match_stats(i,:) = [match_ind match_val curr_prob hamming];
    max_inds(i) = strjoin(string(squeeze(max_ind)),'');
    
    if mod(i,2500) == 0
        disp(sprintf('%s: Processed %d entries',sec2time(toc),i))
    end
    
end

match_table = array2table((1:size(all_peaks,1))');
match_table{:,2:5} = spatial_locations;
match_table(:,6:15) = umi_table(match_stats(:,1),:);
match_table{:,16} = match_stats(:,2);
match_table{:,17} = max_inds(:);
match_table{:,18} = match_stats(:,3);
match_table{:,19} = match_table{:,16}-match_table{:,18};
match_table{:,20} = match_stats(:,4);

warning('on')
disp(sprintf('%s: Quantified all spots',sec2time(toc)));

%% Save initial match results

dlmwrite(sprintf('%s/intensity_mat_%s.txt',cell_dir,run),intensity_mat);
dlmwrite(sprintf('%s/all_spots_%s.txt',cell_dir,run),all_spots);
writetable(match_table,sprintf('%s/match_table_%s.txt',cell_dir,run),'QuoteStrings',true)

disp(sprintf('%s: Wrote table with %d entries to file',sec2time(toc),size(match_table,1)))

%% Filter duplicates

sort_table = sortrows(match_table,18);
[C,ia,ic] = unique(sort_table{:,11});
sort_table = sort_table(ia,:);
filter_table = sortrows(sort_table,18);

%% Add genome info

chr_order_file = 'mm10_chr_order.txt';
chr_order = readtable(sprintf('%s',chr_order_file),'FileType','text');
chr_order = table2array(chr_order);

filter_table.Var6(contains(filter_table.Var6,"Un_")) = {'unplaced'};
filter_table.Var6(contains(filter_table.Var6,"_random")) = {'unlocalized'};
filter_table.Var6(contains(filter_table.Var6,"_alt")) = {'unlocalized'};

chr_index = zeros(length(filter_table.Var6),1);
for chr=1:length(chr_order)
   chr_index(strcmp(filter_table.Var6,chr_order(chr))) = chr;
end
filter_table.Var6 = chr_index;

%% Match filtering

vis_slope = 2;
vis_int = 4;
vis_high_qual = 1.5;

radius = 1; d = datadensity(filter_table{:,16},filter_table{:,18},radius);

f = figure('visible','off');
scatter(filter_table{:,16},filter_table{:,18},10,d); hold on;
colormap(gca,jet)
axis([0 10 0 10]); hold on;

x = linspace(0,10); y = x./x*vis_high_qual;
plot(x,y); hold on;

x = linspace(0,10); y = vis_slope*x - vis_int;
plot(x,y); hold on;

bot_left = sum((filter_table{:,18}<vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18})); text(1,.5,sprintf('n=%d',bot_left),'Color','black'); hold on;
top_left = sum((filter_table{:,18}>vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18})); text(1,7,sprintf('n=%d',top_left),'Color','black'); hold on;
bot_right = sum((filter_table{:,18}<vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)>filter_table{:,18})); text(8,.5,sprintf('n=%d',bot_right),'Color','black'); hold on;
top_right = sum((filter_table{:,18}>vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)>filter_table{:,18})); text(8,7,sprintf('n=%d',top_right),'Color','black'); hold on;

disp(sprintf('%s: Mapped %d of %d (%d%%) high quality rolonies and %d total',sec2time(toc),bot_left,bot_left+bot_right,round(bot_left/(bot_left+bot_right)*100),bot_left+top_left));

title(sprintf('%d high quality, %d%% match rate',bot_left+bot_right,round(bot_left/(bot_left+bot_right)*100)))
ylabel('Spot purity score')
xlabel('UMI mapping score')

fig = gcf;
pos = get(fig,'Position');
fig.PaperPosition = [0 0 pos(3)*2 pos(4)*2];
set(fig,'PaperPositionMode','Auto','Units','Inches','PaperUnits','Inches','PaperSize',[pos(3)*2, pos(4)*2])

fig_dir = sprintf('figures/%s/match_filter_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));

sfig(2) = gca;

%% Matched reads visualization

f = figure('visible','off');

matches = (filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18};
pass_table = filter_table(matches,:);
%pass_table.Properties.VariableNames = {'index','x','y','z','region_size','umi_seq','hg19_chr','hg19_pos','hg19_dups','hg19_ind','hg38_chr','hg38_pos','hg38_dups','hg38_ind','umi_cs','match_score','umi_rol','purity_score','diff_score','hamming_dist'};
writetable(pass_table,sprintf('%s/pass_table_%s.txt',cell_dir,run),'QuoteStrings',true)

x = pass_table{:,2} * 221.87 / 2048;
y = pass_table{:,3} * 221.87 / 2048;
z = pass_table{:,4} * 0.4;

scatter3(x,y,z,15,colors(pass_table{:,6},:),'filled');view(70,75); hold on

xy_dim = max(cell_bounds(5),cell_bounds(4)) * 221.87 / 2048;
z_dim = int16(xy_dim*0.3*2048/221.87) * 221.87 / 2048;
set(gca, 'Xlim',[0 xy_dim],'YLim',[0 xy_dim],'ZLim',[0 xy_dim]); hold on;
xlabel('\mum'); ylabel('\mum'); zlabel('\mum');
title(sprintf('Cell %d from %s, FOV %d: %d reads',cell,condition,fov,length(x)))

fig = gcf;
pos = get(fig,'Position');

fig_dir = sprintf('figures/%s/matched_reads_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));

sfig(1) = gca;

%% Spot overlay xy

f = figure('visible','off');
imshow(capImage(max(max(norm_stack(:,:,:,:,1),[],3),[],4),99,'prc'),[]); hold on;

high_qual_mapped = filter_table{(filter_table{:,18}<vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18}),1};
plot(all_peaks(high_qual_mapped,2),all_peaks(high_qual_mapped,1),'rx')

low_qual_mapped = filter_table{(filter_table{:,18}>vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18}),1};
plot(all_peaks(low_qual_mapped,2),all_peaks(low_qual_mapped,1),'gx')

high_qual_unmapped = filter_table{(filter_table{:,18}<vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)>filter_table{:,18}),1};
plot(all_peaks(high_qual_unmapped,2),all_peaks(high_qual_unmapped,1),'bx','MarkerSize',3)

%low_qual_unmapped = find((filter_table{:,18}>vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)>filter_table{:,18}));
%plot(all_peaks(low_qual_unmapped,2),all_peaks(low_qual_unmapped,1),'yx')

fig = gcf;
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','Position',[0 0 pos(3)*2, pos(4)*2])

fig_dir = sprintf('figures/%s/spot_overlay_xy_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));

sfig(3) = gca;

%% Spot overlay xz

f = figure('visible','off');
imshow(capImage(squeeze(max(max(norm_stack(:,:,:,:,1),[],1),[],4))',99,'prc'),[]); hold on;

high_qual_mapped = filter_table{(filter_table{:,18}<vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18}),1};
plot(all_peaks(high_qual_mapped,2),all_peaks(high_qual_mapped,3),'rx')

low_qual_mapped = filter_table{(filter_table{:,18}>vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)<filter_table{:,18}),1};
plot(all_peaks(low_qual_mapped,2),all_peaks(low_qual_mapped,3),'gx')

high_qual_unmapped = filter_table{(filter_table{:,18}<vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)>filter_table{:,18}),1};
plot(all_peaks(high_qual_unmapped,2),all_peaks(high_qual_unmapped,3),'bx','MarkerSize',3)

%low_qual_unmapped = find((filter_table{:,18}>vis_high_qual) & ((filter_table{:,16}*vis_slope-vis_int)>filter_table{:,18}));
%plot(all_peaks(low_qual_unmapped,1),all_peaks(low_qual_unmapped,3),'yx')

fig = gcf;
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','Position',[0 0 pos(3)*2, pos(4)*2])

fig_dir = sprintf('figures/%s/spot_overlay_z_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));

sfig(4) = gca;

%% Save mapping summary

fnew = figure('visible','off');

subplot(2,2,1,copyobj(sfig(1),fnew));
subplot(2,2,2,copyobj(sfig(2),fnew));
subplot(2,2,3,copyobj(sfig(3),fnew));
subplot(2,2,4,copyobj(sfig(4),fnew));

fig = gcf;
pos = get(fig,'Position');
fig.PaperUnits = 'inches';
set(fig,'PaperPositionMode','Auto','Position',[0 0 pos(3)*2, pos(4)*2])

fig_dir = sprintf('figures/%s/cell_summary_%s',exp,run);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/%s_fov%03d_cell%03d.png',fig_dir,condition,fov,cell));
close all

disp(sprintf('%s: Saved mapping summary',sec2time(toc)));

