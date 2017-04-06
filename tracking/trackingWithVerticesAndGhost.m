%% Automatic cell tracking algorithm.

function trackingWithVerticesAndGhost(param,output_path)
%tracking Tracks cells of image sequence in time dependent manner.
%   tracking(param,output_path) Tracking starts from the first time point
%   which is segmented in independent manner. All the subsequent time
%   points are tracked and segmented based on the previous time point.
%   Input parameters are provided in param struct with following fields:
%       
%   o_filename      Input image sequence path.
%   file_number_increment       Every file_number_increment image is tracked.
%   image_interval      ['first time index' 'last time index'].
%   image_interval_start_in_original_seq    First file number index to use.
%   seq_name        Name of tracking sequence.
%   clear_boundaries        Boolean variable for storing or not storing boundary pixels of tracked cells.
%   f_img_path      Path of filtered images. If this is empty filtering is performed in this script.
%   piv_path        Path of velocity field file. If this is empty PIV is computed in this script.
%   bandpass_filter_path     Path to ImageJ script that is used for filtering image sequence before tracking.
%   imagej_path         Path to ImageJ executable.
%
%   output_path defines output folder. Tracklet structure (param) and 
%   segmented images are stored into this folder.

path2segments = 'Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testOutputTracking\segments\img_'

tic
%% Read velocity field from file if requested (otherwise initialise for computations).
if isempty(param.piv_path)
    param.velocity_field=[];
    param.vel_field_scale_factor=1;

% Read in velocities if requested. The velocity field has to cover the
% whole field of view.
else
    load(param.piv_path,'u_original','v_original','x','y');
    param.vel_field_x_coords=x{param.image_interval_start_in_original_seq};
    param.vel_field_y_coords=y{param.image_interval_start_in_original_seq};
    
    param.vel_field_scale_factor=1;
    param.velocity_field=zeros([size(param.vel_field_x_coords) 2 diff(param.image_interval)]);
    for time=(param.image_interval(1):param.image_interval(2))+param.image_interval_start_in_original_seq-1
        u_original{time}=inpaint_nans(u_original{time});
        v_original{time}=inpaint_nans(v_original{time});
        param.velocity_field(:,:,1,time-param.image_interval_start_in_original_seq+1)=...
            u_original{time}/param.vel_field_scale_factor;
        param.velocity_field(:,:,2,time-param.image_interval_start_in_original_seq+1)=...
            v_original{time}/param.vel_field_scale_factor;
    end
    clear u_original v_original x y
end

init_times=toc;tic
%% Filter first two images.
% Image size.
info=imfinfo([param.o_filename num2str(param.image_interval_start_in_original_seq,'%.04d') '.tif']);
param.img_s=[info(1).Width info(1).Height];

temp_path=[output_path 'temp_i_f.tif'];
% Filter and read in first two images if required.
if isempty(param.f_img_path)
    evaluate_imagej_script(param.bandpass_filter_path,param.imagej_path,...
        [param.o_filename  num2str((param.image_interval_start_in_original_seq-1)*...
        param.file_number_increment+1,'%.04d') '.tif'],temp_path);
    f_imgs(:,:,1)=imread(temp_path);

    evaluate_imagej_script(param.bandpass_filter_path,param.imagej_path,...
        [param.o_filename  num2str((param.image_interval_start_in_original_seq)*...
        param.file_number_increment+1,'%.04d') '.tif'],temp_path);
    f_imgs(:,:,2)=imread(temp_path);
% Read in filtered images directly.
else
	f_imgs(:,:,1)=imread([param.f_img_path num2str((param.image_interval_start_in_original_seq-1)*...
        param.file_number_increment+1,'%.04d') '.tif']);
	f_imgs(:,:,2)=imread([param.f_img_path num2str((param.image_interval_start_in_original_seq)*...
        param.file_number_increment+1,'%.04d') '.tif']);
end

init_times(end+1)=toc;tic
%% Segment first image in independent manner.

% Segment directory.
segment_path=[output_path 'segments' filesep];
if ~exist(segment_path,'dir')
    mkdir(segment_path);
end

f_img=f_imgs(:,:,1);
f_img_m=imhmin(f_img,10); %remove non important minima

init_times(end+1)=toc;tic
% Watershed
segment_temp=watershed(f_img_m);
% segment_temp(segment_temp==1) = max(segment_temp(:))


init_times(end+1)=toc;tic
%% Extract cell parameters from segmented image (first image).
s_img=uint8(~logical(segment_temp))*255;
% Collect segment information from segmented image.
info_temp=regionprops(segment_temp,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');        

init_times(end+1)=toc;tic
% Assign new cell parameters.
param.tracks=struct('t',uint16(1),'A',uint32(info_temp(1).Area),'cent',uint16(info_temp(1).Centroid),...
    'ellipse',int16([info_temp(1).MajorAxisLength info_temp(1).MinorAxisLength info_temp(1).Orientation]),...
    'perim',uint32(info_temp(1).Perimeter),'birth',int32(-1),'death',int16(0),'daughters',uint32([0 0]),...
    'neighs',cell(1),'bounds',cell(1));
param.tracks(1).neighs=cell(1);param.tracks(1).bounds=cell(1);
for cell_i=2:size(info_temp,1)
    param.tracks(cell_i)=struct('t',uint16(1),'A',uint32(info_temp(cell_i).Area),...
        'cent',uint16(info_temp(cell_i).Centroid),'ellipse',int16([info_temp(cell_i).MajorAxisLength ...
        info_temp(cell_i).MinorAxisLength info_temp(cell_i).Orientation]),...
        'perim',uint32(info_temp(cell_i).Perimeter),'birth',int32(-1),'death',int16(0),...
        'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1));
    param.tracks(cell_i).neighs=cell(1);param.tracks(cell_i).bounds=cell(1);
end

% Define the ghost cell. 1st cell in param.tracks & param.tracks.t = 0;
param.tracks=struct('t',uint16(0),'A',uint32(0),'cent',uint16(0),...
    'ellipse',int16([0 0 0]),...
    'perim',uint32(0),'birth',int32(0),'death',int16(0),'daughters',uint32([0 0]),...
    'neighs',cell(1),'bounds',cell(1));
param.tracks(1).neighs=cell(1);param.tracks(1).bounds=cell(1);

for cell_i=1:size(info_temp,1)
    param.tracks(cell_i+1)=struct('t',uint16(1),'A',uint32(info_temp(cell_i).Area),...
        'cent',uint16(info_temp(cell_i).Centroid),'ellipse',int16([info_temp(cell_i).MajorAxisLength ...
        info_temp(cell_i).MinorAxisLength info_temp(cell_i).Orientation]),...
        'perim',uint32(info_temp(cell_i).Perimeter),'birth',int32(-1),'death',int16(0),...
        'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1));
    param.tracks(cell_i+1).neighs=cell(1);
    param.tracks(cell_i+1).bounds=cell(1);
end

%% Find cell neighbours and pixels of each junction from segmenetation (first image).
edge_ind=max(segment_temp(:))+1;
segment_temp1=zeros(size(segment_temp)+[4 4]);
segment_temp1(3:end-2,3:end-2)=segment_temp;
segment_temp1([1 end],:)=edge_ind;
segment_temp1(:,[1 end])=edge_ind;

m_inds=find(segment_temp1==0);
m_points=[segment_temp1(m_inds-param.img_s(2)-4-1) segment_temp1(m_inds-param.img_s(2)-4) segment_temp1(m_inds-param.img_s(2)-4+1)...
    segment_temp1(m_inds-1)  segment_temp1(m_inds) segment_temp1(m_inds+1)...
    segment_temp1(m_inds+param.img_s(2)+4-1) segment_temp1(m_inds+param.img_s(2)+4) segment_temp1(m_inds+param.img_s(2)+4+1)];
m_sort=sort(m_points,2);
m_diff=diff(sort(m_points,2),1,2)>0;

% Search double membrane points.
ind2=find(sum(m_diff,2)==2);
[x_ind, y_ind]=find(m_diff(ind2,:)');
n_pairs=reshape(m_sort(sub2ind([size(m_diff,1) 9],ind2(y_ind),x_ind+1)),[2 size(ind2,1)])';
sort_pairs=sortrows([n_pairs [floor(m_inds(ind2)/(param.img_s(2)+4))+1-2 mod(m_inds(ind2)-1,(param.img_s(2)+4))+1-2]],[1 2]);
[uni_pairs, f_locs]=unique(sort_pairs(:,1:2),'rows','first');
[~, l_locs]=unique(sort_pairs(:,1:2),'rows','last');

% Save edge pixels and neighbours to the tracklet struct.
for e_ind=1:size(uni_pairs,1)
    % First half.
    if uni_pairs(e_ind,1)~=edge_ind
        if isempty(param.tracks(uni_pairs(e_ind,1)).bounds{1})
            param.tracks(uni_pairs(e_ind,1)).bounds{1}=cell(0);
        end
        % Mark neighbours (zero if neighbour to edge of the image).
        param.tracks(uni_pairs(e_ind,1)).neighs{1}(end+1)=uint32((uni_pairs(e_ind,2)~=edge_ind)*uni_pairs(e_ind,2));
        % Save membrane edge.
        param.tracks(uni_pairs(e_ind,1)).bounds{1}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
    end

    % Second half.
    if uni_pairs(e_ind,2)~=edge_ind
        if isempty(param.tracks(uni_pairs(e_ind,2)).bounds{1})
            param.tracks(uni_pairs(e_ind,2)).bounds{1}=cell(0);
        end
        % Mark neighbours (zero if neighbour to edge of the image).
        param.tracks(uni_pairs(e_ind,2)).neighs{1}(end+1)=uint32((uni_pairs(e_ind,1)~=edge_ind)*uni_pairs(e_ind,1));
        % Save membrane edge.
        param.tracks(uni_pairs(e_ind,2)).bounds{1}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
    end
end

% Search quadruple and triple membrane points (membrane vertices).
% Quadruple.
ind4=find(sum(m_diff,2)==4);
[x_ind, y_ind]=find(m_diff(ind4,:)');
n_quadruple=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind4(y_ind),x_ind+1)),[4 size(ind4,1)])' ...
    [floor(m_inds(ind4)/(param.img_s(2)+4))+1-2 mod(m_inds(ind4)-1,(param.img_s(2)+4))+1-2]];

% Triple.
ind3=find(sum(m_diff,2)==3);
[x_ind, y_ind]=find(m_diff(ind3,:)');
n_triples=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind3(y_ind),x_ind+1)),[3 size(ind3,1)])' ...
    [floor(m_inds(ind3)/(param.img_s(2)+4))+1-2 mod(m_inds(ind3)-1,(param.img_s(2)+4))+1-2]];

% Combine quadruple points with the triple points.
n_triples=[n_triples;n_quadruple(:,[1 2 3 5 6]);n_quadruple(:,[1 2 4 5 6]);n_quadruple(:,[1 3 4 5 6]);];

for row_ind=1:size(n_triples,1)
    for col_ind=2:3
        % First half.
        if n_triples(row_ind,1)~=edge_ind
            n_ind=find(param.tracks(n_triples(row_ind,1)).neighs{end}==((n_triples(row_ind,col_ind)~=edge_ind)*n_triples(row_ind,col_ind)));
            if isempty(n_ind)
                param.tracks(n_triples(row_ind,1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(n_triples(row_ind,1)).neighs{end}(end+1)=uint32(n_triples(row_ind,col_ind));
            else
                param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1,1);
                param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,1)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
        % Second half.
        if n_triples(row_ind,col_ind)~=edge_ind
            n_ind=find(param.tracks(n_triples(row_ind,col_ind)).neighs{end}==((n_triples(row_ind,1)~=edge_ind)*n_triples(row_ind,1)));
            if isempty(n_ind)
                param.tracks(n_triples(row_ind,col_ind)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(n_triples(row_ind,col_ind)).neighs{end}(end+1)=uint32(n_triples(row_ind,1));
            else
                param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1,1);
                param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,col_ind)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
    end
end

for row_ind=1:size(n_triples,1)
    % First half.
    if n_triples(row_ind,2)~=edge_ind
        n_ind=find(param.tracks(n_triples(row_ind,2)).neighs{end}==((n_triples(row_ind,3)~=edge_ind)*n_triples(row_ind,3)));
        if isempty(n_ind)
            param.tracks(n_triples(row_ind,2)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
            param.tracks(n_triples(row_ind,2)).neighs{end}(end+1)=uint32(n_triples(row_ind,3));
        else
            param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1,1)+1;
            ind2add=param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1,1);
            param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,2)).bounds{end}{n_ind}((ind2add+1):end,:)];
        end
    end
    % Second half.
    if n_triples(row_ind,3)~=edge_ind
        n_ind=find(param.tracks(n_triples(row_ind,3)).neighs{end}==((n_triples(row_ind,2)~=edge_ind)*n_triples(row_ind,2)));
        if isempty(n_ind)
            param.tracks(n_triples(row_ind,3)).bounds{end}{end+1}=uint16([1 0 ; n_triples(row_ind,4:5)]);
            param.tracks(n_triples(row_ind,3)).neighs{end}(end+1)=uint32(n_triples(row_ind,2));
        else
            param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1,1)=param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1,1)+1;
            ind2add=param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1,1);
            param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}=[param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(n_triples(row_ind,3)).bounds{end}{n_ind}((ind2add+1):end,:)];
        end
    end
end

% 2x2 pixel squares are problem to the endpoint detection, thus these
% points must be marked as separate end points.
square_inds=find(conv2(single(~logical(segment_temp1)),[1 1;1 1],'same')==4);
square_inds=[square_inds;square_inds+1;square_inds+param.img_s(2)+4;square_inds+param.img_s(2)+4+1];

m_points=[segment_temp1(square_inds-param.img_s(2)-4-1) segment_temp1(square_inds-param.img_s(2)-4) segment_temp1(square_inds-param.img_s(2)-4+1)...
    segment_temp1(square_inds-1)  segment_temp1(square_inds) segment_temp1(square_inds+1)...
    segment_temp1(square_inds+param.img_s(2)+4-1) segment_temp1(square_inds+param.img_s(2)+4) segment_temp1(square_inds+param.img_s(2)+4+1)];
m_sort=sort(m_points,2);
m_diff=diff(sort(m_points,2),1,2)>0;

[x_ind, y_ind]=find(m_diff');
n_pairs=[reshape(m_sort(sub2ind([size(m_diff,1) 9],y_ind,x_ind+1)),[2 size(m_sort,1)])' [floor(square_inds/(param.img_s(2)+4))+1-2 mod(square_inds-1,(param.img_s(2)+4))+1-2]];

% Save edge pixels and neighbours to the tracklet struct.
for e_ind=1:size(n_pairs,1)
    % first half
    n_ind=find(param.tracks(n_pairs(e_ind,1)).neighs{end}==n_pairs(e_ind,2));
    param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1,1)=param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1,1)+1;
    ind2add=param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1,1);
    param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}=[param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(n_pairs(e_ind,1)).bounds{end}{n_ind}((ind2add+1):end,:)];

    % second half
    n_ind=find(param.tracks(n_pairs(e_ind,2)).neighs{end}==n_pairs(e_ind,1));
    param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1,1)=param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1,1)+1;
    ind2add=param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1,1);
    param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}=[param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(n_pairs(e_ind,2)).bounds{end}{n_ind}((ind2add+1):end,:)];
end

% Write segmentation image.
imwrite(s_img,[segment_path param.seq_name '_' num2str(1,'%.03d') '.tif'],'tif','compression','none');

% Print running times of first time step.
disp(['init times ' num2str([init_times sum(init_times)])]);
init_times(end+1)=toc;
tic
param  =  produceVertices(1,param,path2segments);
%% Track and segment all the subsequent time points in time dependent manner.
time_slots=[];
for time=1:(param.image_interval(2)-1)

    disp(time)
    %% (time loop) Compute velocity if not preloaded.
    tic;slot_i=0;

    % Check if already computed and loaded. The loaded velocity field has
    % to be padded from edges to cover the whole field of view.
    if isempty(param.piv_path)
        settings=[];
        settings.preset=1;
        [vel_field_x_coords,vel_field_y_coords,u_filtered,v_filtered,~]=...
            PIVlab_vel_field_settings(imresize(f_imgs(:,:,1),param.vel_field_scale_factor),...
            imresize(f_imgs(:,:,2),param.vel_field_scale_factor),settings);

        param.vel_field_x_coords=vel_field_x_coords/param.vel_field_scale_factor;
        param.vel_field_y_coords=vel_field_y_coords/param.vel_field_scale_factor;
        param.velocity_field(:,:,1,time)=u_filtered/param.vel_field_scale_factor;
        param.velocity_field(:,:,2,time)=v_filtered/param.vel_field_scale_factor;    
    end

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic

    %% (time loop) Clear previous cell junction coordinates from memory if requested.
    if param.clear_boundaries
        for i_ind=1:length(param.tracks)
            for j_ind=1:(length(param.tracks(i_ind).bounds)-2)
                param.tracks(i_ind).bounds{j_ind}=[];
            end
        end
    end
    
    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic    
    
    %% (time loop) Terminate ingressing cells. Area threshold or less than 3 neighbours.
    inds_alive=zeros(length(param.tracks),3);
    inds_alive_ind=1;
    for cell_i=1:length(param.tracks)
        cell_i
        if param.tracks(cell_i).death==0
            % Is area big enough? Is number of neighbours more than 2?
            if param.tracks(cell_i).A(end)>=15 & length(param.tracks(cell_i).bounds{end})>2
                inds_alive(inds_alive_ind,1:3)=[cell_i 0 0];
                inds_alive_ind=inds_alive_ind+1;
            % If does not have enough neighbours kill the tracklet.
            elseif param.tracks(cell_i).A(end)>=15
                param.tracks(cell_i).death=-7;
            % If area is not big enough kill the tracklet.
            else
                param.tracks(cell_i).death=-1;
            end
        end
    end
    inds_alive(inds_alive_ind:end,:)=[];
    n_a_segs=size(inds_alive,1);

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic
    
    %% (time loop) Find inflow regions and add seed points.

    % Bottom.
    boundary=[(1:size(f_imgs,2))' ones(size(f_imgs,2),1)*size(f_imgs,1)];
    boundary(:,3)=boundary(:,1)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,1,time),boundary(:,1),boundary(:,2));
    boundary(:,4)=boundary(:,2)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,2,time),boundary(:,1),boundary(:,2));        

    inner_bound=floor(min(boundary(:,4)));

    f_img_m_indpendent=imhmin(f_imgs(inner_bound:end,:,2),10);
    segment_temp_independent=watershed(f_img_m_indpendent);
    info_temp_independent=regionprops(segment_temp_independent,'centroid');
    independent_cents = round(reshape(struct2array(info_temp_independent),[2 length(info_temp_independent)])');
    independent_cents(:,2)=independent_cents(:,2)+inner_bound-1;
    
    independent_cents((independent_cents(:,2))<=boundary(independent_cents(:,1),4),:)=[];
    independent_cents_to_collect=independent_cents;

    % Top.
    boundary=[(1:size(f_imgs,2))' ones(size(f_imgs,2),1)*1];
    boundary(:,3)=boundary(:,1)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,1,time),boundary(:,1),boundary(:,2));
    boundary(:,4)=boundary(:,2)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,2,time),boundary(:,1),boundary(:,2));        

    inner_bound=ceil(max(boundary(:,4)));

    f_img_m_indpendent=imhmin(f_imgs(1:inner_bound,:,2),10);
    segment_temp_independent=watershed(f_img_m_indpendent);
    info_temp_independent=regionprops(segment_temp_independent,'centroid');
    independent_cents = round(reshape(struct2array(info_temp_independent),[2 length(info_temp_independent)])');

    independent_cents((independent_cents(:,2))>=boundary(independent_cents(:,1),4),:)=[];
    independent_cents_to_collect=[independent_cents_to_collect;independent_cents;];

    % Left.
    boundary=[ones(size(f_imgs,1),1)*1 (1:size(f_imgs,1))'];
    boundary(:,3)=boundary(:,1)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,1,time),boundary(:,1),boundary(:,2));
    boundary(:,4)=boundary(:,2)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,2,time),boundary(:,1),boundary(:,2));        

    inner_bound=ceil(max(boundary(:,3)));

    f_img_m_indpendent=imhmin(f_imgs(:,1:inner_bound,2),10);
    segment_temp_independent=watershed(f_img_m_indpendent);
    info_temp_independent=regionprops(segment_temp_independent,'centroid');
    independent_cents = round(reshape(struct2array(info_temp_independent),[2 length(info_temp_independent)])');

    independent_cents((independent_cents(:,1))>=boundary(independent_cents(:,2),3),:)=[];
    independent_cents_to_collect=[independent_cents_to_collect;independent_cents;];

    % Right.
    boundary=[ones(size(f_imgs,1),1)*size(f_imgs,2) (1:size(f_imgs,1))'];
    boundary(:,3)=boundary(:,1)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,1,time),boundary(:,1),boundary(:,2));
    boundary(:,4)=boundary(:,2)+interp2(param.vel_field_x_coords,param.vel_field_y_coords,param.velocity_field(:,:,2,time),boundary(:,1),boundary(:,2));        

    inner_bound=floor(min(boundary(:,3)));

    f_img_m_indpendent=imhmin(f_imgs(:,inner_bound:end,2),10);
    segment_temp_independent=watershed(f_img_m_indpendent);
    info_temp_independent=regionprops(segment_temp_independent,'centroid');
    independent_cents = round(reshape(struct2array(info_temp_independent),[2 length(info_temp_independent)])');
    independent_cents(:,1)=independent_cents(:,1)+inner_bound-1;

    independent_cents(independent_cents(:,1)<=boundary(independent_cents(:,2),3),:)=[];
    independent_cents_to_collect=[independent_cents_to_collect;independent_cents;];            

    % Collect all cells that are to be tracked.
    inds_alive=[inds_alive; [zeros(size(independent_cents_to_collect,1),1) round(independent_cents_to_collect)]];

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic
    
    % Make new segment mask.
    % Estimate locations of live tracks in frame t+1 using the velocity
    % field.
    segment_mask=logical(zeros(fliplr(param.img_s)));
    s_factor=0.4;
    % Collect all the points that are tracked further.
    points=zeros(n_a_segs,2);
    for a_ind=1:n_a_segs
        points(a_ind,:)=double(param.tracks(inds_alive(a_ind,1)).cent(end,:));
    end

    %% (time loop) Use PIV to estimate locations of cells in next time point.
    points=points+round([interp2(param.vel_field_x_coords,param.vel_field_y_coords,...
        param.velocity_field(:,:,1,time),points(:,1),points(:,2),'linear') ...
        interp2(param.vel_field_x_coords,param.vel_field_y_coords,...
        param.velocity_field(:,:,2,time),points(:,1),points(:,2),'linear')]);    
    % Add points from inflowed boundary.
    points = [points; inds_alive((n_a_segs+1:end),2:3)];
    
    points(:,1)=max(points(:,1), 2);
    points(:,1)=min(points(:,1), param.img_s(1)-1);
    points(:,2)=max(points(:,2), 2);
    points(:,2)=min(points(:,2), param.img_s(2)-1);

    %% (time loop) Create segmentation seed points into locations.
    % Create segmentation mask point by point.
    for a_ind=1:size(inds_alive,1)
        cell_i=inds_alive(a_ind,1);
        point=points(a_ind,:);

        if a_ind<=n_a_segs
            cell_cent=double(param.tracks(cell_i).cent(end,1));
            
            e_points=param.tracks(cell_i).bounds{end};
            for e_ind=1:length(e_points)
                e_points{e_ind}(1,:)=[];
            end
            e_points=double(cell2mat(e_points'));

            e_points=[e_points(:,1)+(-double(param.tracks(cell_i).cent(end,1))+point(1))...
                e_points(:,2)+(-double(param.tracks(cell_i).cent(end,2))+point(2))];                                    
        else
            cell_cent=point;
            e_points=[-1 -1 -1 0 1 1 1 0;-1 0 1 1 1 0 -1 -1]';
            e_points(:,1)=e_points(:,1)+point(1);
            e_points(:,2)=e_points(:,2)+point(2);
        end
        
        e_points(e_points(:,1)<1,1)=1;
        e_points(e_points(:,2)<1,2)=1;
        e_points(e_points(:,1)>param.img_s(1),1)=param.img_s(1);
        e_points(e_points(:,2)>param.img_s(2),2)=param.img_s(2);
        
        e_points(:,1)=round((e_points(:,1)-point(1))*s_factor+double(point(1)));
        e_points(:,2)=round((e_points(:,2)-point(2))*s_factor+double(point(2)));
        
        segment_mask(e_points(:,2)+(e_points(:,1)-1)*param.img_s(2))=1;        
        inds_alive(a_ind,2:3)=e_points(1,:);
    end
    
    segment_mask=imfill(segment_mask,'holes');
    l_img=bwlabel(segment_mask);    
    
    already_found=zeros(size(inds_alive,1),1);
    
    overlapping_blobs=zeros(10000,2);
    overlapping_blobs_ind=1;
    for a_ind=1:size(inds_alive,1)
        if already_found(l_img(inds_alive(a_ind,3),inds_alive(a_ind,2)))==1            
            overlapping_blobs(overlapping_blobs_ind,:)=[inds_alive(a_ind,1) a_ind];
            overlapping_blobs_ind=overlapping_blobs_ind+1;
        else            
            already_found(l_img(inds_alive(a_ind,3),inds_alive(a_ind,2)))=1;
        end
    end
    overlapping_blobs((overlapping_blobs_ind):end,:)=[];    
    
    if ~isempty(overlapping_blobs)
        % Terminate overlapping tracks.
        for term_ind=1:size(overlapping_blobs,1)
            % Check that overlap did not occur in inflow boundary segment.
            if overlapping_blobs(term_ind,1)>0
                % Set death=-6 for cells where overlapping occured.
                param.tracks(overlapping_blobs(term_ind,1)).death=-6;
            end
        end
        inds_alive(overlapping_blobs(:,2),:)=[];
    end

    % Apply mask.
    f_img=f_imgs(:,:,2);
    f_img_m=f_img/2+2^7+1;
    mod_f_img_m=f_img_m.*uint8(~segment_mask);

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic
    
    %% (time loop) Segment using seed point mask.
    mod_f_img_m=imhmin(mod_f_img_m,2^7);
    
    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic

    % Watershed segment frame t+1 using seeded image.
    segment_temp=watershed(mod_f_img_m);
    
    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic    
    
    %% (time loop) Extract cell parameters from segmented image.
    s_img(:,:,2)=uint8(~logical(segment_temp))*255;

    % Collect segment information from segmented image.
    info_temp=regionprops(segment_temp,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic    
    
    % Assign new cell parameters to existing struct.
    number_of_tracked_cells=length(param.tracks);

    for cell_i=1:size(inds_alive,1)
        % Cell that is being tracked already.
        if inds_alive(cell_i,1)>0
            info_ind=segment_temp(inds_alive(cell_i,3),inds_alive(cell_i,2));
            inds_alive(cell_i,4)=info_ind;
            param.tracks(inds_alive(cell_i,1)).t(end+1)=time+1;
            param.tracks(inds_alive(cell_i,1)).A(end+1)=uint32(info_temp(info_ind).Area);
            param.tracks(inds_alive(cell_i,1)).cent(end+1,1:2)=uint16(info_temp(info_ind).Centroid);
            param.tracks(inds_alive(cell_i,1)).ellipse(end+1,1:3)=int16([info_temp(info_ind).MajorAxisLength info_temp(info_ind).MinorAxisLength info_temp(info_ind).Orientation]);
            param.tracks(inds_alive(cell_i,1)).perim(end+1)=uint32(info_temp(info_ind).Perimeter);
            param.tracks(inds_alive(cell_i,1)).bounds{end+1}=cell(0);
            param.tracks(inds_alive(cell_i,1)).neighs{end+1}=zeros(0);
        % Cell that appears due to a positive inflow on boundary.
        else
            info_ind=segment_temp(inds_alive(cell_i,3),inds_alive(cell_i,2));
            inds_alive(cell_i,4)=info_ind;
            param.tracks(end+1).t(1)=time+1;
            param.tracks(end).A(1)=uint32(info_temp(info_ind).Area);
            param.tracks(end).cent(1,1:2)=uint16(info_temp(info_ind).Centroid);
            param.tracks(end).ellipse(1,1:3)=int16([info_temp(info_ind).MajorAxisLength info_temp(info_ind).MinorAxisLength info_temp(info_ind).Orientation]);
            param.tracks(end).perim(1)=uint32(info_temp(info_ind).Perimeter);
            param.tracks(end).bounds{1}=cell(0);
            param.tracks(end).neighs{1}=zeros(0);

            param.tracks(end).birth=-42;
            param.tracks(end).death=0;
            param.tracks(end).daughters=[0 0];

            inds_alive(cell_i,1)=length(param.tracks);
        end
    end

    %% (time loop) Find cell neighbours and pixels of each junction from segmenetation.
    edge_ind=max(segment_temp(:))+1;
    segment_temp1=zeros(size(segment_temp)+[4 4]);
    segment_temp1(3:end-2,3:end-2)=segment_temp;
    segment_temp1([1 end],:)=edge_ind;
    segment_temp1(:,[1 end])=edge_ind;

    m_inds=find(segment_temp1==0);
    m_points=[segment_temp1(m_inds-param.img_s(2)-4-1) segment_temp1(m_inds-param.img_s(2)-4) segment_temp1(m_inds-param.img_s(2)-4+1)...
        segment_temp1(m_inds-1)  segment_temp1(m_inds) segment_temp1(m_inds+1)...
        segment_temp1(m_inds+param.img_s(2)+4-1) segment_temp1(m_inds+param.img_s(2)+4) segment_temp1(m_inds+param.img_s(2)+4+1)];
    m_sort=sort(m_points,2);
    m_diff=diff(sort(m_points,2),1,2)>0;

    % Search double membrane points.
    ind2=find(sum(m_diff,2)==2);
    [x_ind, y_ind]=find(m_diff(ind2,:)');
    n_pairs=reshape(m_sort(sub2ind([size(m_diff,1) 9],ind2(y_ind),x_ind+1)),[2 size(ind2,1)])';
    sort_pairs=sortrows([n_pairs [floor(m_inds(ind2)/(param.img_s(2)+4))+1-2 mod(m_inds(ind2)-1,(param.img_s(2)+4))+1-2]],[1 2]);
    [uni_pairs, f_locs]=unique(sort_pairs(:,1:2),'rows','first');
    [~, l_locs]=unique(sort_pairs(:,1:2),'rows','last');

    inds_alive=sortrows(inds_alive,4);
    inds_alive(end+1,1)=0;

    % Save edge pixels and neighbours to the tracklet struct.
    for e_ind=1:size(uni_pairs,1)
        % First half.
        if uni_pairs(e_ind,1)~=edge_ind
            % Mark neighbours (zero if neighbour to edge of the image).
            param.tracks(inds_alive(uni_pairs(e_ind,1),1)).neighs{end}(end+1)=uint32(inds_alive(uni_pairs(e_ind,2),1));
            % Save membrane edge.
            param.tracks(inds_alive(uni_pairs(e_ind,1),1)).bounds{end}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
        end

        % Second half.
        if uni_pairs(e_ind,2)~=edge_ind
            % Mark neighbours (zero if neighbour to edge of the image).
            param.tracks(inds_alive(uni_pairs(e_ind,2),1)).neighs{end}(end+1)=uint32(inds_alive(uni_pairs(e_ind,1),1));
            % Save membrane edge.
            param.tracks(inds_alive(uni_pairs(e_ind,2),1)).bounds{end}{end+1}=uint16([0 0;sort_pairs(f_locs(e_ind):l_locs(e_ind),3:4)]);
        end
    end

    % Search quadruple and triple membrane points (membrane vertices).
    % Quadruple.
    ind4=find(sum(m_diff,2)==4);
    [x_ind, y_ind]=find(m_diff(ind4,:)');
    n_quadruple=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind4(y_ind),x_ind+1)),[4 size(ind4,1)])' [floor(m_inds(ind4)/(param.img_s(2)+4))+1-2 mod(m_inds(ind4)-1,(param.img_s(2)+4))+1-2]];

    % Triple.
    ind3=find(sum(m_diff,2)==3);
    [x_ind, y_ind]=find(m_diff(ind3,:)');
    n_triples=[reshape(m_sort(sub2ind([size(m_diff,1) 9],ind3(y_ind),x_ind+1)),[3 size(ind3,1)])' [floor(m_inds(ind3)/(param.img_s(2)+4))+1-2 mod(m_inds(ind3)-1,(param.img_s(2)+4))+1-2]];

    % Combine quadruple points to the triple points.
    n_triples=[n_triples;n_quadruple(:,[1 2 3 5 6]);n_quadruple(:,[1 2 4 5 6]);n_quadruple(:,[1 3 4 5 6]);];

    for row_ind=1:size(n_triples,1)
        for col_ind=2:3
            % First half.
            if n_triples(row_ind,1)~=edge_ind
                n_ind=find(param.tracks(inds_alive(n_triples(row_ind,1),1)).neighs{end}==inds_alive(n_triples(row_ind,col_ind),1));
                if isempty(n_ind)
                    param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                    param.tracks(inds_alive(n_triples(row_ind,1),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,col_ind),1));
                else
                    param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1,1)+1;
                    ind2add=param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1,1);
                    param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,1),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
                end
            end
            % Second half.
            if n_triples(row_ind,col_ind)~=edge_ind
                n_ind=find(param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).neighs{end}==inds_alive(n_triples(row_ind,1),1));
                if isempty(n_ind)
                    param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                    param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,1),1));
                else
                    param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1,1)+1;
                    ind2add=param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1,1);
                    param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,col_ind),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
                end
            end
        end
    end

    for row_ind=1:size(n_triples,1)
        % First half.
        if n_triples(row_ind,2)~=edge_ind
            n_ind=find(param.tracks(inds_alive(n_triples(row_ind,2),1)).neighs{end}==inds_alive(n_triples(row_ind,3),1));
            if isempty(n_ind)
                param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(inds_alive(n_triples(row_ind,2),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,3),1));
            else
                param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1,1);
                param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,2),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
        % Second half.
        if n_triples(row_ind,3)~=edge_ind
            n_ind=find(param.tracks(inds_alive(n_triples(row_ind,3),1)).neighs{end}==inds_alive(n_triples(row_ind,2),1));
            if isempty(n_ind)
                param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{end+1}=uint16([1 0;n_triples(row_ind,4:5)]);
                param.tracks(inds_alive(n_triples(row_ind,3),1)).neighs{end}(end+1)=uint32(inds_alive(n_triples(row_ind,2),1));
            else
                param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1,1)+1;
                ind2add=param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1,1);
                param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_triples(row_ind,4:5)) ; param.tracks(inds_alive(n_triples(row_ind,3),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
            end
        end
    end
    
    % 2x2 pixel squares are problem to the endpoint detection, thus these
    % points must be marked as separate end points.
    square_inds=find(conv2(single(~logical(segment_temp1)),[1 1;1 1],'same')==4);
    square_inds=[square_inds;square_inds+1;square_inds+param.img_s(2)+4;square_inds+param.img_s(2)+4+1];

    m_points=[segment_temp1(square_inds-param.img_s(2)-4-1) segment_temp1(square_inds-param.img_s(2)-4) segment_temp1(square_inds-param.img_s(2)-4+1)...
        segment_temp1(square_inds-1)  segment_temp1(square_inds) segment_temp1(square_inds+1)...
        segment_temp1(square_inds+param.img_s(2)+4-1) segment_temp1(square_inds+param.img_s(2)+4) segment_temp1(square_inds+param.img_s(2)+4+1)];
    m_sort=sort(m_points,2);
    m_diff=diff(sort(m_points,2),1,2)>0;

    [x_ind, y_ind]=find(m_diff');
    n_pairs=[reshape(m_sort(sub2ind([size(m_diff,1) 9],y_ind,x_ind+1)),[2 size(m_sort,1)])' [floor(square_inds/(param.img_s(2)+4))+1-2 mod(square_inds-1,(param.img_s(2)+4))+1-2]];

    % Save edge pixels and neighbours to the tracklet struct.
    for e_ind=1:size(n_pairs,1)
        % First half.
        n_ind=find(param.tracks(inds_alive(n_pairs(e_ind,1),1)).neighs{end}==inds_alive(n_pairs(e_ind,2),1));
        param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1,1)+1;
        ind2add=param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1,1);
        param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(inds_alive(n_pairs(e_ind,1),1)).bounds{end}{n_ind}((ind2add+1):end,:)];

        % Second half.
        n_ind=find(param.tracks(inds_alive(n_pairs(e_ind,2),1)).neighs{end}==inds_alive(n_pairs(e_ind,1),1));
        param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1,1)=param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1,1)+1;
        ind2add=param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1,1);
        param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}=[param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}(1:ind2add,:) ; uint16(n_pairs(e_ind,3:4)) ; param.tracks(inds_alive(n_pairs(e_ind,2),1)).bounds{end}{n_ind}((ind2add+1):end,:)];
    end

    inds_alive(end,:)=[];

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic

    %% (time loop) Terminate cell with less than 3 neighbours.
    % Validate if sufficient number of neighbours. Terminate tracks with 
    % too few (<3) neighbours.

    rm_inds_alive=[];
    for track_i=1:size(inds_alive,1)
        if length(param.tracks(inds_alive(track_i,1)).neighs{end})<3
            param.tracks(inds_alive(track_i,1)).death=-7;
            rm_inds_alive(end+1)=track_i;
        end
    end

    inds_alive(rm_inds_alive,:)=[];

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic
    
    %% (time loop) Detect cell divisions (see details below).    
    % 
    % # Consider each cell individually
    % # Initial condition for division: Cell exists in first time point,
    % cell flowed in or cell track is more than 5 h (100 time points) long.
    % # Secondary conditions for division: aspect ratio of cell (length 
    % along major axis divided by length along minor axis) is greater than
    % 2.5 and area of cell meets following conditions 250 < A < 2500 (detect dumbbells).
    % # If above conditions are met cell is segmented into two by adding
    % two new seeds onto the major axis and segment.
    % # If one segment is 50 % larger than the other 
    % division is not induced.
    % # If division is detected then check eccentricities of both new
    % segments if smaller than 0.7 division is induced.
    % # Alternatively if average pixel intensity on the newly created 
    % segment boundary is greater than average pixel intensity on all the 
    % other junctions division is induced.
    
    % Detect division and erroneous segmentation.
    % Loop through all alive cells.
    for track_i=1:size(inds_alive,1)
        division_ok=0;
        split_ok=1;

        % Decide if cell might be dividing.
        cell_i=inds_alive(track_i);

        % Skip inflown cells.
        if cell_i>number_of_tracked_cells
            continue;
        end

        % Skip cells on boundary of field of view.
        zero_neighs=0;
        for b_i=1:length(param.tracks(cell_i).bounds{end})
            if param.tracks(cell_i).bounds{end}{b_i}(1,1)==0
                zero_neighs=1;
            end
        end
        if zero_neighs
            continue;
        end

        % Check if cell is old enough to divide.
        if ((param.tracks(cell_i).t(end)-param.tracks(cell_i).t(1)+1)>100)
            division_ok=1;
        end
        
        % Cell that exists in first frame or was created on inflow boundary may divide.      
        if param.tracks(cell_i).birth==-1 | param.tracks(cell_i).birth==-42
            division_ok=1;
        end
        
        % An at least three time points long split cell may divide.
        if param.tracks(cell_i).birth>0 & param.tracks(param.tracks(cell_i).birth).death==-8 & ...
                ((param.tracks(cell_i).t(end)-param.tracks(cell_i).t(1)+1)>3)
            division_ok=1;
        end
        
        % A cell with less than two neighbours may not divide.
        if size(param.tracks(cell_i).t,2)<2
            division_ok=0;
        end

        % Compute aspect ratio of dividing cell (= length along major axis 
        % divided by length along minor axis).
        % If the ratio is greater than 2.5 the cell may divide. Following
        % area condition must also be met 250 < A < 2500.
        if division_ok
            t_ind=find(param.tracks(cell_i).t==(time+1));
            if ~isempty(t_ind)
                bound_points=[];

                for i_ind=1:length(param.tracks(cell_i).bounds{t_ind})
                    bound_points=[bound_points; param.tracks(cell_i).bounds{t_ind}{i_ind}(2:end,:)];
                end
                max_points=double(max(bound_points));
                min_points=double(min(bound_points));

                max_ray=max(ceil(sqrt(([min_points(1) min_points(1) max_points(1) max_points(1)]-double(param.tracks(cell_i).cent(t_ind,1))).^2+...
                    ([min_points(2) max_points(2) min_points(2) max_points(2)]-double(param.tracks(cell_i).cent(t_ind,2))).^2)));

                line_points1=uint16([(0:max_ray)'*cosd(-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,1)) (0:max_ray)'*sind(-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,2))]);
                i1_points=intersect(bound_points,line_points1,'rows');

                line_points2=uint16([(0:max_ray)'*cosd(180-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,1)) (0:max_ray)'*sind(180-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,2))]);
                i2_points=intersect(bound_points,line_points2,'rows');

                line_points1=uint16([(0:max_ray)'*cosd(90-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,1)) (0:max_ray)'*sind(90-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,2))]);
                i3_points=intersect(bound_points,line_points1,'rows');

                line_points2=uint16([(0:max_ray)'*cosd(-90-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,1)) (0:max_ray)'*sind(-90-double(param.tracks(cell_i).ellipse(t_ind,3)))+double(param.tracks(cell_i).cent(t_ind,2))]);
                i4_points=intersect(bound_points,line_points2,'rows');

                if ~isempty(i1_points) & ~isempty(i2_points) & ~isempty(i3_points) & ~isempty(i4_points)
                    MM(1)=sqrt(sum((double(i1_points(1,:))-double(i2_points(1,:))).^2));                    

                    MM(2)=sqrt(sum((double(i3_points(1,:))-double(i4_points(1,:))).^2));

                    if MM(1)/MM(2)>2.5 & (param.tracks(cell_i).A(end)>= (1000/4)) & (param.tracks(cell_i).A(end)<=(10000/4))
                        division_ok=1;
                    else
                        division_ok=0;
                    end
                else
                    division_ok=0;                                        
                end
            else
                division_ok=0;                                    
            end
        end

        % Detect errors in segmentation, i.e. segment need to be split into two.
        % Check whether cell should split.
        score=[];
        t_split=time-double(param.tracks(cell_i).t(1))+1;

        % Long enough track (>3 time points).
        if time>1 & ((param.tracks(cell_i).t(end)-param.tracks(cell_i).t(1)+1)>3)

            e_points=[];
            for b_i=1:length(param.tracks(cell_i).bounds{t_split})
                e_points=[e_points; param.tracks(cell_i).bounds{t_split}{b_i}(2:end,:)];
            end
            e_points=double(e_points);

            e_points1=[];
            for b_i=1:length(param.tracks(cell_i).bounds{t_split+1})
                e_points1=[e_points1; param.tracks(cell_i).bounds{t_split+1}{b_i}(2:end,:)];
            end
            e_points1=double(e_points1);

            % Find median displacement of neighbouring cells.
            v_m=[0 0];
            for n_i=1:length(param.tracks(cell_i).neighs{t_split})
                n_cell_i=param.tracks(cell_i).neighs{t_split}(n_i);
                if n_cell_i==0
                    continue
                end

                t_n=time-param.tracks(n_cell_i).t(1)+1;

                if t_n>(length(param.tracks(n_cell_i).t)-1)
                    continue
                end

                v_m(end+1,:)=diff(double(param.tracks(n_cell_i).cent(t_n:t_n+1,:)));
            end
            v_m=round(median(v_m,1));

            % Overlap old segment with new one.
            e_points1=e_points1-repmat(v_m,[length(e_points1) 1]);

            d_min=min([e_points ; e_points1]);
            d_max=max([e_points ; e_points1]);

            dims=d_max-d_min+1;

            fill_img=zeros(fliplr(dims));
            fill_img(sub2ind(fliplr(dims),e_points(:,2)-d_min(2)+1,e_points(:,1)-d_min(1)+1))=1;
            fill_img=imfill(fill_img,'holes');

            fill_img1=zeros(fliplr(dims));
            fill_img1(sub2ind(fliplr(dims),e_points1(:,2)-d_min(2)+1,e_points1(:,1)-d_min(1)+1))=1;
            fill_img1=imfill(fill_img1,'holes');
            
            % Compute how similar the segments are: ("sum of areas not
            % overlapping")/("sum of areas overlapping)").
            score=sum(sum((fill_img(:)+fill_img1(:))==1))/sum(sum((fill_img(:)+fill_img1(:))==2));
            
            % Skip the split if the score is smaller than 0.5 or if segment
            % is decreasing in size or if new area is smaller than 250.
            if score<0.5 | sum(fill_img1(:))<sum(fill_img(:)) | (sum(fill_img1(:))<(2000/4))
                split_ok=0;
            end
        else
            split_ok=0;
        end

        % If division seems to occur test it, use watershed algorithm to divide the mother
        % candidate segment into two.
        if division_ok
            point=param.tracks(cell_i).cent(end,:);
            
            % Find points of shrinken segments.
            e_points=double(cell2mat(param.tracks(cell_i).bounds{end}'));
            % Remove lengths of boundary points.
            e_points_l=1;
            for b_ind=1:(length(param.tracks(cell_i).bounds{end})-1)
                e_points_l(end+1)=size(param.tracks(cell_i).bounds{end}{b_ind},1);
            end
            % Create image for the watershedder.
            %   - Background is zeros.
            %   - Segment is outlined with boundary values 2^8-1.
            %   - Two seed points are marking locations new daughter cells.

            e_points(cumsum(e_points_l),:)=[];
            min_e_points=min(e_points);
            max_e_points=max(e_points);
            div_w_s=fliplr(max_e_points-min_e_points+1);

            mask=logical(zeros(div_w_s+[2 2]));
            mask(sub2ind(div_w_s+[2 2],e_points(:,2)-min_e_points(2)+2,e_points(:,1)-min_e_points(1)+2))=1;
            mask1=imfill(mask,'holes');
            mask2=mask1-mask;
            [y_c, x_c]=find(mask2);

            div_img=zeros(div_w_s+[2 2]);
            div_img(sub2ind(div_w_s+[2 2],y_c,x_c))=f_img_m(sub2ind(fliplr(param.img_s),y_c-2+min_e_points(2),x_c-2+min_e_points(1)));
            div_img(div_img==(2^8-1))=2^8-2;

            [y_c_bound, x_c_bound]=find(mask);
            y_c_bound=y_c_bound-2+min_e_points(2);
            x_c_bound=x_c_bound-2+min_e_points(1);

            y_c_bound(y_c_bound<1)=1;
            x_c_bound(x_c_bound<1)=1;
            y_c_bound(y_c_bound>param.img_s(2))=param.img_s(2);
            x_c_bound(x_c_bound>param.img_s(1))=param.img_s(1);

            b_points=f_img_m(sub2ind(fliplr(param.img_s),y_c_bound,x_c_bound));

            div_img(find(mask))=2^8-1;

            % Seed daughter cells.
            d_seed1=uint16(double(param.tracks(cell_i).cent(end,:))+[cosd(-double(param.tracks(cell_i).ellipse(end,3))) sind(-double(param.tracks(cell_i).ellipse(end,3)))]*double(param.tracks(cell_i).ellipse(end,1))/4);
            d_seed2=uint16(double(param.tracks(cell_i).cent(end,:))-[cosd(-double(param.tracks(cell_i).ellipse(end,3))) sind(-double(param.tracks(cell_i).ellipse(end,3)))]*double(param.tracks(cell_i).ellipse(end,1))/4);

            radius=double(param.tracks(cell_i).ellipse(end,1))/8;

            daughter_mask=zeros(size(mask2));
            daughter_mask=ellipseMatrix(double(d_seed1(2)-min_e_points(2)+2),double(d_seed1(1)-min_e_points(1)+2),radius,radius,0,daughter_mask,1);
            daughter_mask=ellipseMatrix(double(d_seed2(2)-min_e_points(2)+2),double(d_seed2(1)-min_e_points(1)+2),radius,radius,0,daughter_mask,1);
            daughter_mask=daughter_mask&mask2;

            div_img(find(daughter_mask))=0;
            % Segment.
            div_img_mod=imhmin(div_img,2^7);
            div_segment=watershed(div_img_mod);

            [y_c_mid, x_c_mid]=find(~logical(div_segment) & ~mask);
            mid_points=f_img_m(sub2ind(fliplr(param.img_s),y_c_mid-2+min_e_points(2),x_c_mid-2+min_e_points(1)));

            % Collect results.
            div_info=regionprops(div_segment,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');

            % Check if first seed point is within ROI.
            if (d_seed1(2)-min_e_points(2)+2)<1 | (d_seed1(2)-min_e_points(2)+2)>size(div_segment,1) | (d_seed1(1)-min_e_points(1)+2)<1 | (d_seed1(1)-min_e_points(1)+2)>size(div_segment,2)
                d_ind1=0;
            else
                d_ind1=div_segment(d_seed1(2)-min_e_points(2)+2,d_seed1(1)-min_e_points(1)+2);
            end
            % Check if second seed point is within ROI.
            if (d_seed2(2)-min_e_points(2)+2)<1 | (d_seed2(2)-min_e_points(2)+2)>size(div_segment,1) | (d_seed2(1)-min_e_points(1)+2)<1 | (d_seed2(1)-min_e_points(1)+2)>size(div_segment,2)
                d_ind2=0;
            else
                d_ind2=div_segment(d_seed2(2)-min_e_points(2)+2,d_seed2(1)-min_e_points(1)+2);
            end                
            % Check that two unique seeds are in the image.
            if d_ind1==1 | d_ind1==0 | d_ind2==1 | d_ind2==0 | length(div_info)>3
                division_ok=0;
            % Check areas and eccentricities of daughter cells meet
            % conditions.
            else
                A=sort([div_info(d_ind2).Area div_info(d_ind1).Area]);
                ecc=sqrt(1-[div_info(d_ind1).MinorAxisLength div_info(d_ind2).MinorAxisLength].^2./[div_info(d_ind1).MajorAxisLength div_info(d_ind2).MajorAxisLength].^2);

                if (A(2)/A(1)-1)>0.5 | ecc(1)>0.7 | ecc(2)>0.7
                    division_ok=0;
                end
                if mean(double(mid_points))>mean(double(b_points)) & (A(2)/A(1)-1)<0.5
                    division_ok=1;
                end                    
            end
        end
        % Cause split.
        if (~division_ok) & split_ok
            e_points=double(cell2mat(param.tracks(cell_i).bounds{end}'));
            % Remove lengths of boundary points.
            e_points_l=1;
            for b_ind=1:(length(param.tracks(cell_i).bounds{end})-1)
                e_points_l(end+1)=size(param.tracks(cell_i).bounds{end}{b_ind},1);
            end
            % Create image for the watershedder.
            %   - Background is zeros.
            %   - Segment is outlined with boundary values 2^8-1.
            %   - Two seed points are marking locations new daughter cells.

            e_points(cumsum(e_points_l),:)=[];
            min_e_points=min(e_points);
            max_e_points=max(e_points);
            div_w_s=fliplr(max_e_points-min_e_points+1);

            mask=logical(zeros(div_w_s+[2 2]));
            mask(sub2ind(div_w_s+[2 2],e_points(:,2)-min_e_points(2)+2,e_points(:,1)-min_e_points(1)+2))=1;
            mask1=imfill(mask,'holes');
            mask2=mask1-mask;
            [y_c, x_c]=find(mask2);

            div_img=zeros(div_w_s+[2 2]);
            div_img(sub2ind(div_w_s+[2 2],y_c,x_c))=f_img_m(sub2ind(fliplr(param.img_s),y_c-2+min_e_points(2),x_c-2+min_e_points(1)));
            div_img(div_img==(2^8-1))=2^8-2;

            div_img(find(mask))=2^8-1;

            [y_ind, x_ind]=find(fill_img & fill_img1);
            d_seed1=uint16(mean([x_ind y_ind])+v_m+d_min-1);

            [y_ind2, x_ind2]=find((~fill_img) & fill_img1);
            d_seed2=uint16(mean([x_ind2 y_ind2])+v_m+d_min-1);

            radius=double(param.tracks(cell_i).ellipse(end,1))/8;

            daughter_mask=zeros(size(mask2));

            daughter_mask=ellipseMatrix(double(d_seed1(2)-min_e_points(2)+2),double(d_seed1(1)-min_e_points(1)+2),radius,radius,0,daughter_mask,1);
            daughter_mask=ellipseMatrix(double(d_seed2(2)-min_e_points(2)+2),double(d_seed2(1)-min_e_points(1)+2),radius,radius,0,daughter_mask,1);

            daughter_mask=daughter_mask&mask2;

            div_img(find(daughter_mask))=0;
            % Segment.
            div_img_mod=imhmin(div_img,2^7);
            div_segment=watershed(div_img_mod);

            % Collect results.
            div_info=regionprops(div_segment,'centroid','area','perimeter','orientation','majoraxislength','minoraxislength');

            % Check if first seed point is within ROI.
            d_ind1=0;
            if (d_seed1(2)-min_e_points(2)+2)>0 & (d_seed1(2)-min_e_points(2)+2)<=size(div_segment,1) & (d_seed1(1)-min_e_points(1)+2)>0 & (d_seed1(1)-min_e_points(1)+2)<=size(div_segment,2)
                d_ind1=div_segment(d_seed1(2)-min_e_points(2)+2,d_seed1(1)-min_e_points(1)+2);
            end
            % Check if second seed point is within ROI.
            d_ind2=0;
            if (d_seed2(2)-min_e_points(2)+2)>0 & (d_seed2(2)-min_e_points(2)+2)<=size(div_segment,1) & (d_seed2(1)-min_e_points(1)+2)>0 & (d_seed2(1)-min_e_points(1)+2)<=size(div_segment,2)
                d_ind2=div_segment(d_seed2(2)-min_e_points(2)+2,d_seed2(1)-min_e_points(1)+2);
            end                

            % Check that two unique seeds are in the image.
            if d_ind1==1 | d_ind1==0 | d_ind2==1 | d_ind2==0 | d_ind1==d_ind2 | length(div_info)>3
                split_ok=0;
            end
        end

        % If division or split happened modify tracklet structure.
        if ~division_ok & ~split_ok
            continue
        end
        
        % Collect segment information from segmented image
        % replace mother segment with daughter segments

        % Create daughters.
        param.tracks(end+1)=struct('t',uint16(time+1),'A',uint32(div_info(d_ind1).Area),'cent',uint16(div_info(d_ind1).Centroid+min_e_points-2),'ellipse',int16([div_info(d_ind1).MajorAxisLength div_info(d_ind1).MinorAxisLength div_info(d_ind1).Orientation]),'perim',uint32(div_info(d_ind1).Perimeter),'birth',int32(cell_i),'death',int16(0),'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1),'vertices',cell(1));
        param.tracks(end+1)=struct('t',uint16(time+1),'A',uint32(div_info(d_ind2).Area),'cent',uint16(div_info(d_ind2).Centroid+min_e_points-2),'ellipse',int16([div_info(d_ind2).MajorAxisLength div_info(d_ind2).MinorAxisLength div_info(d_ind2).Orientation]),'perim',uint32(div_info(d_ind2).Perimeter),'birth',int32(cell_i),'death',int16(0),'daughters',uint32([0 0]),'neighs',cell(1),'bounds',cell(1),'vertices',cell(1));

        param.tracks(end-1).bounds=cell(1);
        param.tracks(end-1).bounds{end}=cell(0);
        param.tracks(end-1).neighs=cell(1);

        param.tracks(end).bounds=cell(1);
        param.tracks(end).bounds{end}=cell(0);
        param.tracks(end).neighs=cell(1);

        % Update neighbours outer membrane of "old" mother segment.
        for mom_neigh=1:length(param.tracks(cell_i).bounds{end})
            num_boundary=param.tracks(cell_i).bounds{end}{mom_neigh}(1,1);

            ind1_found=zeros(num_boundary,1);
            ind2_found=zeros(num_boundary,1);
            for end_ind=1:num_boundary
                b_point=double(param.tracks(cell_i).bounds{end}{mom_neigh}(1+end_ind,:))-min_e_points+2;
                b_point_neigh=div_segment(b_point(2)-1:b_point(2)+1,b_point(1)-1:b_point(1)+1);
                ind1_found(end_ind)=~isempty(find(b_point_neigh==d_ind1));
                ind2_found(end_ind)=~isempty(find(b_point_neigh==d_ind2));
            end

            % If second daughter cell does not have any common boundary endpoints with
            % current segment, the whole corresponding boundary belongs to the first
            % daughter cell.
            if isempty(find(ind2_found))
                neighbour_track_i=param.tracks(cell_i).neighs{end}(mom_neigh);
                if neighbour_track_i>0
                    neighbour_bound_i=find(param.tracks(neighbour_track_i).neighs{end}==cell_i);
                    param.tracks(neighbour_track_i).neighs{end}(neighbour_bound_i)=length(param.tracks)-1;
                end
                param.tracks(length(param.tracks)-1).neighs{end}(end+1)=neighbour_track_i;
                param.tracks(length(param.tracks)-1).bounds{end}{end+1}=param.tracks(cell_i).bounds{end}{mom_neigh};

            % If first daughter cell does not have any common boundary endpoints with
            % current segment, the whole corresponding boundary belongs to the second
            % daughter cell.
            elseif isempty(find(ind1_found))
                neighbour_track_i=param.tracks(cell_i).neighs{end}(mom_neigh);
                if neighbour_track_i>0
                    neighbour_bound_i=find(param.tracks(neighbour_track_i).neighs{end}==cell_i);
                    param.tracks(neighbour_track_i).neighs{end}(neighbour_bound_i)=length(param.tracks);
                end
                param.tracks(length(param.tracks)).neighs{end}(end+1)=neighbour_track_i;
                param.tracks(length(param.tracks)).bounds{end}{end+1}=param.tracks(cell_i).bounds{end}{mom_neigh};

            % Otherwise the boundary needs to be split into two.
            else
                bound_points=zeros(size(param.tracks(cell_i).bounds{end}{mom_neigh},1)-1,1);
                for b_p_ind=2:size(param.tracks(cell_i).bounds{end}{mom_neigh},1)
                    b_point=double(param.tracks(cell_i).bounds{end}{mom_neigh}(b_p_ind,:))-min_e_points+2;
                    b_point_neigh=div_segment(b_point(2)-1:b_point(2)+1,b_point(1)-1:b_point(1)+1);
                    ind1_ok=~isempty(find(b_point_neigh==d_ind1));
                    ind2_ok=~isempty(find(b_point_neigh==d_ind2));
                    bound_points(b_p_ind-1)=ind1_ok+ind2_ok*2;
                end

                % Form first half of membrane segment.
                b_segment=param.tracks(cell_i).bounds{end}{mom_neigh}([find(ind1_found)+1 ; find(bound_points==3)+1],:);
                b_segment=[size(b_segment,1) 0 ; b_segment ; param.tracks(cell_i).bounds{end}{mom_neigh}(find((bound_points==1) | (bound_points==0))+1,:)];

                neighbour_track_i=param.tracks(cell_i).neighs{end}(mom_neigh);
                if neighbour_track_i>0
                    neighbour_bound_i=find(param.tracks(neighbour_track_i).neighs{end}==cell_i);
                    param.tracks(neighbour_track_i).neighs{end}(neighbour_bound_i)=length(param.tracks)-1;
                    param.tracks(neighbour_track_i).bounds{end}{neighbour_bound_i}=uint16(b_segment);
                end

                param.tracks(length(param.tracks)-1).bounds{end}{end+1}=uint16(b_segment);
                param.tracks(length(param.tracks)-1).neighs{end}(end+1)=neighbour_track_i;

                % Form second half of membrane segment.
                b_segment=param.tracks(cell_i).bounds{end}{mom_neigh}([find(ind2_found)+1 ; find(bound_points==3)+1],:);
                b_segment=[size(b_segment,1) 0 ; b_segment ; param.tracks(cell_i).bounds{end}{mom_neigh}(find((bound_points==2)  | (bound_points==0))+1,:)];

                neighbour_track_i=param.tracks(cell_i).neighs{end}(mom_neigh);
                if neighbour_track_i>0
                    param.tracks(neighbour_track_i).neighs{end}(end+1)=length(param.tracks);
                    param.tracks(neighbour_track_i).bounds{end}{end+1}=uint16(b_segment);
                end

                param.tracks(length(param.tracks)).bounds{end}{end+1}=uint16(b_segment);
                param.tracks(length(param.tracks)).neighs{end}(end+1)=neighbour_track_i;
            end
        end
        
        % Find and create membrane between two new daughter segments.
        b_segment=[0 0];
        [ind_y, ind_x]=find(div_segment==0);
        for b_ind=1:length(ind_y)
            b_neigh=div_segment(ind_y(b_ind)-1:ind_y(b_ind)+1,ind_x(b_ind)-1:ind_x(b_ind)+1);
            ind1_ok=~isempty(find(b_neigh==d_ind1));
            ind2_ok=~isempty(find(b_neigh==d_ind2));
            indb_ok=~isempty(find(b_neigh==1));
            if ~indb_ok
                b_segment(end+1,:)=[ind_x(b_ind) ind_y(b_ind)];
            elseif ind2_ok & ind1_ok
                b_segment(1,1)=b_segment(1,1)+1;
                b_segment=[b_segment(1:b_segment(1,1),:) ; ind_x(b_ind) ind_y(b_ind) ; b_segment((b_segment(1,1)+1):end,:)];
            end
        end
        b_segment(2:end,:)=b_segment(2:end,:)+repmat(min_e_points-2,[size(b_segment,1)-1 1]);

        param.tracks(length(param.tracks)-1).bounds{end}{end+1}=uint16(b_segment);
        param.tracks(length(param.tracks)-1).neighs{end}(end+1)=length(param.tracks);

        param.tracks(length(param.tracks)).bounds{end}{end+1}=uint16(b_segment);
        param.tracks(length(param.tracks)).neighs{end}(end+1)=length(param.tracks)-1;

        b_segment(find(b_segment(:,1)<1 | b_segment(:,1)>param.img_s(1) | b_segment(:,2)<1 | b_segment(:,2)>param.img_s(2)),:)=[];

        s_img(b_segment(:,2)+(b_segment(:,1)-1)*param.img_s(2)+param.img_s(2)*param.img_s(1))=255;

        % Kill mother.
        param.tracks(cell_i).t(end)=[];
        param.tracks(cell_i).A(end)=[];
        param.tracks(cell_i).cent(end,:)=[];
        param.tracks(cell_i).ellipse(end,:)=[];
        param.tracks(cell_i).perim(end)=[];
        if division_ok
            param.tracks(cell_i).death=-3;
        else
            param.tracks(cell_i).death=-8;
        end
        param.tracks(cell_i).daughters=[length(param.tracks)-1 length(param.tracks)];
        param.tracks(cell_i).neighs(end)=[];
        param.tracks(cell_i).bounds(end)=[];
    end
    
    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    tic    
    
    %% (time loop) Filter next image (or read in if already filtered).
    s_img(:,:,1)=s_img(:,:,2);
    f_imgs(:,:,1)=f_imgs(:,:,2);

    if time<param.image_interval(2)-1                
        % Filter and read in first two images if required.
        if isempty(param.f_img_path)
            evaluate_imagej_script(param.bandpass_filter_path,param.imagej_path,...
                [param.o_filename  num2str((param.image_interval_start_in_original_seq+time)*...
                param.file_number_increment+1,'%.04d') '.tif'],temp_path);
            f_imgs(:,:,2)=imread(temp_path);
        % Read in filtered images directly.
        else
            f_imgs(:,:,2)=imread([param.f_img_path num2str((param.image_interval_start_in_original_seq+time)...
                *param.file_number_increment+1,'%.04d') '.tif']);
        end
    end

    slot_i=slot_i+1;
    time_slots(time,slot_i)=toc;
    disp([time_slots(time,:) sum(time_slots(time,:))])

    %% (time loop) Save backup every 20th time point and save in the end.
    if mod(time,20)==0 | time==(param.image_interval(2)-1)
        param.computing_time=time_slots;
        tic
        param.computing_time=time_slots;
        temp_time=param.image_interval(2);
        param.image_interval(2)=time;

        % Saving.
        file_name=[param.seq_name '_segmentation_parameters_'];

        param1=[];   
        file_num=1;

        track_interval=[1 length(param.tracks); length(param.tracks) length(param.tracks)];
        while 1                    
            param1.tracks=param.tracks(track_interval(file_num,1):track_interval(file_num,2));
            size_param=whos('param1');
            size_param.bytes
            while size_param.bytes>(2e9)
                track_interval(file_num,2)=round(track_interval(file_num,2)*0.9);
                track_interval(file_num+1,1)=track_interval(file_num,2)+1;

                param1.tracks=param.tracks(track_interval(file_num,1):track_interval(file_num,2));
                size_param=whos('param1');                            
            end

            param1.track_interval=track_interval(file_num,:);
            save([output_path file_name num2str(file_num)],'param1');                       

            if track_interval(file_num,2)==length(param.tracks)
                break;
            end

            file_num=file_num+1;
            track_interval(end+1,:)=[length(param.tracks) length(param.tracks)];
        end

        clear('param1');

        f_names=fieldnames(param);
        f_names(strcmp(f_names,'tracks'))=[];        
        param0=[];
        for i_ind=1:length(f_names)
            f_name=f_names(i_ind);
            f_name=f_name{1};        
            param0=setfield(param0,f_name,getfield(param,f_name));
        end    
        save([output_path file_name num2str(0)],'param0');   

        clear('param0'); 


        param.image_interval(2)=temp_time;
        disp(['saving backup: ' num2str(toc)]);
    end
        % Write segmentation image.
    imwrite(s_img(:,:,1),[segment_path param.seq_name '_' num2str(time+1,'%.03d') '.tif'],'tif','compression','none');
    param  =  produceVertices(time+1,param,path2segments);
end % time
save([output_path 'param.mat'],'param');
% Delete temporary file.
if exist(temp_path,'file')
    delete(temp_path);
end
