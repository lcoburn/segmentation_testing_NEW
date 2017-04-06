classdef pivVelocityField < handle
    %PIV_VEL_FIELD_STRUCTURE Class for storing, computing and visualising
    %PIV velocity field.
    %   Detailed explanation goes here
    
    properties
        % Absolute path of the velocity field file stored on disk.
        vel_field_file_path;
        % Time index of first time point.
        first_time_index;
        % Velocity field components along x-direction.
        u;
        % Velocity field components along y-direction.
        v;
        % Velocity field vector positions along x-direction.
        x;
        % Velocity field vector positions along y-direction.
        y;
        % Computation settings of the velocity field.
        settings;
    end
    
    methods
        % Constructor method that computes velocity field. Velocity field
        % is computed from image data speficied by image_stack object (an
        % object of experiment_file_folder class). Velocity field is
        % computed for time points specified by interval parameter.
        % Optional fourth parameter settings for the PIV computation (see
        % PIVlab_vel_field_settings function). Optional fifth
        % parameters specifies number of Matlab workers used for the PIV
        % computations.        
        % If interval parameter is empty image_stack parameter is used to
        % determine interval (maximum number of images).
        % If output_path is empty the PIV field is stored into same folder
        % with input images specified by image_stack parameter.
        
        function obj=pivVelocityField(interval,output_path,image_stack,varargin)
            % Use default output file if requested.
            if isempty(output_path)
                obj.vel_field_file_path=[fileparts(fileparts(image_stack.path_name))...
                    filesep 'PIV_field' filesep 'PIV_field.mat'];
            else
                obj.vel_field_file_path=output_path;
            end
            
            % Check if user requested specific PIV computation settings.
            settings_in.preset=2;
            if nargin>3
                if isstruct(varargin{1})
                    settings_in=varargin{1};
                else
                    settings_in.preset=varargin{1};
                end
            end
            
            % Check if user requested specific number of workers.
            num_workers=1;
            if nargin>4
                num_workers=varargin{2};
            end
            num_workers=pivVelocityField.start_parallel_pool(num_workers);

            % Use default interval in case user did not specify interval.
            if isempty(interval)
                interval=[image_stack.frame_index_min image_stack.frame_index_max-1];
            end                        
            obj.first_time_index=interval(1);            

            % Compute PIV between first two frames of image_stack to get
            % settings_out structure that is saved later.
            [~,~,~,~,settings_out]=PIVlab_vel_field_settings(image_stack.get_image(1),...
                image_stack.get_image(1+1),settings_in);
            
            % Prepare file name cell arrays to be used inside parfor loop.
            % Splitting data in this fashion improves performance.
            file_name1=cell(diff(interval)+1,1);
            file_name2=cell(diff(interval)+1,1);            
            for t_ind=interval(1):interval(2)
                file_name1{t_ind-interval(1)+1}=image_stack.get_filename(t_ind);
                file_name2{t_ind-interval(1)+1}=image_stack.get_filename(t_ind+1);
            end
            
            % Prepare PIV compution settings cell arrays to be used inside 
            % parfor loop.
            settings_in_cell=cell(diff(interval)+1,1);            
            for i_ind=1:length(settings_in_cell)
                settings_in_cell{i_ind}=settings_in;
            end
            
            % Define order in which time points of the whole experiment are
            % processed. Simply using parfor loop would yeild arbitrary
            % order.
            loop_inds=1:(interval(2)-interval(1)+1);            
            time_passed_total=0;
            total_time_points=length(loop_inds);
            
            % Keep looping until all the time points are processed.
            while ~isempty(loop_inds)
                tic                
                % Guarantee correct number of parallel workers.
                num_workers_real=pivVelocityField.start_parallel_pool(num_workers);
                
                % Prepare output cell arrays for the parfor loop.
                x_test1=cell(num_workers_real,1);
                y_test1=cell(num_workers_real,1);
                u_original_test1=cell(num_workers_real,1);
                v_original_test1=cell(num_workers_real,1);
                
                % Do not exceed total number of time points.
                if length(loop_inds)<num_workers_real
                    num_workers_real=length(loop_inds);
                end
                
                % Prepare input cell arrays for the parfor loop.
                img_in_inds1=cell(num_workers_real,1);
                img_in_inds2=cell(num_workers_real,1);
                settings_inds=cell(num_workers_real,1);                
                for i_ind=loop_inds(1:min(num_workers_real,end))
                    img_in_inds1{i_ind-loop_inds(1)+1}=file_name1{i_ind};
                    img_in_inds2{i_ind-loop_inds(1)+1}=file_name2{i_ind};
                    settings_inds{i_ind-loop_inds(1)+1}=settings_in_cell{i_ind};
                end
                
                % Compute velocity fields for num_workers_real+1
                % consecutive time points.
                parfor i_ind=1:num_workers_real
                    % Input images.
                    img1=imread(img_in_inds1{i_ind});
                    img2=imread(img_in_inds2{i_ind});
                    
                    % Do the actual velocity field computation.
                    [x1,y1,u_original1,v_original1,~]=...
                        PIVlab_vel_field_settings(img1,img2,settings_inds{i_ind});
                    
                    % Store output data.
                    x_test1{i_ind}=x1;
                    y_test1{i_ind}=y1;
                    u_original_test1{i_ind}=u_original1;
                    v_original_test1{i_ind}=v_original1;
                end
                
                % Store all newly computed velocity fields into current
                % intance of this object.
                for i_ind=loop_inds(1):loop_inds(num_workers_real)                    
                    obj.x{i_ind}=x_test1{i_ind-loop_inds(1)+1};
                    obj.y{i_ind}=y_test1{i_ind-loop_inds(1)+1};
                    obj.u{i_ind}=u_original_test1{i_ind-loop_inds(1)+1};
                    obj.v{i_ind}=v_original_test1{i_ind-loop_inds(1)+1};                    
                end
                
                % Keep track of the time points that not computed yet.
                loop_inds(1:num_workers_real)=[];
                
                % Estimate and print the remaining computation time.
                time_passed_total=time_passed_total+toc;                
                time_left=time_passed_total/(total_time_points-length(loop_inds))*length(loop_inds);
                days_left=floor(time_left/86000);
                hours_left=floor((time_left-days_left*86400)/60/60);
                mins_left=ceil((time_left-days_left*86400-hours_left*3600)/60);
                disp([num2str(total_time_points-length(loop_inds)) ' / ' ...
                    num2str(total_time_points) ' done, time left: ' ...
                    num2str(days_left) ' d ' num2str(hours_left) ' h ' num2str(mins_left)]);                
            end % while ~isempty(loop_inds)
            
            % Save all velocity field data of current experiment into
            % output file. The format is same as PIVlab uses.
            settings=settings_out;
            u_original=obj.u;
            v_original=obj.v;
            x=obj.x;
            y=obj.y;            
            mkdir(fileparts(obj.vel_field_file_path));
            save(obj.vel_field_file_path,'u_original','v_original','x','y','settings');
            
            % Store used velocity fild settings into this object.
            obj.settings=settings_out;
        end
        
        % Read velocity field from disk and store into this object.
        function obj=load_vel_field(obj)
            load(obj.vel_field_file_path)
            obj.u=u_original;
            obj.v=v_original;
            obj.x=x;
            obj.y=y;
        end
        
        % Clear velocity field of this object (this does not delete data
        % from disk).
        function obj=unload(obj)
            obj.u=[];
            obj.v=[];
            obj.x=[];
            obj.y=[];
        end
        
        % Run PIVlab post processing for PIV fields of the whole experiment
        % using num_workers parallel workers. See more details from the
        % PIVlab's PIVlab_commandline function.
        function obj=postprocess(obj,num_workers)
            
            umin = -60; % minimum allowed u velocity
            umax = 60; % maximum allowed u velocity
            vmin = -60; % minimum allowed v velocity
            vmax = 60; % maximum allowed v velocity
            stdthresh=6; % threshold for standard deviation check
            epsilon=0.15; % epsilon for normalized median test
            thresh=3; % threshold for normalized median test
            
            u=obj.u;
            v=obj.v;
            
            % Guarantee correct number of parallel workers.
            num_workers=pivVelocityField.start_parallel_pool(num_workers);
            
            parfor t_ind=1:length(obj.u)                
                u_filtered=u{t_ind};
                v_filtered=v{t_ind};
                
                u_filtered(u_filtered<umin)=NaN;
                u_filtered(u_filtered>umax)=NaN;
                v_filtered(v_filtered<vmin)=NaN;
                v_filtered(v_filtered>vmax)=NaN;
                % stddev check
                meanu=nanmean(nanmean(u_filtered));
                meanv=nanmean(nanmean(v_filtered));
                std2u=nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
                std2v=nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
                minvalu=meanu-stdthresh*std2u;
                maxvalu=meanu+stdthresh*std2u;
                minvalv=meanv-stdthresh*std2v;
                maxvalv=meanv+stdthresh*std2v;
                u_filtered(u_filtered<minvalu)=NaN;
                u_filtered(u_filtered>maxvalu)=NaN;
                v_filtered(v_filtered<minvalv)=NaN;
                v_filtered(v_filtered>maxvalv)=NaN;
                % normalized median check
                % Westerweel & Scarano (2005): Universal Outlier detection for PIV data
                [J,I]=size(u_filtered);
                % medianres=zeros(J,I);
                normfluct=zeros(J,I,2);
                b=1;
                for c=1:2
                    if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
                    for i=1+b:I-b
                        for j=1+b:J-b
                            neigh=velcomp(j-b:j+b,i-b:i+b);
                            neighcol=neigh(:);
                            neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                            med=median(neighcol2);
                            fluct=velcomp(j,i)-med;
                            res=neighcol2-med;
                            medianres=median(abs(res));
                            normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
                        end
                    end
                end
                info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
                u_filtered(info1==1)=NaN;
                v_filtered(info1==1)=NaN;
                
                u{t_ind}=u_filtered;
                v{t_ind}=v_filtered;
                
            end
            
            obj.u=u;
            obj.v=v;
            
            obj.settings.postprocess=1;
            
        end
        
        % Run inpaint_nans for all the time points of the experiment. This
        % will replace all the nans with numbering values in all the
        % velocity fields.
        function obj=inpaint_nans(obj)
            for t_ind=1:length(obj.u)
                obj.u{t_ind}=inpaint_nans(obj.u{t_ind},4);
                obj.v{t_ind}=inpaint_nans(obj.v{t_ind},4);
            end
            
            obj.settings.inpaint_nans=1;
        end
        
        % Pad edges of all the velocity fields of the whole experiment. The
        % edges are padded with the boundary values.
        function obj=pad_matrix_edges(obj,focused_images)
            
            for t_ind=1:length(obj.u)
                obj.x{t_ind}=repmat([0 obj.x{t_ind}(1,:) ...
                    focused_images.img_width+1],[size(obj.x{t_ind},1)+2 1]);
                obj.y{t_ind}=repmat([0 ; obj.y{t_ind}(:,1) ; ...
                    focused_images.img_height],[1 size(obj.y{t_ind},2)+2]);
                
                u_filtered1=zeros(size(obj.u{t_ind})+[2 2]);
                u_filtered1(2:end-1,2:end-1)=obj.u{t_ind};
                u_filtered1([1 end],2:end-1)=u_filtered1([2 end-1],2:end-1);
                u_filtered1(1:end,[1 end])=u_filtered1(1:end,[2 end-1]);
                
                v_filtered1=zeros(size(obj.v{t_ind})+[2 2]);
                v_filtered1(2:end-1,2:end-1)=obj.v{t_ind};
                v_filtered1([1 end],2:end-1)=v_filtered1([2 end-1],2:end-1);
                v_filtered1(1:end,[1 end])=v_filtered1(1:end,[2 end-1]);
                
                obj.u{t_ind}=u_filtered1;
                obj.v{t_ind}=v_filtered1;
            end
            
            obj.settings.pad_matrix_edges=1;
        end
        
        % Post processing steps PIV that is used for cell tracking.
        function obj=postprocess_for_tracking(obj,focused_images,num_workers)
            obj.postprocess(num_workers);
            obj.inpaint_nans;
            obj.pad_matrix_edges(focused_images);
        end
        
        % Visualise PIV field. The velocity field is overlayed with 
        % a stack of background_images (empty for black background). 
        % Following visualisation parameters can be used to alter the ouput:
        % visualisationSettings:
        %   background_image_scale
        %   vector_scale
        %   spatial_averaging_window_s
        %   temporal_averaging_window_s
        %   vector_show_step
        %   folder_out If this parameter is empty a default folder is
        %       generating in to the working directory.
        %   images_to_store
        %   vector_width    
        %   mask_path This an optional parameter for applying a binary mask
        %       for final result.
        %   filter_using_standard_deviation This is an optinal parameter
        % Fourth parameter is optional and describes element RGB colour of
        % velocity field arrows.
        function obj=visualise_vel_field(obj,settings_in,background_images,varargin)
            % Load PIV data if not already in memory.
            if isempty(obj.u)
                disp('Loading velocity fields...')
                obj.load_vel_field;
            end
            
            % Eliminate nans from the fields.
            obj.inpaint_nans();
            
            % Add support for generating overlays
            if ~isfield(settings_in,'overlay')
                settings_in.overlay=1;
            end
            
            % Check if a binary mask is used for filtering out some of the
            % visualised vector.
            shrink=@(x)x(:);
            if isfield(settings_in,'mask_path')
                % Filter the mask to achieve a smoother result that does
                % not fluctuate significantly over time.
                sum_mask=[];
                for t=1:settings_in.images_to_store(end)
                    try
                        sum_mask(t)=mean(shrink(imread([settings_in.mask_path...
                            num2str(t,'%.04d') '.tif'])));
                    catch
                        break
                    end
                end  
                mask_inds=1:length(sum_mask);                
                inds_temp=find(diff(sum_mask)<-0.02)+1;
                while 1
                    inds_temp=find(diff(sum_mask)<-0.02)+1;
                    sum_mask(inds_temp)=sum_mask(inds_temp-1);
                    % These are final mask indices that are used later to
                    % define which mask is used for each individual frame.
                    mask_inds(inds_temp)=mask_inds(inds_temp-1);
                    if isempty(inds_temp)
                        break;
                    end
                end
            end

            % Precompute arrow head dimensions.
            arrow_head_width=settings_in.vector_width*4;
            arrow_head_height=settings_in.vector_width*2;

            % Generate a default output folder if user did not request a
            % folder.
            if isempty(settings_in.folder_out)
                if settings_in.overlay==1;
                    settings_in.folder_out=['.' filesep ...
                        'vel_field_spat_' num2str(settings_in.spatial_averaging_window_s) ...
                        '_temp_' num2str(settings_in.temporal_averaging_window_s) ...
                        '_step_' num2str(settings_in.vector_show_step) ...
                        '_mag_' num2str(settings_in.vector_scale) filesep];
                else
                    settings_in.folder_out=['.' filesep 'imgs'  filesep];
                end
            end            
            mkdir(settings_in.folder_out);
            
            % Load in standard deviations. This is used in case averaged
            % velocity fields where "error" information may exist. 
            % num_points described how many velocity vectors are used in
            % averaging.
            visualise_standard_deviations=0;
            if isfield(settings_in,'filter_using_standard_deviation')
                load(obj.vel_field_file_path,'covariance_uu','covariance_uv',...
                    'covariance_vv','num_points');
                visualise_standard_deviations=1;
            end
            
            % Generate each time point individually.
            for t_ind=settings_in.images_to_store
                % Ignore not existing time points silently.
                if t_ind>length(obj.u)
                    continue
                end
                
                % Check if user requested background.
                if ~isempty(background_images)
                    img_in_org=imadjust(imresize(background_images.get_image(t_ind),...
                        settings_in.background_image_scale))*0.8;
                % Generate back background if not.
                else
                    img_in_org=zeros(obj.y{t_ind}(end,1)*settings_in.background_image_scale+50,...
                        obj.x{t_ind}(1,end)*settings_in.background_image_scale+50,'uint8');
                end                
                % Make background image a colour image in case coloured
                % arrows are requested.
                img_in=img_in_org;
                if nargin>3
                    img_in(:,:,2)=img_in_org;
                    img_in(:,:,3)=img_in_org;
                end
                
                % Make background image transparent.
                if settings_in.overlay==0;
                    img_in(:)=0;
                end
                
                % Temporal averaging indices.
                temp_inds=(t_ind-floor((settings_in.temporal_averaging_window_s-1)/2)):...
                    (t_ind+floor((settings_in.temporal_averaging_window_s-1)/2));
                temp_inds(temp_inds<1)=[];
                temp_inds(temp_inds>length(obj.u))=[];

                % All spatio-temporal data used for averaging.
                smooth_block_u=zeros([size(obj.u{t_ind}) length(temp_inds)]);
                smooth_block_v=zeros([size(obj.v{t_ind}) length(temp_inds)]);
                for i_ind=temp_inds
                    smooth_block_u(:,:,i_ind-temp_inds(1)+1)=obj.u{i_ind};
                    smooth_block_v(:,:,i_ind-temp_inds(1)+1)=obj.v{i_ind};
                end

                % Average in space and in time simultaneously.
                smooth_u=convn(smooth_block_u,ones(settings_in.spatial_averaging_window_s,...
                    settings_in.spatial_averaging_window_s,length(temp_inds))/...
                    settings_in.spatial_averaging_window_s^2/length(temp_inds),'same');
                smooth_v=convn(smooth_block_v,ones(settings_in.spatial_averaging_window_s,...
                    settings_in.spatial_averaging_window_s,length(temp_inds))/...
                    settings_in.spatial_averaging_window_s^2/length(temp_inds),'same');
                
                % Choose which vectors are visualised.
                smooth_x=obj.x{t_ind}(1:settings_in.vector_show_step:end,...
                    1:settings_in.vector_show_step:end);
                smooth_y=obj.y{t_ind}(1:settings_in.vector_show_step:end,...
                    1:settings_in.vector_show_step:end);
                smooth_u=smooth_u(1:settings_in.vector_show_step:end,...
                    1:settings_in.vector_show_step:end,temp_inds==t_ind);
                smooth_v=smooth_v(1:settings_in.vector_show_step:end,...
                    1:settings_in.vector_show_step:end,temp_inds==t_ind);
                
                % Compute eigen values and eigen vectors or covarience
                % matrix.
                if visualise_standard_deviations
                    % First eigen value.
                    eig1=((covariance_uu{t_ind}+covariance_vv{t_ind})+...
                        sqrt((covariance_uu{t_ind}+covariance_vv{t_ind}).^2-...
                        4*((covariance_uu{t_ind}.*covariance_vv{t_ind})-...
                        covariance_uv{t_ind}.^2)))/2;
                    % Second eigen value.
                    eig2=((covariance_uu{t_ind}+covariance_vv{t_ind})-...
                        sqrt((covariance_uu{t_ind}+covariance_vv{t_ind}).^2-...
                        4*((covariance_uu{t_ind}.*covariance_vv{t_ind})-...
                        covariance_uv{t_ind}.^2)))/2;
                    % Orientation of first eigen value.
                    orient=atan2((eig1-covariance_uu{t_ind})./...
                        covariance_uv{t_ind},ones(size(covariance_uu{t_ind})));

                    % Compute standard deviations from variances.
                    eig1=real(sqrt(eig1));
                    eig2=real(sqrt(eig2));    
                    
                    % Choose which "tensors" are visualised.
                    eig1=eig1(1:settings_in.vector_show_step:end,...
                        1:settings_in.vector_show_step:end);
                    eig2=eig2(1:settings_in.vector_show_step:end,...
                        1:settings_in.vector_show_step:end);
                    orient=orient(1:settings_in.vector_show_step:end,...
                        1:settings_in.vector_show_step:end);
                    % num_points is number of velocity field vectors used
                    % for the averaging.
                    num_points_use=num_points{t_ind}(1:settings_in.vector_show_step:end,...
                        1:settings_in.vector_show_step:end);
                    
                    % Filter vectors based of the standard deviation of
                    % larger eigen value. Currently shows only vectors with
                    % larger magnitude than the corresponding larger
                    % standard deviation.
                    colour_factor=1-double(((eig1*settings_in.filter_using_standard_deviation)>...
                        (sqrt(smooth_u.^2+smooth_v.^2))) | num_points_use<=1)*1.0;
                % Otherwise all the vectors have constant colour.
                else
                    colour_factor=ones(size(smooth_u));
                end                 
                
                % Filter velocity fields based on the optional mask and
                % predetermined mask_inds.
                if isfield(settings_in,'mask_path')
                    mask_img=imread([settings_in.mask_path num2str(mask_inds(t_ind),'%.04d') '.tif']);
                    filter_inds=ones(numel(smooth_x),1);                
                    filter_inds(logical(filter_inds))=mask_img(round(...
                        smooth_x(logical(filter_inds))-1)*size(mask_img,1)...
                        +round(smooth_y(logical(filter_inds))));                
                    smooth_u(~logical(filter_inds))=0;
                    smooth_v(~logical(filter_inds))=0;
                end
                
                % The actual velocity field visualisation. Each vector is
                % drawn individually by modifying the backgroung bitmap.
                for i_ind=1:numel(smooth_x)
                    % Generate an arrow as 7 pairs of x-y points.
                    vecs=[-settings_in.vector_width/2 0; ...
                        settings_in.vector_width/2 0; ...
                        settings_in.vector_width/2 max(sqrt(smooth_u(i_ind)^2+smooth_v(i_ind)^2)*settings_in.vector_scale-arrow_head_height,0); ...
                        arrow_head_width/2 max(sqrt(smooth_u(i_ind)^2+smooth_v(i_ind)^2)*settings_in.vector_scale-arrow_head_height,0); ...
                        0 sqrt(smooth_u(i_ind)^2+smooth_v(i_ind)^2)*settings_in.vector_scale; ...
                        -arrow_head_width/2 max(sqrt(smooth_u(i_ind)^2+smooth_v(i_ind)^2)*settings_in.vector_scale-arrow_head_height,0); ...
                        -settings_in.vector_width/2 max(sqrt(smooth_u(i_ind)^2+smooth_v(i_ind)^2)*settings_in.vector_scale-arrow_head_height,0)]*settings_in.background_image_scale;
                    
                    % Rotate the arrow into correct orientation.
                    vec_angle=atan2(smooth_v(i_ind),smooth_u(i_ind))-pi/2;
                    vecs=vecs*[cos(vec_angle) sin(vec_angle);-sin(vec_angle) cos(vec_angle)];
                    
                    % Translate the arrow into correct location.
                    min_shift=min(vecs);
                    vecs(:,1)=vecs(:,1)-min_shift(1);
                    vecs(:,2)=vecs(:,2)-min_shift(2);
                    
                    % Find set of pixels that are inside of the arrow
                    % polygon.
                    pix_img_s=ceil(max(vecs));
                    [xp, yp]=meshgrid(1:pix_img_s(1),1:pix_img_s(2));                    
                    in_pol=inpolygon(xp,yp,vecs(:,1),vecs(:,2));
                    
                    % Take the output image scale into account.
                    xp=round(xp+min_shift(1)+smooth_x(i_ind)*settings_in.background_image_scale);
                    yp=round(yp+min_shift(2)+smooth_y(i_ind)*settings_in.background_image_scale);
                    
                    % Take care of the boundaries.
                    in_pol(xp<1)=0;
                    in_pol(xp>size(img_in_org,2))=0;                    
                    in_pol(yp<1)=0;
                    in_pol(yp>size(img_in_org,1))=0;
                    
                    % Compute linear indices in to the final image.
                    lin_inds=yp(in_pol(:)==1)+(xp(in_pol(:)==1)-1)*(size(img_in_org,1));

                    % Coloured arrows.
                    if nargin>3
                        img_in(lin_inds)=varargin{1}(1);
                        img_in(lin_inds+(numel(img_in_org*1)))=varargin{1}(2);
                        img_in(lin_inds+(numel(img_in_org*2)))=varargin{1}(3);
                    % Grey scale arrow.
                    else
                        img_in(lin_inds)=round(255*colour_factor(i_ind));
                    end                                        
                end                
                imwrite(img_in,[settings_in.folder_out 'img_' num2str(t_ind,'%.04d') '.tif'],'tif');                                
            end % t_ind time loop
        end % visualise function definition
    end % methods
    
    methods(Static)
        % Start parallel pool with nWorkers workers if not running already.
        function nWorkers=start_parallel_pool(nWorkers)
            pool=gcp('nocreate');
            if ~isempty(pool)
                if pool.NumWorkers~=nWorkers
                    delete(gcp('nocreate'))
                    pool=parpool(nWorkers);        
                end
            else
                pool=parpool(nWorkers);    
            end
            nWorkers=pool.NumWorkers;
        end   
    end
end
