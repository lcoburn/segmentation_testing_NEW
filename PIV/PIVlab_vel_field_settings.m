function [x,y,u_filtered,v_filtered,settings]=...
    PIVlab_vel_field_settings(image1,image2,settings)
%PIVlab_vel_field_settings Compute PIV velocity field between two images
% using given settings.
%   Velocity field is computed from image1 into image2. If settings struct 
%   has preset field one of the below preset sets are used for settings.
%   Velocity field vector positions are returned in x and y. The x and y
%   component are returned in variabled u and v, respectively. The used
%   settings are returned in the settings variable.

% Preset for "Antti's" tracking settings.
if isfield(settings,'preset') & settings.preset==1
    settings.postprocess=1;
    settings.pad_matrix_edges=1;
    settings.additional_preprocess=1;
    settings.inpaint_nans=1;

    s = cell(10,2); % To make it more readable, let's create a "settings table"
    %Parameter                       %Setting           %Options
    s{1,1}= 'Int. area 1';           s{1,2}=64;         % window size of first pass
    s{2,1}= 'Step size 1';           s{2,2}=32;         % step of first pass
    s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';         s{6,2}=1;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';           s{7,2}=32;         % second pass window size
    s{8,1}= 'Int. area 3';           s{8,2}=16;         % third pass window size
    s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
    s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower
    settings.s=s;
    
    %% Standard image preprocessing settings
    p = cell(8,1);
    %Parameter                       %Setting           %Options
    p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
    p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
    p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
    p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
    p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
    p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
    p{7,1}= 'Clipping thresh.';      p{7,2}=0;          % 0-255 clipping threshold
    p{8,1}= 'Intensity Capping';     p{8,2}=0;          % 1 = enable intensity capping, 0 = disable    
    settings.p=p;

% Preset I for "Emil's" velocity field settings.
elseif isfield(settings,'preset') & settings.preset==2
    settings.postprocess=0;
    settings.pad_matrix_edges=0;
    settings.additional_preprocess=0;   
    settings.inpaint_nans=0;
    
    p = cell(8,1);
    %Parameter                       %Setting           %Options
    p{1,1}= 'ROI';                   p{1,2}=[];     % same as in PIV settings
    p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
    p{3,1}= 'CLAHE size';            p{3,2}=20;         % CLAHE window size
    p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
    p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
    p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
    p{7,1}= 'Clipping thresh.';      p{7,2}=0;          % 0-255 clipping threshold
    p{8,1}= 'Intensity Capping';     p{8,2}=0;          % 1 = enable intensity capping, 0 = disable        
    settings.p=p;
    
    s = cell(10,2); % To make it more readable, let's create a "settings table"
    %Parameter                       %Setting           %Options
    s{1,1}= 'Int. area 1';           s{1,2}=128;         % window size of first pass
    s{2,1}= 'Step size 1';           s{2,2}=64;         % step of first pass
%     s{2,1}= 'Step size 1';           s{2,2}=16;         % step of first pass
    s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';         s{6,2}=3;          % 1-4 nr. of passes
%     s{6,1}= 'Nr. of passes';         s{6,2}=1;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';           s{7,2}=64;         % second pass window size
    s{8,1}= 'Int. area 3';           s{8,2}=32;         % third pass window size
    s{9,1}= 'Int. area 4';           s{9,2}=[];         % fourth pass window size
    s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower        
    settings.s=s;    
% Preset II for "Emil's" velocity field settings.
elseif isfield(settings,'preset') & settings.preset==3
    settings.postprocess=0;
    settings.pad_matrix_edges=0;
    settings.additional_preprocess=0;   
    settings.inpaint_nans=0;
    
    p = cell(8,1);
    %Parameter                       %Setting           %Options
    p{1,1}= 'ROI';                   p{1,2}=[];     % same as in PIV settings
    p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
    p{3,1}= 'CLAHE size';            p{3,2}=20;         % CLAHE window size
    p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
    p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
    p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
    p{7,1}= 'Clipping thresh.';      p{7,2}=0;          % 0-255 clipping threshold
    p{8,1}= 'Intensity Capping';     p{8,2}=0;          % 1 = enable intensity capping, 0 = disable        
    settings.p=p;
    
    s = cell(10,2); % To make it more readable, let's create a "settings table"
    %Parameter                       %Setting           %Options
    s{1,1}= 'Int. area 1';           s{1,2}=64;         % window size of first pass
    s{2,1}= 'Step size 1';           s{2,2}=32;         % step of first pass
    s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';         s{6,2}=2;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';           s{7,2}=32;         % second pass window size
    s{8,1}= 'Int. area 3';           s{8,2}=[];         % third pass window size
    s{9,1}= 'Int. area 4';           s{9,2}=[];         % fourth pass window size
    s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower        
    settings.s=s;      
% Preset for image merging    
elseif isfield(settings,'preset') & settings.preset==4
    settings.postprocess=0;
    settings.pad_matrix_edges=1;
    settings.additional_preprocess=1;
    settings.inpaint_nans=1;
        
    s = cell(10,2); % To make it more readable, let's create a "settings table"
    %Parameter                       %Setting           %Options
    s{1,1}= 'Int. area 1';           s{1,2}=128;         % window size of first pass
    s{2,1}= 'Step size 1';           s{2,2}=32;         % step of first pass
    s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';         s{6,2}=1;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';           s{7,2}=32;         % second pass window size
    s{8,1}= 'Int. area 3';           s{8,2}=16;         % third pass window size
    s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
    s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower
    settings.s=s;  

    %% Standard image preprocessing settings
    p = cell(8,1);
    %Parameter                       %Setting           %Options
    p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
    p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
    p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
    p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
    p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
    p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
    p{7,1}= 'Clipping thresh.';      p{7,2}=0;          % 0-255 clipping threshold
    p{8,1}= 'Intensity Capping';     p{8,2}=0;          % 1 = enable intensity capping, 0 = disable    
    settings.p=p;
    
elseif isfield(settings,'preset') & settings.preset==5
    settings.postprocess=1;
    settings.pad_matrix_edges=1;
    settings.additional_preprocess=1;
    settings.inpaint_nans=1;

    s = cell(10,2); % To make it more readable, let's create a "settings table"
    %Parameter                       %Setting           %Options
    s{1,1}= 'Int. area 1';           s{1,2}=128;         % window size of first pass
    s{2,1}= 'Step size 1';           s{2,2}=32;         % step of first pass
    s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';         s{6,2}=1;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';           s{7,2}=32;         % second pass window size
    s{8,1}= 'Int. area 3';           s{8,2}=16;         % third pass window size
    s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
    s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower
    settings.s=s;
    
    %% Standard image preprocessing settings
    p = cell(8,1);
    %Parameter                       %Setting           %Options
    p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
    p{2,1}= 'CLAHE';                 p{2,2}=1;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
    p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
    p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
    p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
    p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
    p{7,1}= 'Clipping thresh.';      p{7,2}=0;          % 0-255 clipping threshold
    p{8,1}= 'Intensity Capping';     p{8,2}=0;          % 1 = enable intensity capping, 0 = disable    
    settings.p=p;    
end

% Set default values for various processing steps.
if ~isfield(settings,'additional_preprocess');
    settings.additional_preprocess=0;
end
if ~isfield(settings,'postprocess') 
    settings.postprocess=0;
end
if ~isfield(settings,'pad_matrix_edges') 
    settings.pad_matrix_edges=0;
end
if ~isfield(settings,'inpaint_nans') 
    settings.inpaint_nans=0;
end

%% Use standard PIV Settings if no settings were provided by the user.
if ~isfield(settings,'s') 
    s = cell(10,2); % To make it more readable, let's create a "settings table"
    %Parameter                       %Setting           %Options
    s{1,1}= 'Int. area 1';           s{1,2}=128;         % window size of first pass
    s{2,1}= 'Step size 1';           s{2,2}=64;         % step of first pass
    s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
    s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
    s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
    s{6,1}= 'Nr. of passes';         s{6,2}=3;          % 1-4 nr. of passes
    s{7,1}= 'Int. area 2';           s{7,2}=64;         % second pass window size
    s{8,1}= 'Int. area 3';           s{8,2}=32;         % third pass window size
    s{9,1}= 'Int. area 4';           s{9,2}=[];         % fourth pass window size
    s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower
else
    s=settings.s;
end

%% Preprocessing settings
% Scale both input images into mutual min/max range.
if settings.additional_preprocess    
    max_val=max([image1(:);image2(:)]);
    min_val=min([image1(:);image2(:)]);
    image1=uint8((image1-min_val)*(255/double(max_val-min_val+1)));
    image2=uint8((image2-min_val)*(255/double(max_val-min_val+1)));
end

% Preprocessing of PIVlab.
if isfield(settings,'p') 
    p=settings.p;
    image1 = PIVlab_preproc (image1,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
    image2 = PIVlab_preproc (image2,p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});
end

%% Compute velocity field.
assignin('base', 's', s)
[x, y, u_filtered, v_filtered, ~] = piv_FFTmulti_modified(image1,image2,s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});

%% PIVlab postprocessing loop
if settings.postprocess
    % Settings
    
    if ~isfield(settings,'max_length')
        settings.max_length=60;
    end
    
    umin = -settings.max_length; % minimum allowed u velocity
    umax = settings.max_length; % maximum allowed u velocity
    vmin = -settings.max_length; % minimum allowed v velocity
    vmax = settings.max_length; % maximum allowed v velocity
    stdthresh=6; % threshold for standard deviation check
    epsilon=0.15; % epsilon for normalized median test
    thresh=3; % threshold for normalized median test

    % % % typevector_filtered=typevector;
    % vellimit check
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
end

% Interpolate missing data (PIVlab).
if settings.inpaint_nans
    u_filtered=inpaint_nans(u_filtered,4);
    v_filtered=inpaint_nans(v_filtered,4);
end

% Pad the edeges.
if settings.pad_matrix_edges
    x=repmat([0 x(1,:) size(image1,2)+1],[size(x,1)+2 1]);
    y=repmat([0 ; y(:,1) ; size(image1,1)],[1 size(y,2)+2]);

    u_filtered1=zeros(size(u_filtered)+[2 2]);
    u_filtered1(2:end-1,2:end-1)=u_filtered;
    u_filtered1([1 end],2:end-1)=u_filtered1([2 end-1],2:end-1);
    u_filtered1(1:end,[1 end])=u_filtered1(1:end,[2 end-1]);

    v_filtered1=zeros(size(v_filtered)+[2 2]);
    v_filtered1(2:end-1,2:end-1)=v_filtered;
    v_filtered1([1 end],2:end-1)=v_filtered1([2 end-1],2:end-1);
    v_filtered1(1:end,[1 end])=v_filtered1(1:end,[2 end-1]);
    
    u_filtered=u_filtered1;
    v_filtered=v_filtered1;
end

% % clearvars -except x y u_filtered v_filtered settings
% % clearvars -global -except x y u_filtered v_filtered settings
