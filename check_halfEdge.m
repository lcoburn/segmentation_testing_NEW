% plot half edge results
clear all, close all, clc

outPath = pwd;

load('half_edge_data_cells_9.mat')
load('half_edge_data_he_9.mat')
load('half_edge_data_vert_9.mat')
load 'Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testOutputTracking\param.mat'

imPath='Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testOutputTracking\segments\img_';

timePointsAfterDead = 5;
t0 = 1; % time point to sample the first half edge
hE2F = datasample(half_edges1(t0).HE(:,1),1); % half edge to follow

% find cells of halg edge to follow
c1 = half_edges1(t0).HE(hE2F,2);
c2 = half_edges1(t0).HE(hE2F,3);

% create the output folder
outPath = [outPath filesep 'edgeTrack_' num2str(hE2F) '_t0_' num2str(t0)];
mkdir(outPath)

alive = true;
for t=1:size(C1,2)
    
    HE = half_edges1(t).HE;
    V = verts1(t).V;
    
    x = [V(:,1) V(:,3)];
    y = [V(:,2) V(:,4)];
    
    centr1 = param.tracks(c1).cent(param.tracks(c1).t==t,:);
    centr2 = param.tracks(c2).cent(param.tracks(c2).t==t,:);
    centrs = [centr1; centr2];
    
    if alive;
        % find coordinates of half edge to follow
        hE2Fx = [V(hE2F,1) V(hE2F,3)];
        hE2Fy = [V(hE2F,2) V(hE2F,4)];
        
        % find hE2F in the next time frame
        hE2F = HE(hE2F,8); % last column
        
        if hE2F == 0; alive = false; end
    else
        timePointsAfterDead = timePointsAfterDead-1;
    end
    
    % create the image and save it
	I = imread([imPath num2str(t,'%03d') '.tif']);
    fig=figure('Visible', 'off'); % overwriting the figure avoids filling the memory

    imshow(I)
    hold on
    plot(x',y', 'lineWidth', 1, 'Color', 'c')
    plot(centrs(:,1), centrs(:,2), '.y', 'markersize', 10)
    
    if alive; plot(hE2Fx',hE2Fy', 'lineWidth', 2, 'Color', 'r'); end
    
    frame = getframe(fig);
    imwrite(frame.cdata,[outPath filesep 'img_' num2str(t,'%04d') '.tif']);
    clear fig
    
    if timePointsAfterDead < 0; break; end
end

disp(['Tracking finished results in: ' outPath])

% generate the movie. ffmpeg must be added to the path
system(['ffmpeg -r 10 -i ' outPath filesep 'img_%04d.tif'...
        ' -vcodec msmpeg4 -vf scale=1920:-1 -q:v 8 ' outPath filesep...
        'edgeTrack_' num2str(hE2F) '_t0_' num2str(t0) '.avi']);

 load handel; sound(y,Fs) % program finished 