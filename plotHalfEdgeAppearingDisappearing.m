% plot half edge on data
% plot half edge results
clear all, close all, clc

outPath = pwd;

load('half_edge_data_cells_9.mat')
load('half_edge_data_he_9.mat')
load('half_edge_data_vert_9.mat')
load 'Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testOutputTracking\param.mat'

% imPath='Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testOutputTracking\segments\img_';

outPath = ['Y:\mDrives\storage4\Guillermo\segmentation_testing\output\testHalfEdge\img\'];
mkdir(outPath)

eR = expReader('X:\mDrives\storage4\Guillermo\segmentation_testing\fabulousTestData\');
for t=1:size(C1,2)-1

    HEa = half_edges1(t).HE;
    HEb = half_edges1(t+1).HE;
    
    Va = verts1(t).V;
    Vb = verts1(t+1).V;
    
    V1 = Va(HEa(:,8) == 0,:); %disappearing edges
    V2 = Vb(HEb(:,7) == 0,:); %appearing edges
    
    xa = [Va(:,1) Va(:,3)];
    ya = [Va(:,2) Va(:,4)];
    
	xb = [Vb(:,1) Vb(:,3)];
    yb = [Vb(:,2) Vb(:,4)];
   
    x1 = [V1(:,1) V1(:,3)];
    y1 = [V1(:,2) V1(:,4)];
    
	x2 = [V2(:,1) V2(:,3)];
    y2 = [V2(:,2) V2(:,4)];
    % create the image and save it
% 	I = imread([imPath num2str(t,'%03d') '.tif']);
    fig=figure('Visible', 'off','units','normalized','outerposition',[0 0 1 1]); % overwriting the figure avoids filling the memory
%     fig = figure(1);
    subplot(1,2,1)
    imshow(eR.currentImage);
    hold on
    plot(xa',ya',':','lineWidth', 1, 'Color', 'c')
    plot(x1',y1',':','lineWidth', 1, 'Color', 'm')
    
    subplot(1,2,2)
    eR.step
    imshow(eR.currentImage);
    hold on
    plot(xb',yb',':','lineWidth', 1, 'Color', 'c')
    plot(x2',y2',':','lineWidth', 1, 'Color', 'y')
    
    frame = getframe(fig);
    imwrite(frame.cdata,[outPath filesep 'img_' num2str(t,'%04d') '.tif']);

    clear fig
end

disp(['Tracking finished results in: ' outPath])

% generate the movie. ffmpeg must be added to the path
system(['ffmpeg -r 10 -i ' outPath filesep 'img_%04d.tif'...
        ' -vcodec msmpeg4 -vf scale=1920:-1 -q:v 8 ' outPath filesep...
        'edgeTrack_' num2str(hE2F) '_t0_' num2str(t0) '.avi']);

 load handel; sound(y,Fs) % program finished 