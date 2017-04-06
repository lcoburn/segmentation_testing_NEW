function runTracking
clc
%testTracking Test function for tracking an image sequence in time
%dependent manner.

% Remember the current search paths configuration and restore to default.
pathStore=path;
restoredefaultpath;

% Add misc functions into the Matlab path.
addpath([pwd filesep 'misc'])
% PIV.
addpath([pwd filesep 'PIV'])
% Tracking.
addpath([pwd filesep 'tracking'])

% Reset output folder.
outputFolderName=[pwd filesep 'output' filesep 'testOutputTracking'];
% if exist(outputFolderName,'dir')
%     rmdir(outputFolderName,'s')
% end
mkdir(outputFolderName)

% Input parameters. See description from the tracking function.
param=struct('o_filename',[pwd filesep 'fabulousTestData' filesep 'img_'],...
    'file_number_increment',1,...
    'image_interval',[1 3],...
    'image_interval_start_in_original_seq',1,...
    'seq_name','img',...
    'clear_boundaries',1,...
    'f_img_path',[],...
    'piv_path',[],...
    'bandpass_filter_path',[pwd filesep 'misc' filesep 'bandpass_filter_z-stack.ijm'],...
    'imagej_path',[pwd filesep 'Fiji.app' filesep 'Contents/MacOS/ImageJ-macosx']);

% Check if ImageJ is found. Terminate if not.
if ~isempty(param.imagej_path)
    if ~exist(param.imagej_path,'file')
        disp(['ImageJ not found from: ' param.imagej_path '. Stopping the script.']);
        return;
    end
end

output_path=[outputFolderName filesep];

% Run tracking function.
trackingWithVerticesAndGhost(param,output_path);
% tracking(param,output_path)

% save param
load(['output' filesep 'testOutputTracking' filesep 'param.mat'],'param');

% Remove misc functions from the Matlab path.
rmpath([pwd filesep 'misc'])
% PIV.
rmpath([pwd filesep 'PIV'])
% Tracking.
rmpath([pwd filesep 'tracking'])

% Revert back to user's own search paths.
path(path,pathStore);
