function runPivVelocityField
%testComputeAndVisualisePiv Test function for computing and visualising PIV
%   velocity field for test data set.

% Remember the current search paths configuration and restore to default.
pathStore=path;
restoredefaultpath;

% PIV.
addpath([pwd filesep 'PIV'])
% Helper functions.
addpath([pwd filesep 'misc'])

% Reset output folder.
outputFolderName=[pwd filesep 'output' filesep 'testPivVelocityField'];
if exist(outputFolderName,'dir')
    rmdir(outputFolderName,'s')
end
mkdir(outputFolderName)

% Generate a file object for the input image data.
inputFileName=[pwd filesep 'fabulousTestData' filesep];
nTimePoints=61;
focusedImages=experiment_file_folder('normal','tif',2,{'img_' '.tif'},...
    inputFileName,1,nTimePoints,'%.04d',1,[],[],[],[],1);

% Compute PIV velocity field.
timeInterval=[1 61-1];
nCores=6;
settingsPresetGroup=2;
outputPath=[outputFolderName filesep 'PIV_field' filesep ...
    'PIV_field.mat'];
pivField=pivVelocityField(timeInterval,outputPath,focusedImages,...
    settingsPresetGroup,nCores);

% Visualise velocity field.
visualisationSettings=[];
visualisationSettings.background_image_scale=1;
visualisationSettings.vector_scale=10;
visualisationSettings.spatial_averaging_window_s=5;
visualisationSettings.temporal_averaging_window_s=6;
visualisationSettings.vector_show_step=5;
visualisationSettings.folder_out=[outputFolderName filesep 'overlay' filesep];
visualisationSettings.images_to_store=1:2000;
visualisationSettings.vector_width=5;
pivField.visualise_vel_field(visualisationSettings,focusedImages);

system(['ffmpeg -r 10 -i ' visualisationSettings.folder_out 'img_%04d.tif'...
    ' -vcodec msmpeg4 -vf scale=1920:-1 -q:v 8 ' outputFolderName filesep...
    'velocity_field.avi']);

% Save velocity field.
pivField.unload;
save([outputFolderName filesep 'velocityFieldObject'],'pivField')

% PIV
rmpath([pwd filesep 'PIV'])
% Helper functions.
rmpath([pwd filesep 'misc'])

% Revert back to user's own search paths.
path(path,pathStore);
