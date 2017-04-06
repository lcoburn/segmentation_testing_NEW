function guillermo_evaluate_imagej_script_windows(script_file,varargin)
imagej_path='X:\mDrives\storage2\Software\Fiji.app\ImageJ-win64.exe';
options='--no-splash -batch';

args=[];
%     if ~isempty(varargin)
%         args=varargin{1};
%     end
l=length(varargin);
if l>0
    args='"';
    for i=1:length(varargin)
%         args=[args varargin{i} '\t'];
        args=[args varargin{i} ' '];
    end
    args=[args '"'];
end
command=[imagej_path ' ' options ' ' script_file ' ' args];
%     disp(command)

system(command);
% ./ImageJ-linux64 --no-splash -batch /mnt/cjwsmb/Antti/time_sequences_2d_tracking/macros/filter_stack_20X_linux.ijm /mnt/cjwsmb/Antti/time_sequences_2d_tracking/No122_241111/segmentation_revisited/No122_241111KS20X_F131-150_S7.tif

% antti_evaluate_imagej_script('/mnt/cjwsmb/Antti/time_sequences_2d_tracking/macros/filter_stack_20X_linux.ijm',...
%     '/mnt/cjwsmb/Antti/time_sequences_2d_tracking/No122_241111/segmentation_revisited/No122_241111KS20X_F131-150_S7.tif');

% antti_evaluate_imagej_script('/mnt/cjwsmb/Antti/time_sequences_2d_tracking/macros/filter_stack_20X_linux.ijm',...
%     '/mnt/cjwsmb/Antti/time_sequences_2d_tracking/No122_241111/whole_pipeline_F101-300_S7/No122_241111KS20X_F101-300_S7.tif');

end

%/opt/Fiji.app/ImageJ-linux64 --no-splash -batch /mnt/LocalStorage/TwoTerabyte/antti/experiments/epi-fluorescent_NEW/exp18_10x/0001_stack_focusing_imagej/focuse_many_stacks.ijm
