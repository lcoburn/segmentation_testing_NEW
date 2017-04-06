function evaluate_imagej_script(script_file,imagej_path,varargin)
% Options for ImageJ. Do not show graphical user interface.
options='--no-splash -batch';

% All variable arguments are passed as input arguments for the ImageJ
% script.
args=[];
l=length(varargin);
if l>0
    args='"';
    for i=1:length(varargin)
        args=[args varargin{i} ' '];
    end
    args=[args '"'];
end

% Form ImageJ command and run in the system.
command=[imagej_path ' ' options ' ' script_file ' ' args];
system(command);
end
