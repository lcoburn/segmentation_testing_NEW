function run_imagej_script_windows(imagej_path,options,script_file,varargin)

if ~exist('imagej_path','var') || isempty(imagej_path);
    imagej_path='X:\mDrives\storage2\Software\Fiji.app\ImageJ-win64.exe';
end

if ~exist('options','var') || isempty(options);
    options='--headless --no-splash -batch';
end

args=[];
l=length(varargin);
if l>0
    args='"';
    for i=1:length(varargin)
        args=[args varargin{i} ' '];
    end
    args=[args '"'];
end
args = strrep(args,'\','\\');
command=[imagej_path ' ' options ' ' script_file ' ' args];
disp(command)
system(command);
end