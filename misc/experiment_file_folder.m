classdef experiment_file_folder 
    %EXPERIMENT_FILE_FOLDER helper class for dealing with time point files 
    %and z-stack files of an experiment.
    %   This class enables convienient file name generation and image
    %   reading.
    
    properties
        description;
        image_format;
        number_file_pieces;
        file_name_pieces;
        path_name;
        frame_index_min;
        frame_index_max;
        frame_index_format;
        frame_index_step;        
        z_index_min;
        z_index_max;
        z_index_format;
        z_index_step;    
        index_locs_in_name; % time first, then z (e.g. [1 2])
        img_width;
        img_height;
        velocity_field=[];
        stitch_information=[];
    end
    
    methods
        % Constructor. Store input parameters into the object.
        function obj=experiment_file_folder(description,image_format,number_file_pieces,file_name_pieces,path_name,frame_index_min,frame_index_max,frame_index_format,frame_index_step,z_index_min,z_index_max,z_index_format,z_index_step,index_locs_in_name)
            obj.description=description;
            obj.image_format=image_format;
            obj.number_file_pieces=number_file_pieces;
            obj.file_name_pieces=file_name_pieces;
            obj.path_name=path_name;
            obj.frame_index_min=frame_index_min;
            obj.frame_index_max=frame_index_max;
            obj.frame_index_format=frame_index_format;
            obj.frame_index_step=frame_index_step;            
            obj.z_index_min=z_index_min;
            obj.z_index_max=z_index_max;
            obj.z_index_format=z_index_format;
            obj.z_index_step=z_index_step;            
            obj.index_locs_in_name=index_locs_in_name;
                
            if isempty(z_index_min)
                file_name=[obj.path_name obj.file_name_pieces{1} num2str(frame_index_min,obj.frame_index_format) obj.file_name_pieces{2}];
            else
                if obj.index_locs_in_name(1)==1
                    file_name_numbers{1}=num2str(frame_index_min,obj.frame_index_format);
                    file_name_numbers{2}=num2str(z_index_min,obj.z_index_format);
                else
                    file_name_numbers{2}=num2str(frame_index_min,obj.frame_index_format);
                    file_name_numbers{1}=num2str(z_index_min,obj.z_index_format);
                end
                file_name=[obj.path_name obj.file_name_pieces{1} file_name_numbers{1} obj.file_name_pieces{2} file_name_numbers{2} obj.file_name_pieces{3}];
            end                                    

            if ~strcmp(obj.image_format,'raw')
                info=imfinfo(file_name);
                obj.img_width=info.Width;
                obj.img_height=info.Height;
            end
        end
        
        % Construct file name string based on time and z indices.
        function file_name=get_filename(varargin) % time first, then z
            obj=varargin{1};
            time=(varargin{2}-1)*obj.frame_index_step+obj.frame_index_min;
            if nargin==2
                file_name=[obj.path_name obj.file_name_pieces{1} num2str(time,obj.frame_index_format) obj.file_name_pieces{2}];
            elseif nargin==3
                if obj.index_locs_in_name(1)==1
                    file_name_numbers{1}=num2str(time,obj.frame_index_format);
                    file_name_numbers{2}=num2str(varargin{3},obj.z_index_format);
                else
                    file_name_numbers{2}=num2str(time,obj.frame_index_format);
                    file_name_numbers{1}=num2str(varargin{3},obj.z_index_format);                    
                end
                file_name=[obj.path_name obj.file_name_pieces{1} file_name_numbers{1} obj.file_name_pieces{2} file_name_numbers{2} obj.file_name_pieces{3}];
            end                        
        end        
        
        % Read in and return image using a file name.
        function img=get_image(varargin)
            obj=varargin{1};
            if nargin==2
                filename=obj.get_filename(varargin{2});
            elseif nargin==3
                filename=obj.get_filename(varargin{2},varargin{3});
            end
            % Sometimes imread fails when reading over network.
            read_ok=0;
            while(~read_ok)
                try
                    img=imread(filename);
                    read_ok=1;
                catch ME                    
                    disp(['error in experiment_file_folder method in reading image:'  ME.message])
                    pause(1)
                end
            end
        end
        
    end
    
end

