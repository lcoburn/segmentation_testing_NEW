function outputValue=findMatchingNumber(filename,matching_string,varargin)
%findMatchingNumber Read value of a text attribute from a text file.
%   outputValue=findMatchingNumber(filename,matching_string,varargin) find
%   a parameter value from file filename corresponding to the attribute
%   matching_string{1}. Parameter type is specified in matching_string{2}
%   using same syntax as Matlab function sscanf uses. First match is
%   returned by default. Third argument can be used to define non-default
%   matching index. An empty array is returned if match is not found.

% Check if users want first match.
if nargin==3
    count=varargin{1};
% Chose first match by default.
else
    count=1;
end

fid=fopen(filename);
% If no match is found return empty array.
outputValue=[];
while ~feof(fid)
    % Read file line by line
    tline = fgetl(fid);
    % Find location of the attribute.
    str_loc=strfind(tline,matching_string{1});
    if ~isempty(str_loc);
        tline=tline(str_loc:end);
        % Extract parameter value.
        outputValue=sscanf(tline,[matching_string{1} matching_string{2}]);
        % Update count.
        count=count-1;
    end        
    % Stop searching when correct number of matches is found. 
    if count==0
        break;
    end
end
fclose(fid);
