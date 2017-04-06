function timeStep=findAverageTimeStep(infoFileName)
%findAverageTimeStep Find average time step of a DSLM experiment.
%   timeStep=findAverageTimeStep(infoFileName) find and return average time
%   step from an info.txt log file.

% Read in full content of the info file.
text=fileread(infoFileName);

% Find all instances of where scan is finished.
start_inds=strfind(text,' 	Scan Ready for time point: ');

% Find time and index of second time point.
start_vals=readOneInstance(text, start_inds, 2);
% Find time and index of last time point.
end_vals=readOneInstance(text, start_inds, length(start_inds));

% Compute the average time.
timeStep=[3600 60 1 1/1000]*(end_vals(1:4)-start_vals(1:4))/(end_vals(5)-start_vals(5));

end

% Find a time stamp and index number of index number text_ind.
function valuesOfLine=readOneInstance(text, start_inds, text_ind)
% Find index of the last white space character before the time stamp.
start=start_inds(text_ind)-1;
while(1)
    if isstrprop(text(start),'wspace')
        break;
    end
    start=start-1;
end
% Add one to start from the first character of the time stamp.
start=start+1;

% Find index of the first white space character after the time point index.
end_ind=start_inds(text_ind)+29;
while(1)
    if isstrprop(text(end_ind),'wspace')
        break;
    end
    end_ind=end_ind+1;
end
% Subtract one to finish to the last character of the line.
end_ind=end_ind-1;

valuesOfLine=sscanf(text(start:end_ind),'%d:%d:%d.%d 	Scan Ready for time point: %d');
end

