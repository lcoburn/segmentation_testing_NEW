function comb=mergeColourChannels(red,green,blue,grey)
%mergeColourChannels merge RGB colour channels with grey scale image.
%   comb=mergeColourChannels(red,green,blue,grey) red, green, blue and grey
%   are single channel images that are combined to comb that is ouput.

% Grey.
if ~isempty(grey)
    comb1=grey;
    comb2=grey;
    comb3=grey;
else
    comb1=zeros(max([size(red);size(green);size(blue)]));
    comb2=comb1;
    comb3=comb1;    
end

% Red.
if ~isempty(red)
    comb1=comb1+red;
end

% Green
if ~isempty(green)
    comb2=comb2+green;
end

% Blue.
if ~isempty(blue)
    comb3=comb3+blue;
end

% Combine into as output image.
comb=zeros([size(comb1) 3]);
comb(:,:,1)=comb1;
comb(:,:,2)=comb2;
comb(:,:,3)=comb3;

% If all input data types are uint8, use uint8 for output as well.
if (isa(red,'uint8') | isempty(red)) & (isa(green,'uint8') | isempty(green)) ...
        & (isa(blue,'uint8') | isempty(blue)) & (isa(grey,'uint8') | isempty(grey))
    comb=uint8(comb);
end
