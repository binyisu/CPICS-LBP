function [result,lbpimage] = icslbp(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.3
% Authors: binyi su

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a
% predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.


% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};
d_image=double(image);

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);

    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
end

% Determine the dimensions of the input image.
[ysize xsize] = size(image);



miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1; %3
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1; %3

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0)); %2
origx=1-floor(min(minx,0)); %2

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex; %22
dy = ysize - bsizey; %22

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx); %2:24,2:24
d_C = double(C); 

% bins = 2^5;
bins = 2^(neighbors/2+1);  %256

% Initialize the result matrix with zeros.
result=zeros(dy+1,dx+1); %23*23
sumN = zeros(dy+1,dx+1);
%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N(1:dy+1,1:dx+1,i) = image(ry:ry+dy,rx:rx+dx);
    D = N(1:dy+1,1:dx+1,i) >= C; 
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = roundn((1 - tx) * (1 - ty),-6);
    w2 = roundn(tx * (1 - ty),-6);
    w3 = roundn((1 - tx) * ty,-6) ;
    % w4 = roundn(tx * ty,-6) ;
    w4 = roundn(1 - w1 - w2 - w3, -6);
            
    % Compute interpolated pixel values
    N(1:dy+1,1:dx+1,i) = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    N(1:dy+1,1:dx+1,i) = roundn(N(1:dy+1,1:dx+1,i),-4);
%     D = N >= d_C; 
  end 
  % Update the result matrix.
%   v = 2^(i-1);
%   result = result + v*D;
end
for k = 1:neighbors
    sumN = sumN + N(:,:,i);
end
averageN = sumN/neighbors;
for j = 1:neighbors/2
    D = N(:,:,j)>=N(:,:,j+(neighbors/2));
    v = 2^j;
    result = result + v*D;
end
image_o = d_image(fy:fy+dy,fx:fx+dx);
result = result + (averageN > image_o)*2^0;
%   result = result + ((N(1:dy+1,1:dx+1,1)+N(1:dy+1,1:dx+1,5))/2<d_C)*2^0+((N(1:dy+1,1:dx+1,2)+N(1:dy+1,1:dx+1,6))/2<d_C)*2^1+...
%       ((N(1:dy+1,1:dx+1,3)+N(1:dy+1,1:dx+1,7))/2<d_C)*2^2+((N(1:dy+1,1:dx+1,4)+N(1:dy+1,1:dx+1,8))/2<d_C)*2^3;
%% 原始ICS-LBP
%     resultaverge = (N(1:dy+1,1:dx+1,1)+N(1:dy+1,1:dx+1,2)+N(1:dy+1,1:dx+1,3)+N(1:dy+1,1:dx+1,4)+N(1:dy+1,1:dx+1,5)+...
%         N(1:dy+1,1:dx+1,6)+N(1:dy+1,1:dx+1,7)+N(1:dy+1,1:dx+1,8))/8;
%         image_o = d_image(fy:fy+dy,fx:fx+dx);
%     result = result + (resultaverge(1:dy+1,1:dx+1)>image_o(1:dy+1,1:dx+1))*2^0+(N(1:dy+1,1:dx+1,1)>N(1:dy+1,1:dx+1,5))*2^1+(N(1:dy+1,1:dx+1,2)>N(1:dy+1,1:dx+1,6))*2^2+...
%       (N(1:dy+1,1:dx+1,3)>N(1:dy+1,1:dx+1,7))*2^3+(N(1:dy+1,1:dx+1,4)>N(1:dy+1,1:dx+1,8))*2^4;

  lbpimage = result; %126*126原图（128*128）
%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(result,1) %23
        for j = 1:size(result,2) %23
            result(i,j) = mapping.table(result(i,j)+1);  %mapping之后的结果
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    result=hist(result(:),0:(bins-1));
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
else
    %Otherwise return a matrix of unsigned integers
    if ((bins-1)<=intmax('uint8'))
        result=uint8(result);
    elseif ((bins-1)<=intmax('uint16'))
        result=uint16(result);
    else
        result=uint32(result);
    end
end

end

function x = roundn(x, n)

error(nargchk(2, 2, nargin, 'struct'))
validateattributes(x, {'single', 'double'}, {}, 'ROUNDN', 'X')
validateattributes(n, ...
    {'numeric'}, {'scalar', 'real', 'integer'}, 'ROUNDN', 'N')

if n < 0
    p = 10 ^ -n;
    x = round(p * x) / p;
elseif n > 0
    p = 10 ^ n;
    x = p * round(x / p);
else
    x = round(x);
end

end


