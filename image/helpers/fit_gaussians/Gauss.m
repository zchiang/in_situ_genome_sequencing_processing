% Originally used in Buenrostro et al. NBT 2014
% Author: Jason Buenrostro, Greenleaf lab, Stanford University
% Gaussian fit function

function fitImage = Gauss(subFullFit,xMax,yMax,c)
% generate empty grid
fitImage = zeros(xMax, yMax);

% loop through clusters and add to fit image
for i = 1:size(subFullFit,1);
    amp = subFullFit(i,1);
    xCent = subFullFit(i,3)-c;
    yCent = subFullFit(i,2)-c;
    sigma = subFullFit(i,4);
    
    % if not 0's add to image
    if sum(subFullFit(i,1)) ~= 0;
        for x = 1:xMax;
            for y = 1:yMax;
                xComp = ((x-xCent)^2)/(2*(sigma^2));
                yComp = ((y-yCent)^2)/(2*(sigma^2));
                fitImage(x,y) = fitImage(x,y)+amp*exp(-1*(xComp+yComp));
            end
        end
    end
end
end