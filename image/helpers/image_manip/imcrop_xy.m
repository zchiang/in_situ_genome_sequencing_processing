% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Crops an image in 2D using a rectangle (x_start, y_start, x_len, y_len) 
% regardless of the number of dimensions

function[cropped_img] = imcrop_xy(img, rect_2d)

    if ndims(img) == 3
        cropped_img = img(rect_2d(2):rect_2d(2)+rect_2d(4)-1, rect_2d(1):rect_2d(1)+rect_2d(3)-1,:);
    elseif ndims(img) == 4
        cropped_img = img(rect_2d(2):rect_2d(2)+rect_2d(4)-1, rect_2d(1):rect_2d(1)+rect_2d(3)-1,:,:);
    elseif ndims(img) == 2
        cropped_img = img(rect_2d(2):rect_2d(2)+rect_2d(4)-1, rect_2d(1):rect_2d(1)+rect_2d(3)-1);
    elseif ndims(img) == 5
        cropped_img = img(rect_2d(2):rect_2d(2)+rect_2d(4)-1, rect_2d(1):rect_2d(1)+rect_2d(3)-1,:,:,:);
    end
    
end 
