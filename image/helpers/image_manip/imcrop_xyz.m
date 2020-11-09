% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Crops an image in 3D using a rectangle (x_start, y_start, z_start, x_len,
% y_len. z_len) regardless of the number of dimensions

function[cropped_img] = imcrop_xyz(img, rect_3d)

    if ndims(img) == 3
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(5)-1, rect_3d(1):rect_3d(1)+rect_3d(4)-1, rect_3d(3):rect_3d(3)+rect_3d(6)-1);
    elseif ndims(img) == 4
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(5)-1, rect_3d(1):rect_3d(1)+rect_3d(4)-1, rect_3d(3):rect_3d(3)+rect_3d(6)-1,:);
    elseif ndims(img) == 5
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(5)-1, rect_3d(1):rect_3d(1)+rect_3d(4)-1, rect_3d(3):rect_3d(3)+rect_3d(6)-1,:,:);
    end
        
end 

