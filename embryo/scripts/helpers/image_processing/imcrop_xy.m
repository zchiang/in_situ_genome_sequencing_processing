function[cropped_img] = imcrop_xy(img, rect_3d)

    if ndims(img) == 3
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(4)-1, rect_3d(1):rect_3d(1)+rect_3d(3)-1,:);
    elseif ndims(img) == 4
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(4)-1, rect_3d(1):rect_3d(1)+rect_3d(3)-1,:,:);
    elseif ndims(img) == 2
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(4)-1, rect_3d(1):rect_3d(1)+rect_3d(3)-1);
    elseif ndims(img) == 5
        cropped_img = img(rect_3d(2):rect_3d(2)+rect_3d(4)-1, rect_3d(1):rect_3d(1)+rect_3d(3)-1,:,:,:);
    end
end 
