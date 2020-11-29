% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Uses cross correlation to register one stack to another in x, y, and z

function[stackBr] = simple_offset_xyz(stackA, stackB)
    
    % handle stacks with different numbers of z-planes

    if size(stackA,3) < size(stackB,3)
        new_stackA = zeros(size(stackB),'uint16');
        new_stackA(:,:,1:size(stackA,3)) = stackA;
        stackA = new_stackA;
        stackA(stackA == 0) = median(stackA(:));
        clearvars new_stackA
    end

    stackBr = zeros(size(stackB),'uint16');

    % take maximum intensity projections
    
    imgA = squeeze(max(max(stackA,[],3),[],4));
    imgB = squeeze(max(max(stackB,[],3),[],4));
    
    imgC = squeeze(max(max(stackA,[],2),[],4));
    imgD = squeeze(max(max(stackB,[],2),[],4));
    
	% perform cross correlation
    
    c = normxcorr2(imgB, imgA);
    [max_c, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(imgB,1)) (ypeak-size(imgB,2))];
    xoffset = corr_offset(1);
    yoffset = corr_offset(2);
    
    c = normxcorr2(imgD, imgC);
    [max_c, imax] = max(abs(c(:)));
    [xpeak, zpeak] = ind2sub(size(c),imax(1));
    corr_offset = [(xpeak-size(imgD,1)) (zpeak-size(imgD,2))];
    zoffset = corr_offset(2);
    
	% calculate offset bounds
    
    xbegin_A = max(1+xoffset, 1);
    xend_A   = min(size(stackB,2), size(stackB,2)+xoffset);
    xbegin_B = max(1-xoffset, 1);
    xend_B   = min(size(stackB,2), size(stackB,2)-xoffset);
    
    ybegin_A = max(1+yoffset, 1);
    yend_A   = min(size(stackB,2), size(stackB,2)+yoffset);
    ybegin_B = max(1-yoffset, 1);
    yend_B   = min(size(stackB,2), size(stackB,2)-yoffset);
    
    zbegin_A = max(1+zoffset, 1);
    zend_A   = min(size(stackB,3), size(stackB,3)+zoffset);
    zbegin_B = max(1-zoffset, 1);
    zend_B   = min(size(stackB,3), size(stackB,3)-zoffset);

    % apply offsets    
    
    stackBr(ybegin_A:yend_A,xbegin_A:xend_A,zbegin_A:zend_A,:) = stackB(ybegin_B:yend_B,xbegin_B:xend_B,zbegin_B:zend_B,:);
    stackBc(stackBr == 0) = median(stackB(:));
    
end