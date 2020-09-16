function[stackCr] = simple_offset_xyz_applied(stackC,offsets)
   
    stackCr = uint16(zeros(size(stackC)));

    xoffset = offsets(1);
    yoffset = offsets(2);
    zoffset = offsets(3);
    
    xbegin_A = max(1+xoffset, 1);
    xend_A   = min(size(stackC,2), size(stackC,2)+xoffset);
    xbegin_B = max(1-xoffset, 1);
    xend_B   = min(size(stackC,2), size(stackC,2)-xoffset);
    
    ybegin_A = max(1+yoffset, 1);
    yend_A   = min(size(stackC,2), size(stackC,2)+yoffset);
    ybegin_B = max(1-yoffset, 1);
    yend_B   = min(size(stackC,2), size(stackC,2)-yoffset);
    
    zbegin_A = max(1+zoffset, 1);
    zend_A   = min(size(stackC,3), size(stackC,3)+zoffset);
    zbegin_B = max(1-zoffset, 1);
    zend_B   = min(size(stackC,3), size(stackC,3)-zoffset);
           
    stackCr(ybegin_A:yend_A,xbegin_A:xend_A,zbegin_A:zend_A,:) = stackC(ybegin_B:yend_B,xbegin_B:xend_B,zbegin_B:zend_B,:);
    stackCr(stackCr == 0) = median(stackC(:));
    
end