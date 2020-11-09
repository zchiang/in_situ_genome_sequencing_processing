% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Sets the maximum pixel value of an image based on either an absolute or
% percentile based threshold

function[O] = capImage(I,val,type)

    if type == "abs"
        cap = val;
    elseif type == "prc"
        cap = prctile(reshape(I,[],1),val);
    end
    
    O = I;
    O(O>cap) = cap;

end