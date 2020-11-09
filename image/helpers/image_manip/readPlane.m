% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Uses BioFormats bfGetPlane to get an image given a reader, series,
% z-slice, channel, and timepoint

function[img] = readPlane(reader, series, z, c, t)

    reader.setSeries(series - 1);
    plane = reader.getIndex(z - 1, c -1, t - 1) + 1;
    img = bfGetPlane(reader, plane);

end