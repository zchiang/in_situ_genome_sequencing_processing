function chr = getTrainingData(dataname)

%List of human chromosome lengths
chrLen = [249698942 242508799 198450956 190424264 181630948 170805979 ...
    159345973 145138636 138688728 133797422 135186938 133275309 114364328 ...
    108136338 102439437 92211104 83836422 80373285 58617616 64444167 ...
    46709983 51857516 156040895 57264655];

%Read the training data
data = readtable(dataname);
nFOV = max(data.fov);
nCell = max(data.cell);
nChr = max(data.hg38_chr);
nCluster = max(data.cluster);


%Loop through chromosomes
for kk = 1:nChr
    
    distP_all = [];
    distG_all = [];
    
    for ii = 1:nFOV
        for jj = 1:nCell
            for ll = 1:nCluster
                %Collect genomic and physical distances
                [distP,distG] = genomicAndPhysicalDist(data, ii, jj, kk, ll);
                
                distP_all = [distP_all; distP];
                distG_all = [distG_all; distG];
            end
        end
    end
    
    chr(kk).distP = distP_all;
    chr(kk).distG = distG_all;
    chr(kk).length = chrLen(kk);
    

end


