% Originally used in Buenrostro et al. NBT 2014
% Author: Jason Buenrostro, Greenleaf lab, Stanford University
% Uses lsqcurvefit to fit 2D Gaussians over sub-images

function FullFit = QuantifyImage(I,adjX,adjY,zeroPix,adjVals)
% Background and cluster centers
bgEst = mode(I(:));

% Define subtile size
stepSize = 24;  % Stepsize (< tile size) was 16
if strcmp(adjVals,'allCluster') == 1
     bound = 6;
else
     bound = 4;
end
subsetSize = stepSize+(bound*2);  % subset
startX = bound; startY = bound;  % start val
[ySize, xSize] = size(I);  % image size
endX = xSize - stepSize - bound; endY = ySize - stepSize - bound;  % end

% Make X Y coordinates for sub region (as column vectors)
[X,Y] = meshgrid(1:(subsetSize),1:(subsetSize)); %y x-y coordinates
x(1:(subsetSize*subsetSize),1) = X(:); % x = first column
x(1:(subsetSize*subsetSize),2) = Y(:); % y = second column

%define function
fun = @N2DGauss_JACOBIAN;

% initialize output
FullFit = zeros(length(adjX),5);

% Generate matrix to store fit results
if strcmp(adjVals,'allCluster') == 1
    % read seq input and extract x/y positions
    FullFit(:, 2) = adjX;
    FullFit(:, 3) = adjY;
    
    % define fit function options
    options = optimset('TolX',0.05, 'TolFun', 0.0001,'MaxFunEvals', 200,...
    'Display', 'off', 'jacobian', 'on', 'DerivativeCheck', 'off',...
    'PrecondBandWidth', 0, 'LargeScale', 'on');
else
    % read seq input and extract x/y positions
    fid = fopen(adjVals);
    C = textscan(fid, '%f %f %f %f %f');
    xVals = C{2}; yVals = C{3}; sigVals = C{4};
    fclose(fid);
    
    % define fit function options
    options = optimset('TolX',0.05, 'TolFun', 0.0001,'MaxFunEvals', 200,...
    'Display', 'off', 'jacobian', 'on', 'DerivativeCheck', 'off',...
    'PrecondBandWidth', 0, 'LargeScale', 'on');
end

%%% UNCOMMENT FOR DEBUGGING  %%%
%startX = 900; startY = 900;  % start position
%endX = startX+200; endY = startY+200;  % end position

%tic
% Loop over subimages
for xLoc = startX:stepSize:endX
    for yLoc = startY:stepSize:endY
        % clear previous
        clear c0 ub lb cc xEst yEst
        
        % subimage dimensions
        minX = xLoc-bound+1; maxX = xLoc+stepSize+bound;
        minY = yLoc-bound+1; maxY = yLoc+stepSize+bound;
        
        % subimage
        Isub = I(minY:maxY,minX:maxX);
        
        % get cluster centers within fit window
        fIdx = find(adjY >= minY & adjX >= minX & adjX <= maxX & adjY <= maxY);
        n = length(fIdx);
        aEst = zeros(1,n); 
        
        % if frame overlaps with junk then pass
        cPos = sub2ind(size(I),round(adjY(fIdx)),round(adjX(fIdx)));
        idx2 = ismember(cPos,zeroPix);

        % if clusters
        if n && sum(idx2)==0
            % initialize arrays
            nsz = 4*n+3;
            c0 = zeros(1,nsz);
            ub = zeros(1,nsz);
            lb = zeros(1,nsz);
            
            % get first guess amplitude
            for i = 1:n
                a = round(adjY(fIdx(i)));
                b = round(adjX(fIdx(i)));
                aEst(i) = I(a,b);
            end
            aEst = aEst - bgEst;
            aEst(aEst<=0)=1;
            
            % first guess x/y position
            xEst = adjX(fIdx) - xLoc + bound;
            yEst = adjY(fIdx) - yLoc + bound;
            
            % define common fit parameters
            ampPos = 1:4:nsz-3; xCentPos = 2:4:nsz-3; yCentPos = 3:4:nsz-3; sigPos = 4:4:nsz-3;
            %ub(ampPos) = 3*aEst + bgEst; lb(ampPos) = 0;  % amplitude
            ub(ampPos) = 3*aEst; lb(ampPos) = aEst;  % amplitude
            c0(end-2) = 0.0; ub(end-2) = 2.0; lb(end-2) = -2.0;  % x offset
            c0(end-1) = 0.0; ub(end-1) = 2.0; lb(end-1) = -2.0;  % y offset
            c0(end) = bgEst; ub(end) = 2*bgEst; lb(end) = 0;  % background
            c0(sigPos) = 1.5; ub(sigPos) = 3.0; lb(sigPos) = 0;  % sigma
	    
            % define specific fit parameters
            if strcmp(adjVals,'allCluster') == 1
                c0(ampPos) = aEst;
                c0(xCentPos) = xEst; ub(xCentPos) = xEst + 1.0; lb(xCentPos) = xEst - 1.0;  % x center
                c0(yCentPos) = yEst; ub(yCentPos) = yEst + 1.0; lb(yCentPos) = yEst - 1.0;  % y center
            else
                c0(ampPos) = bgEst;
                c0(sigPos) = sigVals(fIdx);  % sigma
                sig0idx = sigPos(c0(sigPos)==0); c0(sig0idx) = 1.3;  %sigma
                c0(xCentPos) = xEst; ub(xCentPos) = xEst + 0.5; lb(xCentPos) = xEst - 0.5;  % x center
                c0(yCentPos) = yEst; ub(yCentPos) = yEst + 0.5; lb(yCentPos) = yEst - 0.5;  % y center
            end
            
            % fit
            cc = lsqcurvefit(fun,c0,x,Isub(:),lb,ub,options);
            
            % get cluster centers within report window
            xCent = adjX(fIdx);
            yCent = adjY(fIdx);
            minxR = xLoc; maxxR = xLoc+stepSize;
            minyR = yLoc; maxyR = yLoc+stepSize;
            rIdx = find(xCent >= minxR & yCent >= minyR & xCent <= maxxR & yCent <= maxyR);
            
            % output to FullFit
            FullFit(fIdx(rIdx), 1) = cc(ampPos(rIdx));  % amp
            FullFit(fIdx(rIdx), 2) = cc(xCentPos(rIdx)) + xLoc + cc(end-2) - bound;  % xCent
            FullFit(fIdx(rIdx), 3) = cc(yCentPos(rIdx)) + yLoc + cc(end-1) - bound;  % yCent
            FullFit(fIdx(rIdx), 4) = cc(sigPos(rIdx));  % sigma
            FullFit(fIdx(rIdx), 5) = cc(end);  % background
        end
    end
    %toc
end
%toc
end

%%% Jacobian for 2D guassian fit %%%
function [F,J] = N2DGauss_JACOBIAN(c, xyArray)
%Generates a function of N gaussians where N is (lenght(c)-1)/5)
F = 0;
J = zeros(length(xyArray), length(c));

% loop through clusters
for i = 0:4:length(c)-4
    % get values
    amp = c(i+1);
    xCluster = c(i+2);yCluster = c(i+3);
    xGlobal = c(end-2);yGlobal = c(end-1);
    sigma = c(i+4);

    % calculate x/y components
    xPos = xyArray(:,1)-xCluster-xGlobal; yPos = xyArray(:,2)-yCluster-yGlobal;
    xComp = xPos.^2/(2*sigma.^2); yComp = yPos.^2/(2*sigma.^2);
    
    % add to fit
    result = amp*exp(-1*(xComp+yComp));
    F = F + result;
    
    % jacobians
    if nargout > 1
        J(:, i+1) = result./amp;  % amplitude
        J(:, i+2) = result.*xPos./(sigma.^2);  % local x
        J(:, i+3) = result.*yPos./(sigma.^2);  % local y
        J(:, i+4) = result.*(yPos.^2 + xPos.^2)./sigma.^3;  % sigma
    end
    J(:,end-2) = J(:,end-2) + result.*xPos./(sigma.^2);
    J(:,end-1) = J(:,end-1) + result.*yPos./(sigma.^2);
end

% Add background
F = F + c(end);
if nargout > 1
    J(:, end)= 1;
end
end
