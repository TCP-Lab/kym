function peak(x,coimask,thr,timeVec,freqVec,mark,minDist,noPlateau)

%
%--------------------------------------------------------------------------------
% CWT Peaks Detection
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% peak(x,coimask,thr,timeVec,freqVec,mark,minDist,noPlateau)
%
% INPUT        TYPE       MEANING
% -----        ----       -------
% x         -> matrix  -> Continuous Wavelet Modulus
% coimask   -> matrix  -> COI Mask
% thr       -> scalar  -> Threshold Percentage
% timeVec   -> array   -> Time Vector
% freqVec   -> array   -> Frequency Vector
% mark      -> array   -> Time Marker Set (Step Number)
% minDist   -> array   -> Minimum Distance Between 2 Peaks
% noPlateau -> boolean -> Recognize Points with the Same Value as Peaks
%
% OUTPUT       TYPE       MEANING
% ------       ----       -------
% -none-    -> plot    -> Plot Resulting from Analysis
%

% In case minimum distance isn't defined for all of x dimensions
% the first value is used as the default for all of the dimensions
dimX = length(size(x));
if (length(minDist) != dimX)
	minDist = minDist(ones(dimX,1));
endif
% Validity checks
minDist = ceil(minDist);
minDist = max([minDist(:)' ; ones(1,length(minDist))]);
minDist = min([minDist ; size(x)]);

n = size(x)(2);
nscale = size(x)(1);
nvoice = (size(x)(1))/(length(freqVec)-1);
lowest = freqVec(end)/(2**(length(freqVec)-1));

% Plateau Handling
% Without this code points with the same hight will be recognized as peaks
if (noPlateau)
	temp = sort(x(:));
	dY = diff(temp);
	% Finding the minimum step in the data
	minimumDiff = min( dY(dY != 0) );
	% Adding noise which won't affect the peaks
	x = x + rand(size(x))*minimumDiff;
endif

% Peaks Detection
se = ones(minDist);
X = imdilate(x,se);
maxima = (x == X);

% Threshold
minimox = min(min(x));
x = x - minimox;
massimox = max(max(x));
x = x / massimox;

threshold = thr/100;

xt = (x >= threshold);
maxima = maxima .* xt .* coimask;

% Merge with original scalogram and mark the peaks
[rx,cx] = find(maxima);
index = find(maxima);

maxmap = (x .* xt .* coimask);
maxmap(rx,:) = ones(length(rx),n);

maxmap(index) = 0.5;
indexx = index(find(index > 3*nscale));
for k = 1:3
	maxmap(indexx-k*nscale) = 0.5;
endfor
indexx = index(find(index <= (n*nscale)-3*nscale));
for k = 1:3
	maxmap(indexx+k*nscale) = 0.5;
endfor

% Markers
maxmap(:,mark) = 0.5;

% Display Results
imagesc(timeVec,[],maxmap)
set(gca,'XTick',[get(gca,'XTick'),floor(timeVec(end))]);
set(gca,'YDir','normal');
set(gca,'YTick',[0:nvoice:nscale]); % Maximum number of displayable ticks on ordinates
set(gca,'YTick',[get(gca,'YTick'),1]); % In order to display freqVec(1) value
set(gca,'YTickLabel',freqVec);
ylabel('Frequency (mHz)','FontSize',18)
xlabel('Time (s)','FontSize',18)

axes('YAxisLocation','right','YColor','r');
set(gca,'XTick',[]);
ylim([0,nscale]);
set(gca,'YTick',rx);

maxVec = lowest*(2.**(rx/nvoice));
maxVec = floor(maxVec*10)/10;
set(gca,'YTickLabel',maxVec);


%---------------------------------------------------------------------%
%                                                                     %
% A.A. 2009/2010 - 2010/2011                                          %
% Original code by Federico Alessandro Ruffinatti                     %
% Università degli Studi di Torino - Italy - DBAU - Scienze MFN       %
% Scuola di Dottorato in Neuroscienze - XXV ciclo                     %
%                                                                     %
% Wavelet computation is regarded as a time convolution and it is     %
% implemented as a product in the Fourier transformed domain.         %
% A standard code for this algorithm can be found, for instance,      %
% in WaveLab850 - http://www-stat.stanford.edu/~wavelab/              %
%                                                                     %
% Peaks detection uses a technique that is based on images dilation.  %
% See, for instance, localMaximum.m m-file by Yonathan Nativ          %
% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/  %
%                                                                     %
%---------------------------------------------------------------------%