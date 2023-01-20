function [maxmap1,maxmap2,maxima,maximaNT] = paths(x,y,coimask,thr,mark,minDist,noPlateau)

%
%--------------------------------------------------------------------------------
% CWT Frequency Paths
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [maxmap1,maxmap2,maxima,maximaNT] = paths(x,y,coimask,thr,mark,minDist,noPlateau)
%
% INPUT        TYPE       MEANING
% -----        ----       -------
% x         -> matrix  -> Continuous Wavelet Modulus
% y         -> matrix  -> Continuous Wavelet Real Part
% coimask   -> matrix  -> COI Mask
% thr       -> scalar  -> Threshold Percentage
% mark      -> array   -> Time Marker Set (Step Number)
% minDist   -> scalar  -> Minimum Distance Between 2 Peaks
% noPlateau -> boolean -> Recognize Points with the Same Value as Peaks
%
% OUTPUT       TYPE        MEANING
% ------       ----        -------
% maxmap1   -> matrix   -> Continuous Wavelet Modulus with Paths
% maxmap2   -> matrix   -> Continuous Wavelet Real Part with Paths
% maxima    -> matrix   -> Thresholded Maxima Coordinate
% maximaNT  -> matrix   -> Non-Thresholded Maxima Coordinate
%

n = size(x)(2);
nscale = size(x)(1);

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

% North-South Directional Maxima
maxima = zeros(nscale,n);
minDist = minDist(1);
se = ones(minDist,1);
for k = 1:n
	X = imdilate(x(:,k),se);
	maxima(:,k) = (x(:,k) == X);
endfor

% x Threshold
minimox = min(min(x));
x = x - minimox;
massimox = max(max(x));
x = x / massimox;

threshold = thr/100;

xt = (x >= threshold);
maximaNT = maxima .* coimask;
maxima = maxima .* xt .* coimask;

% Merge with original scalogram and Mark the peaks in red
index = find(maxima);
[rx,cx] = find(maxima);

maxmap1 = (x .* xt .* coimask);
[maxmap1,b] = gray2ind(maxmap1);
maxmap1 = ind2rgb(maxmap1,bone());

%maxmap1(index) = 0.85; % Constant Path Intensity
maxmap1(index) = 0.7 * maxmap1(index) + 0.3; % Alternative Color Gradient
maxmap1(index + n*nscale) = maxmap1(index + 2*n*nscale) = 0;

for  h = 1:length(index) 
	if (mod(index(h),nscale) != 1)
		%maxmap1(index(h)-1) = 0.55;
		maxmap1(index(h)-1) = 0.8 * maxmap1(index(h)-1) + 0.2;
		maxmap1(index(h)-1 + n*nscale) = maxmap1(index(h)-1 + 2*n*nscale) = 0;
	endif
	if (mod(index(h),nscale) != 0)
		%maxmap1(index(h)+1) = 0.55;
		maxmap1(index(h)+1) = 0.8 * maxmap1(index(h)+1) + 0.2;
		maxmap1(index(h)+1 + n*nscale) = maxmap1(index(h)+1 + 2*n*nscale) = 0;
	endif
endfor

% Markers
maxmap1(:,mark,1) = 0;
maxmap1(:,mark,2) = 0.85;
maxmap1(:,mark,3) = 0.85;

% y Threshold
minimoy = min(min(y));
y = y - minimoy;
massimoy = max(max(y));
y = y / massimoy;

yt = (y >= threshold);
%yt = xt; %Per usare gli stessi perimetri di soglia del modulo

% Merge with original scalogram and Mark the peaks in red
maxmap2 = (y .* yt .* coimask);
[maxmap2,b] = gray2ind(maxmap2);
maxmap2 = ind2rgb(maxmap2,bone());

%maxmap2(index) = 0.85; % Constant Path Intensity
maxmap2(index) = 0.7 * maxmap2(index) + 0.3; % Alternative Color Gradient
maxmap2(index + n*nscale) = maxmap2(index + 2*n*nscale) = 0;

for  h = 1:length(index) 
	if (mod(index(h),nscale) != 1)
		%maxmap2(index(h)-1) = 0.55;
		maxmap2(index(h)-1) = 0.8 * maxmap2(index(h)-1) + 0.2;
		maxmap2(index(h)-1 + n*nscale) = maxmap2(index(h)-1 + 2*n*nscale) = 0;
	endif
	if (mod(index(h),nscale) != 0)
		%maxmap2(index(h)+1) = 0.55;
		maxmap2(index(h)+1) = 0.8 * maxmap2(index(h)+1) + 0.2;
		maxmap2(index(h)+1 + n*nscale) = maxmap2(index(h)+1 + 2*n*nscale) = 0;
	endif
endfor

% Markers
maxmap2(:,mark,1) = 0;
maxmap2(:,mark,2) = 0.85;
maxmap2(:,mark,3) = 0.85;


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