function peak(x,coimask,thr,unit,US,timeVec,freqVec,mark,minDist,noPlateau)

%
%--------------------------------------------------------------------------------
% NU CWT Peaks Detection
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% peak(x,coimask,thr,unit,US,timeVec,freqVec,mark,minDist,noPlateau)
%
% INPUT        TYPE       MEANING
% -----        ----       -------
% x         -> matrix  -> Continuous Wavelet Modulus
% coimask   -> matrix  -> COI Mask
% thr       -> scalar  -> Threshold Percentage
% unit      -> string  -> Time Unit: 's' or 'ms'
% US        -> scalar  -> US-Upsampling factor
% timeVec   -> array   -> Time Vector
% freqVec   -> array   -> Frequency Vector
% mark      -> array   -> Time Marker Set (Step Number)
% minDist   -> array   -> Minimum Distance Between 2 Peaks
% noPlateau -> boolean -> Recognize Points with the Same Value as Peaks
%
% OUTPUT       TYPE       MEANING
% ------       ----       -------
% -none-    -> plot    -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
[nscale,n] = size(x);
nvoice = nscale/(length(freqVec)-1);
lowest = freqVec(end)/(2^(length(freqVec)-1));

% If minimum distance isn't defined for all of x dimensions
% the first value is used as the default for all of the dimensions
if (length(minDist) ~= 2)
	minDist = minDist(ones(1,2));
end

% Plateau handling
% Without this code points with the same hight will be recognized as peaks
if (noPlateau)
	temp = sort(x(:));
	dY = diff(temp);
	% Find the minimum step in the data
	minimumDiff = min(dY(dY ~= 0));
	% Add noise which won't affect the peaks
	x = x + rand(size(x))*minimumDiff;
end

% Peaks detection
se = ones(minDist);
X = imdilate(x,se);
maxima = (x == X);

% Rescale to [0,1]
wtabs = (x-min(min(x)))/max(max(x-min(min(x))));

% Threshold, COI and "ultra-Nyquist zone"
threshold = thr/100;
wtabst = (wtabs >= threshold);
maxima = maxima .* wtabst .* coimask .* [ones(nscale-log2(US)*nvoice,n);zeros(log2(US)*nvoice,n)];

% Find peaks
[rx,cx] = find(maxima);
index1 = find(maxima);

% Convert to RGB
[wtabsind,b] = gray2ind(wtabs,128);
wtabs = ind2rgb(wtabsind,jet(128));
gswtabs = ind2rgb(wtabsind,gray(128));

% Shade thresholded regions and the "ultra-Nyquist zone" with gray
index2 = find(coimask == 0 | wtabst == 0);
wtabs(index2) = gswtabs(index2);
wtabs(index2 + n*nscale) = gswtabs(index2 + n*nscale);
wtabs(index2 + 2*n*nscale) = gswtabs(index2 + 2*n*nscale);
wtabs(end-log2(US)*nvoice+1:end,:,:) = gswtabs(end-log2(US)*nvoice+1:end,:,:);

% Blue Nyquist limit
wtabs(end-log2(US)*nvoice,:,1) = 0;
wtabs(end-log2(US)*nvoice,:,2) = 0.5;
wtabs(end-log2(US)*nvoice,:,3) = 1;

% Blue markers - If length(mark)==0 this does nothing
wtabs(:,mark,1) = 0;
wtabs(:,mark,2) = 0.5;
wtabs(:,mark,3) = 1;

% Red maxima projections
wtabs(rx,:,1) = 1;
wtabs(rx,:,2) = 0;
wtabs(rx,:,3) = 0;

% White maxima
wtabs(index1) = 1;
wtabs(index1 + n*nscale) = 1;
wtabs(index1 + 2*n*nscale) = 1;

% Increase maxima thickness
indexx = index1(find(index1 > 3*nscale));
for k = 1:3
	wtabs(indexx - k*nscale) = 1;
	wtabs(indexx - k*nscale + n*nscale) = 1;
	wtabs(indexx - k*nscale + 2*n*nscale) = 1;
end
indexx = index1(find(index1 <= (n*nscale)-3*nscale));
for k = 1:3
	wtabs(indexx + k*nscale) = 1;
	wtabs(indexx + k*nscale + n*nscale) = 1;
	wtabs(indexx + k*nscale + 2*n*nscale) = 1;
end

% Display scalogram
image(timeVec,1:nscale,wtabs), hold on
xlim([timeVec(1),timeVec(end)])
ylim([0.5,nscale+0.5])
set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'YTickLabel',num2str(freqVec',4));
if (strcmp(unit,'ms'))
	ylabel('Frequency (Hz)','FontSize',s2)
else
	ylabel('Frequency (mHz)','FontSize',s2)
end
xlabel('Time (s)','FontSize',s2)

% Label maxima
set(gca,'Position',[0.12,0.11,0.80,0.815],'Color','none')
axes('Position',[0.12,0.11,0.80,0.815],'XAxisLocation','bottom','YAxisLocation','right','Color','none','YColor','r');
xlim([timeVec(1),timeVec(end)])
ylim([0.5,nscale+0.5])
set(gca,'FontSize',s1,'XTick',[]);
set(gca,'YDir','normal','YTick',unique(rx));
maxVec = lowest*(2.^(unique(rx)/nvoice));
set(gca,'YTickLabel',num2str(maxVec,4));


%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%
%%                                                                                                      %%
%% KYM Project                                                                                          %%
%% -----------                                                                                          %%
%% First Released in 2010                                                                               %%
%% Original code by Federico Alessandro Ruffinatti                                                      %%
%%                                                                                                      %%
%% UNIVERSITY OF TORINO                                                                                 %%
%% DOCTORAL SCHOOL IN LIFE AND HEALTH SCIENCES                                                          %%
%% Neurosciences PhD - Experimental Neurosciences - XXV Cycle                                           %%
%% Department of Life Sciences and Systems Biology                                                      %%
%% Laboratory of Cellular Neurophysiology                                                               %%
%% Via Accademia Albertina 13 10123 Torino                                                              %%
%%                                                                                                      %%
%% Acknowledgements:                                                                                    %%
%% -----------------                                                                                    %%
%% CWT convolution is implemented as a product in the Fourier transformed domain.                       %%
%% In particular, the code for CWT computation is a refinement of WaveLab850 dyadic algorithm.          %%
%% http://www-stat.stanford.edu/~wavelab/                                                               %%
%%                                                                                                      %%
%% A technique based on image dilation has been used for the detection of peaks and maxima.             %%
%% This idea comes from Yonathan Nativ's localMaximum.m m-file.                                         %%
%% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/                                   %%
%%                                                                                                      %%
%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%