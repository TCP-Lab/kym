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
% -none-    -> plot    -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% If minimum distance isn't defined for all of x dimensions
% the first value is used as the default for all of the dimensions
if (length(minDist) ~= 2)
	minDist = minDist(ones(1,2));
end

n = size(x,2);
nscale = size(x,1);
nvoice = size(x,1)/(length(freqVec)-1);
lowest = freqVec(end)/(2^(length(freqVec)-1));

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

% Threshold and COI
threshold = thr/100;
wtabst = (wtabs >= threshold);
maxima = maxima .* wtabst .* coimask;

% Find and mark the peaks
[rx,cx] = find(maxima);
index1 = find(maxima);

% Convert to RGB and shade filtered regions gray - If length(index)==0 this does nothing
index2 = find(coimask == 0 | wtabst == 0);
[wtabsind,b] = gray2ind(wtabs,128);
wtabs = ind2rgb(wtabsind,jet(128));
gswtabs = ind2rgb(wtabsind,gray(128));
wtabs(index2) = gswtabs(index2);
wtabs(index2 + n*nscale) = gswtabs(index2 + n*nscale);
wtabs(index2 + 2*n*nscale) = gswtabs(index2 + 2*n*nscale);

% Colour markers blue - If length(mark)==0 this does nothing
wtabs(:,mark,1) = 0;
wtabs(:,mark,2) = 0.5;
wtabs(:,mark,3) = 1;

% Colour maxima projections red
wtabs(rx,:,1) = 1;
wtabs(rx,:,2) = 0;
wtabs(rx,:,3) = 0;

% Colour maxima white
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

% Display results
image(timeVec(1):timeVec(end),1:size(wtabs,1),wtabs), hold on
xlim([timeVec(1),timeVec(end)])
ylim([0.5,size(wtabs,1)+0.5])
set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'YTickLabel',num2str(freqVec'));
ylabel('Frequency (mHz)','FontSize',s2)
xlabel('Time (s)','FontSize',s2)

% Resize -> Sintax Template: set(gca,'Position',[left bottom width height])
set(gca,'Position',[0.12,0.11,0.80,0.815],'Color','none')
axes('Position',[0.12,0.11,0.80,0.815],'XAxisLocation','bottom','YAxisLocation','right','Color','none','YColor','r');
xlim([timeVec(1),timeVec(end)])
ylim([0.5,size(wtabs,1)+0.5])
set(gca,'FontSize',s1,'XTick',[]);
set(gca,'YDir','normal','YTick',unique(rx));

maxVec = lowest*(2.^(unique(rx)/nvoice));
maxVec = floor(maxVec*10)/10;

set(gca,'YTickLabel',num2str(maxVec));


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
%% Neurosciences Ph.D. - Experimental Neurosciences - XXV Cycle                                         %%
%% Department of Life Sciences and Systems Biology                                                      %%
%% Laboratory of Cellular Neurophysiology                                                               %%
%% Via Accademia Albertina 13 10123 Torino                                                              %%
%%                                                                                                      %%
%% Acknowledgements:                                                                                    %%
%% -----------------                                                                                    %%
%% Wavelet Transform computation is here implemented as a product in the Fourier transformed domain.    %%
%% A standard code for this algorithm can be found, for instance, in WaveLab850.                        %%
%% http://www-stat.stanford.edu/~wavelab/                                                               %%
%%                                                                                                      %%
%% Peaks detection uses a technique that is based on images dilation.                                   %%
%% See, for instance, localMaximum.m m-file by Yonathan Nativ.                                          %%
%% http://www.mathworks.com/matlabcentral/fileexchange/authors/26510/                                   %%
%%                                                                                                      %%
%%------------------------------------------------------------------------------------------------------%%
%%------------------------------------------------------------------------------------------------------%%