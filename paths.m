function [maxima,maximaNT] = paths(x,y,coimask,thr,timeVec,freqVec,mark,minDist,noPlateau,grafunc)

%
%--------------------------------------------------------------------------------
% CWT Frequency Paths
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [maxima,maximaNT] = paths(x,y,coimask,thr,timeVec,freqVec,mark,minDist,noPlateau,grafunc)
%
% INPUT        TYPE       MEANING
% -----        ----       -------
% x         -> matrix  -> Continuous Wavelet Modulus
% y         -> matrix  -> Continuous Wavelet Real Part
% coimask   -> matrix  -> COI Mask
% thr       -> scalar  -> Threshold Percentage
% timeVec   -> array   -> Time Vector
% freqVec   -> array   -> Frequency Vector
% mark      -> array   -> Time Marker Set (Step Number)
% minDist   -> scalar  -> Minimum Distance Between 2 Peaks
% noPlateau -> boolean -> Recognize Points with the Same Value as Peaks
% grafunc   -> boolean -> Enable Graphical Functions
%
% OUTPUT       TYPE        MEANING
% ------       ----        -------
% maxima    -> matrix   -> Thresholded Maxima Map
% maximaNT  -> matrix   -> Non-Thresholded Maxima Map
% -none-    -> plot     -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

n = size(x,2);
nscale = size(x,1);
nvoice = size(x,1)/(length(freqVec)-1);

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

% North-South directional maxima
se = ones(minDist,1);
X = imdilate(x,se);
maxima = (x == X);

% Rescale to [0,1]
wtabs = (x-min(min(x)))/max(max(x-min(min(x))));
wtreal = (y-min(min(y)))/max(max(y-min(min(y))));

% Threshold and COI
threshold = thr/100;
wtabst = (wtabs >= threshold);
wtrealt = (wtreal >= threshold);
%wtrealt = wtabst; % To use the same threshold contours used for modulus
maximaNT = maxima .* coimask;
maxima = maxima .* wtabst .* coimask;

% Enable graphical functions
if (grafunc)
	
	% Find maxima
	%[rx,cx] = find(maxima);
	index1 = find(maxima);
	
	%---------------------------------------------------------------------------
	% Modulus plot preparation
	%---------------------------------------------------------------------------
	
	% Convert to RGB and colour filtered regions black - If length(index)==0 this does nothing
	index2 = find(coimask == 0 | wtabst == 0);
	[wtabs,b] = gray2ind(wtabs,128);
	wtabs = ind2rgb(wtabs,bone(128));
	wtabs(index2) = 0;
	wtabs(index2 + n*nscale) = 0;
	wtabs(index2 + 2*n*nscale) = 0;
	
	% Colour paths red
	wtabs(index1) = 0.7*wtabs(index1) + 0.3; % Transparency effect
	%wtabs(index1) = 0.85; % Constant path intensity
	wtabs(index1 + n*nscale) = 0;
	wtabs(index1 + 2*n*nscale) = 0;
	
	% Increase paths thickness (dark red)
	
	indexx = index1(find(mod(index1,nscale) ~= 1));
	wtabs(indexx - 1) = 0.8*wtabs(indexx - 1) + 0.2;
	%wtabs(indexx - 1) = 0.55;
	wtabs(indexx - 1 + n*nscale) = 0;
	wtabs(indexx - 1 + 2*n*nscale) = 0;
	
	indexxx = index1(find(mod(index1,nscale) ~= 0));
	wtabs(indexxx + 1) = 0.8*wtabs(indexxx + 1) + 0.2;
	%wtabs(indexxx + 1) = 0.55;
	wtabs(indexxx + 1 + n*nscale) = 0;
	wtabs(indexxx + 1 + 2*n*nscale) = 0;
	
	% Markers - If length(mark)==0 this does nothing
	% Colour markers blue
	wtabs(:,mark,1) = 0;
	wtabs(:,mark,2) = 0.5;
	wtabs(:,mark,3) = 1;
	
	%---------------------------------------------------------------------------
	% Real part plot preparation
	%---------------------------------------------------------------------------
	
	% Convert to RGB and colour filtered regions black
	index3 = find(coimask == 0 | wtrealt == 0);
	[wtreal,b] = gray2ind(wtreal,128);
	wtreal = ind2rgb(wtreal,bone(128));
	wtreal(index3) = 0;
	wtreal(index3 + n*nscale) = 0;
	wtreal(index3 + 2*n*nscale) = 0;
	
	% Colour paths red
	wtreal(index1) = 0.7*wtreal(index1) + 0.3; % Transparency effect
	%wtreal(index1) = 0.85; % Constant path intensity
	wtreal(index1 + n*nscale) = 0;
	wtreal(index1 + 2*n*nscale) = 0;
	
	% Increase paths thickness (dark red)
	
	wtreal(indexx - 1) = 0.8*wtreal(indexx - 1) + 0.2;
	%wtreal(indexx - 1) = 0.55;
	wtreal(indexx - 1 + n*nscale) = 0;
	wtreal(indexx - 1 + 2*n*nscale) = 0;
	
	wtreal(indexxx + 1) = 0.8*wtreal(indexxx + 1) + 0.2;
	%wtreal(indexxx + 1) = 0.55;
	wtreal(indexxx + 1 + n*nscale) = 0;
	wtreal(indexxx + 1 + 2*n*nscale) = 0;
	
	% Markers - If length(mark)==0 this does nothing
	% Colour markers blue
	wtreal(:,mark,1) = 0;
	wtreal(:,mark,2) = 0.5;
	wtreal(:,mark,3) = 1;
	
	%---------------------------------------------------------------------------
	% Plot paths
	%---------------------------------------------------------------------------
	
	subplot(2,1,2), hold on
		
		image(timeVec(1):timeVec(end),1:size(wtreal,1),wtreal)
		xlim([timeVec(1),timeVec(end)])
		ylim([0.5,size(wtreal,1)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'YTickLabel',num2str(freqVec'));
		ylabel('Frequency (mHz)','FontSize',s2)
		xlabel('Time (s)','FontSize',s2)
		
	subplot(2,1,1), hold on
		
		image(timeVec(1):timeVec(end),1:size(wtabs,1),wtabs)
		xlim([timeVec(1),timeVec(end)])
		ylim([0.5,size(wtabs,1)+0.5])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
		set(gca,'YTickLabel',num2str(freqVec'));
		ylabel('Frequency (mHz)','FontSize',s2)
		
end


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