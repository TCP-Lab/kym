function maxima = paths(x,y,coimask,thr,unit,US,timeVec,freqVec,mark,minDist,noPlateau,grafunc)

%
%--------------------------------------------------------------------------------
% NU CWT Frequency Paths
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% maxima = paths(x,y,coimask,thr,unit,US,timeVec,freqVec,mark,minDist,noPlateau,grafunc)
%
% INPUT        TYPE       MEANING
% -----        ----       -------
% x         -> matrix  -> Continuous Wavelet Modulus
% y         -> matrix  -> Continuous Wavelet Real Part
% coimask   -> matrix  -> COI Mask
% thr       -> scalar  -> Threshold Percentage
% unit      -> string  -> Time Unit: 's' or 'ms'
% US        -> scalar  -> US-Upsampling factor
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
% -none-    -> plot     -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
[nscale,n] = size(x);
nvoice = nscale/(length(freqVec)-1);

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

% Threshold, COI and "ultra-Nyquist zone"
threshold = thr/100;
wtabst = (wtabs >= threshold);
wtrealt = (wtreal >= threshold);
%wtrealt = wtabst; % To use the same threshold contours used for modulus
maxima = maxima .* wtabst .* coimask .* [ones(nscale-log2(US)*nvoice,n);zeros(log2(US)*nvoice,n)];

% Enable graphical functions
if (grafunc)
	
	% Find maxima
	%[rx,cx] = find(maxima);
	index1 = find(maxima);
	
	%---------------------------------------------------------------------------
	% Modulus plot preparation
	%---------------------------------------------------------------------------
	
	% Convert to RGB and black COI region
	index2 = find(coimask == 0);
	%index2 = find(coimask == 0 | wtabst == 0); % Convert to RGB and black also thresholded regions
	[wtabs,b] = gray2ind(wtabs,128);
	wtabs = ind2rgb(wtabs,bone(128));	
	wtabs(index2) = 0;
	wtabs(index2 + n*nscale) = 0;
	wtabs(index2 + 2*n*nscale) = 0;
	
	% Blue Nyquist limit
	wtabs(end-log2(US)*nvoice,:,1) = 0;
	wtabs(end-log2(US)*nvoice,:,2) = 0.5;
	wtabs(end-log2(US)*nvoice,:,3) = 1;
	
	% Blue markers
	wtabs(:,mark,1) = 0;
	wtabs(:,mark,2) = 0.5;
	wtabs(:,mark,3) = 1;
	
	% Red paths
	wtabs(index1) = 0.5*wtabs(index1) + 0.5; % Transparency
	%wtabs(index1) = 0.85; % Constant intensity
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
	
	%---------------------------------------------------------------------------
	% Real part plot preparation
	%---------------------------------------------------------------------------
	
	% Convert to RGB and black COI region
	index3 = find(coimask == 0);
	%index3 = find(coimask == 0 | wtrealt == 0); % Convert to RGB and black also thresholded regions
	[wtreal,b] = gray2ind(wtreal,128);
	wtreal = ind2rgb(wtreal,bone(128));
	wtreal(index3) = 0;
	wtreal(index3 + n*nscale) = 0;
	wtreal(index3 + 2*n*nscale) = 0;
	
	% Blue Nyquist limit
	wtreal(end-log2(US)*nvoice,:,1) = 0;
	wtreal(end-log2(US)*nvoice,:,2) = 0.5;
	wtreal(end-log2(US)*nvoice,:,3) = 1;
	
	% Blue markers
	wtreal(:,mark,1) = 0;
	wtreal(:,mark,2) = 0.5;
	wtreal(:,mark,3) = 1;
	
	% Red paths
	wtreal(index1) = 0.5*wtreal(index1) + 0.5; % Transparency
	%wtreal(index1) = 0.85; % Constant intensity
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
	
	%---------------------------------------------------------------------------
	% Plot paths
	%---------------------------------------------------------------------------
	
	subplot(2,1,2), hold on
		
		image(timeVec,1:nscale,wtreal)
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
	
	subplot(2,1,1), hold on
		
		image(timeVec,1:nscale,wtabs)
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