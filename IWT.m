function IWT(wt,par,sig,fwt,output)

% 
%--------------------------------------------------------------------------------
% Inverse Wavelet Transform - Reconstruction & Synthesis
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% IWT(wt,par,sig,fwt,output)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% wt       -> matrix    -> 1st WT Output - Continuous Wavelet Transform
% par      -> structure -> 2nd WT Output - Parameters
% sig      -> array     -> 3rd WT Output - Calcium Signal
% fwt      -> matrix    -> FILT Output - Filtered Wavelet Transform
% output   -> boolean   -> Print .eps Output Graphs
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% -none-   -> plot      -> 1 Plot Resulting from Analysis
%
% MEANINGS: Read the plot following this order (MATLAB subplot numbering):
% ========= Subplot 3 shows the starting signal in blue.
%           Subplot 5 shows its related scaleogram.
%           Subplot 3 also shows the signal reconstructed from this
%              scaleogram in red: it can differ from the blue one
%              because of the poor frequency resolution or the discarded
%              regions (low frequencies cutting, COI cutting, PAD cutting).
%           Subplot 6 shows the filtered scaleogram.
%           Subplot 4 shows the signal reconstructed from this in black.
%           Subplot [1 2] shows the three signals superimposed.
%

%---------------------------------------------------------------------------
% Graphic Parameters
%---------------------------------------------------------------------------

s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

%---------------------------------------------------------------------------
% Default Output
%---------------------------------------------------------------------------

ratio = [];

%---------------------------------------------------------------------------
% Variable Assignment
%---------------------------------------------------------------------------

coimask = par.k;
timeVec = par.x;
freqVec = par.y;
periodVec = par.z;
lowest = freqVec(end)/(2^(length(freqVec)-1));
mark = par.m;

% CWT Modulus
%wtabs = abs(wt);
%x = wtabs .* coimask;

% CWT Real part (absolute value)
%wtreal = abs(real(wt));
%y = wtreal .* coimask;

% CWT real part
wtr = real(wt);
z = wtr .* coimask;

% Filtered CWT real part
fwtr = real(fwt);
fz = fwtr .* coimask;

nscale = size(wt,1);
nvoice = size(wt,1)/(length(freqVec)-1);
n = size(wt,2);

%---------------------------------------------------------------------------
% Input Control
%---------------------------------------------------------------------------

% Print plots - Default value = false = Print nothing
if (nargin < 5)
	output = false;
end

% Filtered wavelet - Default value = wt = NO FILTER
if (nargin < 4)
	fwt = wt;
end

%---------------------------------------------------------------------------
% IWT - Inverse Wavelet Transform - Signal Reconstruction
%---------------------------------------------------------------------------

synth1 = sum(z);
synth2 = sum(fz);

%---------------------------------------------------------------------------
% Rescale to [0,1] and convert to RGB
%---------------------------------------------------------------------------

wtr = (wtr-min(min(wtr)))/max(max(wtr-min(min(wtr))));
[wtrind,b] = gray2ind(wtr,128);
wtr = ind2rgb(wtrind,bone(128));

fwtr = (fwtr-min(min(fwtr)))/max(max(fwtr-min(min(fwtr))));
[fwtrind,b] = gray2ind(fwtr,128);
fwtr = ind2rgb(fwtrind,bone(128));

%---------------------------------------------------------------------------
% Apply Gray Mask (Filter + COI) to wtr and fwtr
%---------------------------------------------------------------------------

% Gray scale scaleogram_1
gswtr = ind2rgb(wtrind,gray(128));

% Shade scaleogram_1 masked regions gray
index1 = find(coimask == 0);

% Gray scale - If length(index)==0 this does nothing
wtr(index1) = gswtr(index1);
wtr(index1 + n*nscale) = gswtr(index1 + n*nscale);
wtr(index1 + 2*n*nscale) = gswtr(index1 + 2*n*nscale);

% Delete masked regions from scleogram_2
index2 = find(fwt == 0);

% Light Gray - If length(index)==0 this does nothing
fwtr(index1) = 0.8;
fwtr(index1 + n*nscale) = 0.8;
fwtr(index1 + 2*n*nscale) = 0.8;

fwtr(index2) = 0.8;
fwtr(index2 + n*nscale) = 0.8;
fwtr(index2 + 2*n*nscale) = 0.8;

%---------------------------------------------------------------------------
% Rescale Signals for Superimposition
%---------------------------------------------------------------------------

amp0 = max(sig)-min(sig);
amp1 = max(synth1)-min(synth1);
amp2 = max(synth2)-min(synth2);
synth1 = synth1*(amp0/amp1);
synth2 = synth2*(amp0/amp2);
synth1 = synth1 - mean(synth1) + mean(sig);
synth2 = synth2 - mean(synth2) + mean(sig);

%---------------------------------------------------------------------------
% Plot
%---------------------------------------------------------------------------

figure

subplot(3,2,[1 2]), hold on
	
	plot(timeVec,sig,'-b','LineWidth',1)
	plot(timeVec,synth1,'-r','LineWidth',1)
	plot(timeVec,synth2,'-k','LineWidth',1)
	
	xlim([timeVec(1),timeVec(end)])
	q = ylim;
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Ratio','FontSize',s2)
	title('Inverse Wavelet Transform - Reconstruction and Synthesis','FontSize',s3)
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ', num2str(par.r)],'FontSize',s2,'Color','k')
	
	for k = 1:length(mark) % If length(mark)==0 this does nothing
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
subplot(3,2,3), hold on
	
	plot(timeVec,sig,'-b','LineWidth',1)
	plot(timeVec,synth1,'-r','LineWidth',1)
	
	xlim([timeVec(1),timeVec(end)])
	ylim(q) % Fix y-range to enable comparison among different plots
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Ratio','FontSize',s2)
	
	for k = 1:length(mark) % If length(mark)==0 this does nothing
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
subplot(3,2,5), hold on
	
	image(timeVec(1):timeVec(end),1:size(wtr,1),wtr)
	xlim([timeVec(1),timeVec(end)])
	ylim([0.5,size(wtr,1)+0.5])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'YTickLabel',num2str(freqVec'));
	xlabel('Time (s)','FontSize',s2)
	ylabel('Frequency (mHz)','FontSize',s2)
	
subplot(3,2,6), hold on
	
	image(timeVec(1):timeVec(end),1:size(fwtr,1),fwtr)
	xlim([timeVec(1),timeVec(end)])
	ylim([0.5,size(fwtr,1)+0.5])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'YTickLabel',num2str(freqVec'));
	xlabel('Time (s)','FontSize',s2)
	
subplot(3,2,4), hold on
	
	plot(timeVec,synth2,'-k','LineWidth',1)
	xlim([timeVec(1),timeVec(end)])
	ylim(q) % Fix y-range to enable comparison among different plots
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	
	for k = 1:length(mark) % If length(mark)==0 this does nothing
		plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
% Print output graph
if (output)
	print -depsc reconout.eps
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