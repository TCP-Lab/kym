function energy(x,coimask,timeVec,freqVec,nvoice,mark)

%
%--------------------------------------------------------------------------------
% Wavelet Power Spectrum
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% energy(x,coimask,timeVec,freqVec,nvoice,mark)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% x        -> matrix   -> Continuous Wavelet Modulus
% coimask  -> matrix   -> COI Mask
% timeVec  -> array    -> Time Vector
% freqVec  -> array    -> Frequency Vector
% nvoice   -> scalar   -> Lines of Pixel per Octave
% mark     -> array    -> Time Marker Set (Step Number)
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% -none-   -> plot     -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

dt = mode(diff(timeVec)); % Sampling time mode

nscale = size(x,1);
n = size(x,2);

wtabs = x;
x = x .* coimask;

% Rescale to [0,1] and convert to RGB
wtabs = (wtabs-min(min(wtabs)))/max(max(wtabs-min(min(wtabs))));
[wtabsind,b] = gray2ind(wtabs,128);
wtabs = ind2rgb(wtabsind,jet(128));
gswtabs = ind2rgb(wtabsind,gray(128));

% Shade masked regions gray - If length(index)==0 this does nothing
index = find(coimask == 0);
wtabs(index) = gswtabs(index);
wtabs(index + n*nscale) = gswtabs(index + n*nscale);
wtabs(index + 2*n*nscale) = gswtabs(index + 2*n*nscale);

PWR = (10^5)*sum(x.^2,1)/(nvoice); % Normalized by number of voices per octave
NRG = (10^2)*sum(x.^2,2)*dt; % Normalized by sampling frequency 1/dt

if (length(mark) > 0)
	% Power mean values
	MEAN(1) = (1/mark(1))*sum(PWR(1:mark(1)));
	if (length(mark) > 1)
		for h = 2:length(mark)
			MEAN(h) = (1/(mark(h)-mark(h-1)+1))*sum(PWR(mark(h-1):mark(h)));
		end
	end
	MEAN(length(mark)+1) = (1/(n-mark(end)+1))*sum(PWR(mark(end):end));
	MEAN = floor(MEAN*100)/100; % In order to get only 2 decimal digits
end

% Plotting
subplot(2,2,2), hold on
	
	plot(timeVec,PWR,'-b','LineWidth',1)
	
	xlim([timeVec(1),timeVec(end)])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	set(gca,'XTickLabel','');
	ylabel('10^5 Power P','FontSize',s2)
	title('Wavelet Power Spectrum','FontSize',s3)
	
	if (length(mark) > 0)
		text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(1))],'FontSize',s2,'Color','k')
		if (length(mark) > 1)
			for h = 2:length(mark)
				text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(h))],'FontSize',s2,'Color','k')
			end
		end
		text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(end))],'FontSize',s2,'Color','k')
		
		r = MEAN(2)/MEAN(1);
		r = floor(r*100)/100; % In order to get only 2 decimal digits
		text(min(xlim)-(max(xlim)-min(xlim))*(35/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(r)],'FontSize',s2,'Color','r')
	end
	
	for k = 1:length(mark) % If length(mark)==0 this does nothing
			plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
subplot(2,2,3), hold on
	
	plot(NRG,[1:nscale],'-b','LineWidth',1)
	ylim([1,nscale])
	set(gca,'XDir','reverse');
	set(gca,'FontSize',s1,'YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'YTickLabel',num2str(freqVec'));
	ylabel('Frequency (mHz)','FontSize',s2)
	xlabel('10^2 Energy Density E','FontSize',s2)

subplot(2,2,4), hold on
	
	imagesc(timeVec(1):timeVec(end),1:size(wtabs,1),wtabs)
	xlim([timeVec(1),timeVec(end)])
	ylim([0.5,size(wtabs,1)+0.5])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'YTickLabel','');
	xlabel('Time (s)','FontSize',s2)
	
% Resize Plot
subplot(2,2,2)
set(gca,'Position',[0.35,0.7,0.62,0.2])

subplot(2,2,3)
set(gca,'Position',[0.12,0.11,0.2,0.55])

subplot(2,2,4)
set(gca,'Position',[0.35,0.11,0.62,0.55])


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