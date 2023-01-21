function energy(x,coimask,unit,US,timeVec,freqVec,mark,admis)

%
%--------------------------------------------------------------------------------
% Unitary Wavelet Amplitude Spectral Density (ASD)
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% energy(x,coimask,unit,US,timeVec,freqVec,mark,admis)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% x        -> matrix   -> Continuous Wavelet Modulus
% coimask  -> matrix   -> COI Mask
% unit     -> string   -> Time Unit: 's' or 'ms'
% US       -> scalar   -> US-Upsampling factor
% timeVec  -> array    -> Time Vector
% freqVec  -> array    -> Frequency Vector
% mark     -> array    -> Time Marker Set (Step Number)
% admis    -> scalar   -> Admissibility Condition Factor
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% -none-   -> plot     -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
[nscale,n] = size(x);
nvoice = nscale/(length(freqVec)-1);
lowest = freqVec(end)/(2^(length(freqVec)-1));

% Sampling time mode in s
dt = 1/(2*freqVec(end));
if (strcmp(unit,'s'))
	dt = dt*1000;
end

% Complete Frequency Vector in Hz
freqVecComp = freqVec(end)*(2.^(-[0:1:nscale-1]/nvoice));
freqVecComp = flipud(freqVecComp');
if (strcmp(unit,'s'))
	freqVecComp = freqVecComp/1000;
end

% Complete Scale Vector in s
scaleVecComp = 2*2.^([0:1:nscale-1]/nvoice)*dt;
scaleVecComp = flipud(scaleVecComp');
scaleVecComp = scaleVecComp*ones(1,n);

% Unitary Transform - Energetic Normalization
x = x.*sqrt(scaleVecComp);

% Apply COI mask
wtabs = x;
x = x.*coimask;

% Power Time Course P(t)
% dv = log(2)*v*dlog2(v) with dlog2(v)=1/nvoice
PWR = admis*sum((x.^2).*(freqVecComp*ones(1,n)),1)*(log(2)/nvoice); % Robustness with respect to nvoice
% The factor 2 compensates for missing negative values in the scalogram frequencies
PWR = 2*PWR; % THIS IS CORRECT FOR PURELY REAL SIGNALS ONLY !!!

% Energy Spectral Density (ESD)
NRG = sum(x.^2,2)*dt; % Robustness with respect to sampling time
% The factor 2 compensates for missing negative values in the scalogram frequencies
NRG = 2*NRG; % THIS IS CORRECT FOR PURELY REAL SIGNALS ONLY !!!

% Conservation of energy
ener1 = sum(PWR,2)*dt;
ener2 = admis*sum(NRG.*freqVecComp,1)*(log(2)/nvoice);
fprintf('\n\nConservation of Energy\n');
fprintf('\n\n SUM P(t) dt = %.4f J\n',ener1);
fprintf('\n\n SUM ESD(v) dv = %.4f J\n',ener2);
fprintf('\n\n');

% Power mean values
if (length(mark) > 0)
	MEAN(1) = (1/mark(1))*sum(PWR(1:mark(1)));
	if (length(mark) > 1)
		for h = 2:length(mark)
			MEAN(h) = (1/(mark(h)-mark(h-1)+1))*sum(PWR(mark(h-1):mark(h)));
		end
	end
	MEAN(length(mark)+1) = (1/(n-mark(end)+1))*sum(PWR(mark(end):end));
end

% Rescale to [0,1] and convert to RGB
wtabs = (wtabs-min(min(wtabs)))/max(max(wtabs-min(min(wtabs))));
[wtabsind,b] = gray2ind(wtabs,128);
wtabs = ind2rgb(wtabsind,jet(128));
gswtabs = ind2rgb(wtabsind,gray(128));

% Shade the "ultra-Nyquist zone": the log2(US) highest octaves
wtabs(end-log2(US)*nvoice+1:end,:,:) = gswtabs(end-log2(US)*nvoice+1:end,:,:);

% Shade masked regions gray
index = find(coimask == 0);

% Gray scale - If length(index)==0 this does nothing
wtabs(index) = gswtabs(index);
wtabs(index + n*nscale) = gswtabs(index + n*nscale);
wtabs(index + 2*n*nscale) = gswtabs(index + 2*n*nscale);

% Plotting
subplot(2,2,2), hold on
	
	plot(timeVec,PWR,'-b','LineWidth',1)
	
	xlim([timeVec(1),timeVec(end)])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	set(gca,'XTickLabel','');
	ylabel('Power (W)','FontSize',s2)
	title('Unitary CWT ASD (J/sHz)^{1/2}','FontSize',s3)
	
	if (length(mark) > 0)
		text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(1),'%10.2e')],'FontSize',s2,'Color','k')
		if (length(mark) > 1)
			for h = 2:length(mark)
				text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(h),'%10.2e')],'FontSize',s2,'Color','k')
			end
		end
		text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(end),'%10.2e')],'FontSize',s2,'Color','k')
		
		r = MEAN(2)/MEAN(1);
		text(min(xlim)-(max(xlim)-min(xlim))*(35/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(r,2)],'FontSize',s2,'Color','r')
	end
	
	for k = 1:length(mark) % If length(mark)==0 this does nothing
			plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	end
	
subplot(2,2,3), hold on
	
	plot(NRG,[1:nscale],'-b','LineWidth',1)
	% Nyquist limit: the log2(US) highest octaves
	plot([min(xlim),max(xlim)],[nscale-log2(US)*nvoice,nscale-log2(US)*nvoice],'-k','LineWidth',1)
	ylim([1,nscale])
	set(gca,'XDir','reverse');
	set(gca,'FontSize',s1,'YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
	set(gca,'YTickLabel',num2str(freqVec',4));
	if (strcmp(unit,'ms'))
		ylabel('Frequency (Hz)','FontSize',s2)
	else
		ylabel('Frequency (mHz)','FontSize',s2)
	end
	xlabel('ESD (J/Hz)','FontSize',s2)

subplot(2,2,4), hold on
	
	image(timeVec,1:nscale,wtabs)
	xlim([timeVec(1),timeVec(end)])
	ylim([0.5,nscale+0.5])
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