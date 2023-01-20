function kymfourier(x,coimask,sig,unit,US,freqVec,roi,admis,filter,ylog)

%
%--------------------------------------------------------------------------------
% Fourier Spectrum vs. Wavelet Energy Spectral Density and Hurst Exponent
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% kymfourier(x,coimask,sig,unit,US,freqVec,admis,filter,ylog)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% x        -> matrix    -> Continuous Wavelet Modulus
% coimask  -> matrix    -> COI Mask
% sig      -> array     -> 3rd WT Output - Original Signal
% unit     -> string    -> Time Unit: 's' or 'ms'
% US       -> scalar    -> US-Upsampling factor
% freqVec  -> array     -> Frequency Vector
% roi      -> scalar    -> Last ROI from WT
% admis    -> scalar    -> Admissibility Condition Factor
% filter   -> scalar    -> Filter Type for the Original Signal
% ylog     -> scalar    -> 0 = Linear y-Scale / 1 = Log y-Scale
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% -none-   -> plot      -> Plot Resulting from Analysis
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
logfreqVecComp = log2(freqVecComp);

% Complete Scale Vector in s
scaleVecComp = 2*2.^([0:1:nscale-1]/nvoice)*dt;
scaleVecComp = flipud(scaleVecComp');
scaleVecComp = scaleVecComp*ones(1,n);

%---------------------------------------------------------------------------
% Spectra Computation
%---------------------------------------------------------------------------

% Apply COI mask
x = x.*coimask;

% Alternative, equivalent to ener3
% ener4 = admis*sum(2*sum(x.^2,2)*dt,1)*(log(2)/nvoice); % Robustness with respect to sampling time and nvoice

% Unitary Transform - Energetic Normalization
x = x.*sqrt(scaleVecComp);

% Energy Spectral Density (ESD)
NRG = sum(x.^2,2)*dt; % Robustness with respect to sampling time
% The factor 2 compensates for missing negative values in the scalogram frequencies
NRG = 2*NRG; % THIS IS CORRECT FOR PURELY REAL SIGNALS ONLY !!!

% Remove DC offset and filter
sig = sig - mean(sig);

switch (filter)
	
	case 0
	% No filter
	
	case 1
	% Filtering the signal - Hann window
	hann = 0.5 - 0.5*cos(2*pi*([1:length(sig)]')/length(sig));
	sig = sig.*hann;
	
	case 2
	% Filtering the signal - Tukey window
	hann = 0.5 - 0.5*cos(2*pi*([1:length(sig)]')/length(sig));
	for j = 1:length(sig)
		if (j < length(sig)*(25/100) || j > length(sig)-length(sig)*(25/100))
			sig(j) = sig(j).*hann(j);
		end
	end
	
end

% Fast Fourier Transform
% It must be N >= 2^length(freqVec)
N = 2^(length(freqVec)+5);
FT = fft(sig,N)*dt; % Robustness with respect to sampling time

% For real input data, Hermitian symmetry holds: F(-w)=F*(w) -> Single-Sided Amplitude Spectrum
FT = FT(1:N/2+1); % THIS IS CORRECT FOR PURELY REAL SIGNALS ONLY !!!
% Fourier Spectrum
FT = abs(FT).^2;
% The factor 2 compensates for missing negative values
FT(2:N/2) = 2*FT(2:N/2);

%---------------------------------------------------------------------------
% Parseval's theorem (Conservation of energy)
%---------------------------------------------------------------------------

ener1 = sum(abs(sig).^2,1)*dt;
ener2 = sum(FT,1)*(1/(2*dt*length(FT)));
ener3 = admis*sum(NRG.*freqVecComp,1)*(log(2)/nvoice); % dv = log(2)*v*dlog2(v) with dlog2(v)=1/nvoice
fprintf('\n\nParseval''s Theorem - Conservation of Energy\n');
fprintf('\n\n SUM |f(t)|^2 dt = %.4f J\n',ener1);
fprintf('\n\n SUM |F(v)|^2 dv = %.4f J\n',ener2);
fprintf('\n\n SUM SUM |W(t,v)|^2 dtdv = %.4f J\n',ener3);
%fprintf('\n\n SUM SUM |W(t,v)|^2 dtdv = %.4f J\n',ener4); % Alternative, equivalent to ener3
fprintf('\n\n');

%---------------------------------------------------------------------------
% Log Scale for y-axis and Spectral Index (beta)
%---------------------------------------------------------------------------

if (ylog == 1)
	
	NRG = log2(NRG);
	FT = log2(FT);
	
	tang = hurst(NRG,sig,unit,US,freqVec,roi);
	
end

%---------------------------------------------------------------------------
% Plot Spectra
%---------------------------------------------------------------------------

% Assemble Fourier Spectrum x-axis
xax = ([0:N/2]/(N/2)).*freqVec(end);
if isempty(find((xax == lowest)))
	fprintf('WARNING: No matching between wavelet lowest frequency and Fourier lowest frequency\n');
	fprintf('\n\n');
end
FT = FT(find(xax >= lowest));
xax = xax(find(xax >= lowest));

figure

% Plot Single-Sided Fourier Spectrum vs. Wavelet Energy Spectral Density
plot(log2(xax),FT,'-b','LineWidth',1), hold on
lgn{1} = ['Fourier Power Spectrum']; % Cell Array
plot([log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(NRG)-1):log2(xax(end))],NRG,'-r','LineWidth',2)
lgn{2} = ['Wavelet Energy Spectral Density']; % Cell Array

if (ylog == 1)
	% Plot 1/f trend
	plot([log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(NRG)-1):log2(xax(end))],tang,'-g','LineWidth',1)
end

xlim([log2(xax(1)),log2(xax(end))])
set(gca,'FontSize',s1,'XTick',unique([log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(freqVec)-1):log2(xax(end))]));
set(gca,'XTickLabel',num2str(freqVec',4));

% Nyquist limit: the log2(US) highest octaves
plot([log2(freqVec(end)/US),log2(freqVec(end)/US)],[min(ylim),max(ylim)],'-k','LineWidth',1)

% Log Scale for y-axis: z = log2(y) -> y = 2^z
if (ylog == 1)
	set(gca,'YTickLabel',num2str((2.^(get(gca,'YTick')))',4));
end

if (strcmp(unit,'ms'))
	xlabel('Frequency (Hz)','FontSize',s2)
else
	xlabel('Frequency (mHz)','FontSize',s2)
end
ylabel('Amplitude (J/Hz)','FontSize',s2)

mylegend = legend(lgn,'Location','NorthEast');
set(mylegend,'FontSize',s1);

title(['Fourier Spectrum vs. Wavelet ESD - ROI ',num2str(roi)],'FontSize',s3)


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