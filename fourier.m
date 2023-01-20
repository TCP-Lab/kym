function fourier(x,coimask,sig,freqVec,filter,ylog)

%
%--------------------------------------------------------------------------------
% Fourier Spectrum vs. Wavelet Energy Density
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% fourier(x,coimask,sig,freqVec,filter,ylog)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% x        -> matrix    -> Continuous Wavelet Modulus
% coimask  -> matrix    -> COI Mask
% sig      -> array     -> 3rd WT Output - Original Signal
% freqVec  -> array     -> Frequency Vector
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
lowest = freqVec(end)/(2^(length(freqVec)-1));

dt = 1000/(2*freqVec(end)); % Sampling time mode

nscale = size(x,1);
nvoice = (size(x,1))/(length(freqVec)-1);
n = size(x,2);

e = exp(1); % Essential just for MATLAB

x = x .* coimask;

NRG = (10^3)*sum(x.^2,2)*dt; % Normalized by sampling frequency 1/dt

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
% It must be N > 2^length(freqVec)
N = 2^14;
FT = fft(sig,N)*dt; % Normalized by sampling frequency 1/dt

% Fourier Spectrum
FT = abs(FT).^2;
FT = FT(1:N/2+1);
FT(2:N/2) = 2*FT(2:N/2); % The factor 2 compensates for missing negative values

% Log Scale for y-axis: z = log(y+1)
if (ylog == 1)
	NRG = log(NRG+1);
	FT = log(FT+1);
end

% Assemble Fourier Spectrum x-axis
xax = ([0:N/2]./(N/2)).*freqVec(end);
if isempty(xax == lowest)
	fprintf('\n\nWARNING: No matching between lowest wavelet frequency and Fourier lowest frequency\n');
end
FT = FT(find(xax >= lowest));
xax = xax(find(xax >= lowest));

% Plot Fourier Spectrum vs. Wavelet Energy Density
plot([log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(NRG)-1):log2(xax(end))],NRG,'-r','LineWidth',1), hold on
lgn{1} = ['10^3 Wavelet Energy Density']; % Cell Array
plot(log2(xax),FT,'-b','LineWidth',1)
lgn{2} = ['Fourier Power Spectrum']; % Cell Array

% Alternative Syntax
%[graph,h1,h2] = plotyy([log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(NRG)-1):log2(xax(end))],NRG,log2(xax),FT);

xlim([log2(xax(1)),log2(xax(end))])
set(gca,'FontSize',s1,'XTick',unique([log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(freqVec)-1):log2(xax(end))]));
set(gca,'XTickLabel',num2str(freqVec'));
% Log Scale for y-axis: z = log(y+1) -> y = e^z-1
if (ylog == 1)
	set(gca,'YTickLabel',num2str((floor((e.^(get(gca,'YTick'))-1)*100)/100)')); % In order to get only 1 decimal digit
end

xlabel('Frequency (mHz)','FontSize',s2)
ylabel('Amplitudes','FontSize',s2)

mylegend = legend(lgn,'Location','NorthEast');
set(mylegend,'FontSize',s1);


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