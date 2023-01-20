function fourier(x,sig,freqVec,filter)

%
%---------------------------------------------------------------------------
% Fourier Spectrum vs. Wavelet Energy Density
%---------------------------------------------------------------------------
%
%
% Function Definition
%
% fourier(x,sig,freqVec,filter)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% x        -> matrix    -> Morlet Wavelet Modulus
% sig      -> array     -> 3rd WT Output - Original Signal
% freqVec  -> array     -> Frequency Vector
% filter   -> scalar    -> Filter Type for the Original Signal
%

% Variables Assignment
lowest = freqVec(end)/(2**(length(freqVec)-1));

nscale = size(x)(1);
nvoice = (size(x)(1))/(length(freqVec)-1);
n = size(x)(2);

NRG = (10**3)*sum(x.**2,2);

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
		if (j<length(sig)*(25/100) || j>length(sig)-length(sig)*(25/100))
			sig(j) = sig(j).*hann(j);
		endif
	endfor
	
endswitch

% Fast Fourier Transform
% It must be N > 2**length(freqVec)
%N = 2**12;
N = 2**14;
FT = fft(sig,N);

% Fourier Spectrum - The factor 2 compensates for missing negative values
FT = abs(FT).**2;
FT = FT(1:N/2+1);
FT(2:N/2) = 2*FT(2:N/2);

% Assemble Fourier Spectrum x-axis
xax = [0:N/2]./(N/2).*freqVec(end);
if (sum(xax == lowest) != 1)
	printf("\nWARNING!\nNo matching between lowest wavelet frequency and Fourier lowest frequency\n\n");
endif
FT = FT(find(xax >= lowest));
xax = xax(find(xax >= lowest));

% Plot Fourier Spectrum 
[graph,h1,h2] = plotyy([log2(xax(1)):(log2(xax(end))-log2(xax(1)))./(length(NRG)-1):log2(xax(end))],NRG,log2(xax),FT);

set(graph(1),'XTick',[log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(freqVec)-1):log2(xax(end))]);
set(graph(2),'XTick',[log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(freqVec)-1):log2(xax(end))]);

set(graph(1),'XTickLabel',freqVec);
set(graph(2),'XTickLabel',freqVec);

ylabel(graph(1),'Wavelet Energy Density','FontSize',18)
ylabel(graph(2),'Fourier Spectrum','FontSize',18)
xlabel('Frequency (mHz)')

set(graph(1),'YColor',[1,0,0]);
set(graph(2),'YColor',[0,0,1]);

set(h1,'color',[1,0,0])
set(h2,'color',[0,0,1])

% Alternative Syntax
axes(graph(1))
%set(gca,'XTick',[log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(freqVec)-1):log2(xax(end))]);
%set(gca,'XTickLabel',freqVec);
%ylabel('Wavelet Energy Density','FontSize',18)
%set(gca,'YColor',[1,0,0]);
legend ('Wavelet Energy Density','Location','NorthEast')

axes(graph(2))
%set(gca,'XTick',[log2(xax(1)):(log2(xax(end))-log2(xax(1)))/(length(freqVec)-1):log2(xax(end))]);
%set(gca,'XTickLabel',freqVec);
%ylabel('Fourier Spectrum','FontSize',18)
%set(gca,'YColor',[0,0,1]);
legend('Fourier Spectrum','Location','East') 

% Problemi col ridimensionamento dei grafici plotyy... alcune volte funziona...
% ma anche in quei casi la differenza di 0.005 fra le larghezze dei 2 plot
% è necessaria per avere gli assi allineati!
%set(graph(1),'position',[0.130,0.110,0.77,0.68]);
%set(graph(2),'position',[0.130,0.110,0.775,0.68]);


%---------------------------------------------------------------------%
%                                                                     %
% A.A. 2009 / 2010                                                    %
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