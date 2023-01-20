function tang = hurst(NRG,sig,unit,US,freqVec,roi)

%
%--------------------------------------------------------------------------------
% Hurst Exponent Computation (H) and Diffusion Coefficient (D)
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% tang = hurst(NRG,sig,unit,US,freqVec,roi)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% NRG      -> array     -> Energy Spectral Density (ESD)
% sig      -> array     -> 3rd WT Output - Original Signal
% unit     -> string    -> Time Unit: 's' or 'ms'
% US       -> scalar    -> US-Upsampling factor
% freqVec  -> array     -> Frequency Vector
% roi      -> scalar    -> Last ROI from WT
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% tang     -> array
% -none-   -> plot      -> Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
nscale = length(NRG);
nvoice = nscale/(length(freqVec)-1);
lowest = freqVec(end)/(2^(length(freqVec)-1));

% Sampling time mode
dt = 1000/(2*freqVec(end));
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

%---------------------------------------------------------------------------
% Hurst Exponent Computation (H)
%---------------------------------------------------------------------------

k = 3; % Number of neglected octaves (starting from the Nyquist frequency)
h = 2; % Number of octaves used for interpolation
interval = [nscale-round((log2(US)+k+h)*nvoice):nscale-round((log2(US)+k)*nvoice)];

% Linear fit
[coeff,s] = polyfit(logfreqVecComp(interval),NRG(interval),1);
stdev = sqrt(diag((inv(s.R)*inv(s.R'))*s.normr^2/s.df));
% From MathWorks documentation:
% [p,S] = polyfit(x,y,n)
% S (Error estimation structure) contains the following fields:
% R		-> Triangular factor from a QR decomposition of the Vandermonde matrix of x
% df	-> Degrees of freedom
% normr -> Norm of the residuals
% "[...] an estimate of the covariance matrix of p is (Rinv*Rinv')*normr^2/df, where Rinv is the inverse of R."

coeffgraph = polyfit(interval',NRG(interval),1);
tang = polyval(coeffgraph,[1:nscale]);

fprintf('\n\nSpectral Index 1/f^(beta)\n');
fprintf('\n beta = %.4f +/- %.4f\n',-coeff(1),stdev(1));

return

% alternative computation for uncertainty according to the standard formula
sqrt(sum((polyval(coeff,logfreqVecComp(interval))-NRG(interval)).^2)/(length(interval)-2)/sum((logfreqVecComp(interval)-mean(logfreqVecComp(interval))).^2))




if (-coeff(1)<1)
	
	fprintf('\n\nFractional Gaussian noise\n');
	fprintf('\n H = %.4f +/- %.4f\n',(-coeff(1)+1)/2,0.5*stdev(1)); % According to Propagation of Uncertainty
	
elseif (-coeff(1)>1)
	
	fprintf('\n\nFractional Brownian motion\n');
	fprintf('\n H = %.4f +/- %.4f\n',(-coeff(1)-1)/2,0.5*stdev(1)); % According to Propagation of Uncertainty
	
	%---------------------------------------------------------------------------
	% Diffusion Coefficient Computation (D)
	%---------------------------------------------------------------------------
	
	% Time Averaged Mean Squared Displacement (TAMSD)
	for m = 0:1:sqrt(length(sig))
		msd(m+1) = mean((sig(m+1:end)-sig(1:end-m)).^2);
	end
	delay = dt*[0:1:length(msd)-1];
	
	%delay(2:end)'
	%msd(2:end)'
	%return
	
	% Non-linear fit
	F = @(p,t)(2*p(1)*t.^(-coeff(1)-1));
	p0 = [10];
	[p,resnorm,resid,exitflag,output,lambda,J]= lsqcurvefit(F,p0,delay,msd);
	ci = nlparci(p,resid,'jacobian',J,'alpha',0.3173) % Returns 68.27% (=1 Standard Deviation under Normality assumption) confidence intervals on parameters from lsqcurvefit
	pstdev = (ci(2)-ci(1))/2;
	if (exitflag==1)
		fprintf('<lsqcurvefit> converged to a solution\n');
	else
		fprintf('WARNING: some problem may be occured during interpolation!\n');
		fprintf('Check exitflag value in lsqcurvefit output.\n');
	end
	
	fprintf('\n\nDiffusion Coefficient\n');
	fprintf('\n D = %.4f +/- %.4f\n',p,pstdev);
	
	log2gammaD = log(2*gamma(1-(-coeff(1)-1))*p);
	log2gammaDstdev = sqrt((1/p^2)*pstdev^2+(psi(1-(-coeff(1)-1)))^2*stdev(1)^2); % According to Propagation of Uncertainty
	fprintf('\n log(2GammaD) = %.4f +/- %.4f\n',log2gammaD,log2gammaDstdev);
	
	figure
	hold on
	plot(delay,msd)
	plot(delay,F(p,delay))
	
	xlim([delay(1),delay(end)])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(delay(end))]));
	if (strcmp(unit,'ms'))
		xlabel('Time (ms)','FontSize',s2)
	else
		xlabel('Time (s)','FontSize',s2)
	end
	ylabel('Mean Squared Displacement','FontSize',s2)
	title(['Time Averaged Mean Squared Displacement - ROI ',num2str(roi)],'FontSize',s3)
	hold off

elseif (-coeff(1)==1)
	
	fprintf('\n\nPure Pink (1/f) noise!\n');
	
end

fprintf('\n\n');


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