function [data,timeVec,dt,outroi] = loader(filename,unit,roi,DS)

%
%--------------------------------------------------------------------------------
% Data Loader
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [data,timeVec,dt,outroi] = loader(filename,unit,roi,DS)
%
% INPUT       TYPE      MEANING
% -----       ----      -------
% filename -> string -> File Name or Test type
% unit     -> string -> Time Unit: 's' or 'ms'
% roi      -> array  -> ROI Set ([x:y]==from x to y)
% DS       -> scalar -> Downsampling Rate
%
% OUTPUT      TYPE      MEANING
% ------      ----      -------
% data     -> matrix -> Matrix of Raw Data to be Analyzed
% timeVec  -> array  -> Time Vector
% dt       -> scalar -> Sampling Time
% outroi   -> array  -> Selected ROIs
%

%---------------------------------------------------------------------------
% Default Output
%---------------------------------------------------------------------------

timeVec = [];
dt = [];
outroi = [];

%---------------------------------------------------------------------------
% Data Loading
%---------------------------------------------------------------------------

if (isnumeric(filename))
	
	data = filename;
	
else
	
	switch (filename)
		
		case 'test1' % Plane wave test (for Calibration)
			
			% data = [time_dt=1s wave2_500mHz wave3_250mHz wave4_125mHz wave5_62.5mHz ... wave12_interoctave-frequency]
			
			OctFreq = 1./(2.^[1:10]); % Octave frequency (in Hz); OctFreq(1)=1/2*dt according to Nyquist
			data = (2*pi*OctFreq)'*[1:1/OctFreq(end)]; % (2*pi*nu*t); timeVec(end)=1/OctFreq(end)
			data = data';
			data = exp(i*data);
			data = real(data); % Cosine
			%data = imag(data); % Sine
			data = [[1:1/OctFreq(end)]',data, real(exp(i*(2*pi*(0.044194)*[1:1/OctFreq(end)])'))];
		
		case 'test2' % Square wave test (for Overshoot estimate)
		
			data = square([1:0.01:31]);
			data = [[1:0.01:31]' data'];
		
		case 'test3' % Edge effect masking test
			
			% data = [time_dt=1s Dirac_delta linear_trend]
			
			data = [1:1024];
			data = [data;[1 zeros(1,1023)];(2/data(end))*data-1]';
		
		case 'test4' % Singularity detection test
			
			% data = [time_dt=1s generic_power Dirac_delta step_function benchmark_test noisy_benchmark]
			
			% Samples n must be a power of 2
			data = [0:1:2047]';
			n = length(data);
			
			a = 0.3;
			data1 = data(n/2)^a-(abs(data(1:end)-data(n/2)).^a);
			
			data2 = [zeros(1,n/2) 10 zeros(1,n/2-1)]';
			
			data3 = [zeros(1,n/2) ones(1,n/2)]';
			
			data4a = data(n/8)^0.5-(abs(data(1:n/4)-data(n/8)).^0.5);
			data4b = data(n/8)^0.4-(abs(data(1:n/4)-data(n/8)).^0.4);
			data4c = [zeros(1,n/8) 10 zeros(1,n/8-1)]';
			data4d = [zeros(1,n/16) 5*ones(1,n/16) fliplr(0:5/(n/16-1):5) zeros(1,n/16)]';
			data4 = [data4a;data4b;data4c;data4d];
			
			amp = 0.1;
			data5 = data4 + 2*amp*rand(n,1)-amp;
			
			%%%%%
			data5 = exp((-(data-1023).^2)/((2*30)^2));
			
			data = [data data1 data2 data3 data4 data5];
			
		case 'test5' % White, Pink and Brown(ian) Noise
			
			n = 8192;
			data = [1:n]';
			data1 = rand(n,1); % Uniform White Noise
			data2 = randn(n,1); % Gaussian White Noise
			data3 = wfbm(0,8192)'; % Pink Noise
			data4 = wfbm(0.5,8192)'; % Brown(ian) Noise (Motion)
			
			data = [data data1 data2 data3 data4];
		
		otherwise % Open experimental data file
			
			data = dlmread(filename);
			
	end

end

%---------------------------------------------------------------------------
% Data Selection and Downsampling
%---------------------------------------------------------------------------

% Data size
[n,ntrace] = size(data);

% Input control
if (ntrace < 2)
	fprintf('\n\nWARNING: Less than 2 data columns!\n');
	fprintf('TimeVec or Data or both are missing...\n');
	fprintf('\n\n');
	return
end

% Trace selection
if (length(roi) == 0)
	roi = [2:ntrace];
end
data = data(:,[1 roi(1:end)]);
outroi = roi;

% DS-Downsampling on the run
if (DS > 1)
	data = data([1:DS:n],:);
end

% Samples must be even
if (mod(n,2) ~= 0)
	data = data(1:n-1,:);
end

% Extract time vector and sampling time
timeVec = data(:,1);
data(:,1) = [];
if (strcmp(unit,'ms')) % Convert to seconds
	timeVec = timeVec/1000;
end
dt = mode(diff(timeVec)); % Sampling time mode (in seconds)

% Accepted time error for each sample
timeError = dt/2;
% Check if time steps are equal 
if ~(isempty(find(abs(diff(timeVec) - dt) > timeError)))
	fprintf('\n\nWARNING: Not equal time steps\n');
	fprintf('Artifacts may be introduced in Wavelet Transform computation!\n');
	fprintf('\n\n');
end

% Rebuild an evenly spaced time vector starting from t=0s
timeVec = ([0:1:length(timeVec)-1]')*dt;


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