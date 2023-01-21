function resampler(filename,sub,to,unit,roi,output)

%
%--------------------------------------------------------------------------------
% Trace Resampler
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% resampler(filename,sub,to,unit,roi,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% filename -> string  -> 'filename.csv' -> File Name
% sub      -> scalar  -> 3              -> Traces per Plot
% to       -> scalar  -> 5              -> Final Sampling Time
% unit     -> string  -> 's'            -> Time Unit: 's' or 'ms'
% roi      -> array   -> [2,3,9]        -> ROI Set ([x:y] = from x to y)
% output   -> boolean -> true           -> Print .eps Output Graphs
%
% OUTPUT      TYPE                         MEANING
% ------      ----                         -------
% -none-   -> plot                      -> 1 Plot - Sampling Quality Control
% -none-   -> plot                      -> N Plots - Trace Visualization
%

% Graphic parameters
s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

% Output print control - Default value = false = Print nothing
if (nargin < 6)
	output = false;
end
% ROI control - Default value = All traces
if (nargin < 5)
	roi = [];
end
% Time unit control - Default value = 's'
if (nargin < 4)
	unit = 's';
end
if ~(strcmp(unit,'s') | strcmp(unit,'ms'))
	fprintf('\n\nWARNING: Invalid Time Unit\n');
	fprintf('\n\n');
	return
end
% Final sampling time control - Default value = 5s
if (nargin < 3)
	to = 5;
end
% Subplot control - Default value = 1 plot per window
if (nargin < 2)
	sub = 1;
end

% Open data file
data = dlmread(filename);

% Extract time vector
timeVec = data(:,1);
timeVec = timeVec - timeVec(1); % Start from t=0s
if (strcmp(unit,'ms')) % Measured in s
	timeVec = timeVec/1000;
	to = to/1000;
end
dt = mode(diff(timeVec)); % Starting sampling time mode

% Accepted time error in seconds for each sample
timeError = 0.5;

% Check if time steps are equal (evenly spaced)
if ~(isempty(find(abs(diff(timeVec) - dt) > timeError)))
	fprintf('\n\nWARNING: Not equal time steps\nArtifacts may be introduced in Wavelet Transform computation!\n');
end

fprintf('\n\nStarting Parameters\n');
fprintf('===================\n');
fprintf('\n\nSampling Time (Mode) = %.2f s\n',dt); % "%.2f" stands for "floating point number with only 2 decimal places"

% Sampling quality control
figure
plot(diff(timeVec),'-b','LineWidth',1), hold on

xlim([1,length(timeVec)-1])
set(gca,'FontSize',s1,'XTick',unique([1,get(gca,'XTick'),length(timeVec)-1]));
xlabel('Time Step','FontSize',s2)
ylabel('dt (s)','FontSize',s2)
title('Sampling Quality Control','FontSize',s3)
text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['Sampling Time (Mode) = ',num2str(dt),' s'],'FontSize',s2,'Color','k')

% Print other generic information
fprintf('\n\nNyquist Frequency = %.2f mHz\n',1000/(2*dt));
fprintf('\n\nTime Vector Length = %d\n',length(timeVec));
fprintf('\n\nLowest Frequency = %.2f mHz\n',1000/timeVec(end));
fprintf('\n\nOctaves Number = %d\n',floor(log2(length(timeVec))));

% Data matrix
data(:,1) = [];

% Select data of interest, keeping into account ROI 1 is the background
if (length(roi) == 0)
	roi = [2:size(data,2)+1];
end
data = data(:,roi(1:end)-1);

% Check user input
if (size(data,2) < 1)
	fprintf('\n\nWARNING: This file appears not to have any calcium recordings\n');
	fprintf('\n\n');
	return; %Function breaks
end

% Resample
[resdata,H]=resample(data,dt,to);
[restimeVec,H]=resample(timeVec,dt,to);

fprintf('\n\n');
fprintf('\n\nFinal Parameters\n');
fprintf('================\n');
resdt = round(diff(restimeVec)*100)/100;
fprintf('\n\nSampling Time (Mode) = %.2f s\n',mode(resdt));
fprintf('\n\nSampling Time (Mean) = %.2f s\n',mean(resdt));

% Check resampling on time vector
if (abs(mode(resdt) - to) > timeError)
	fprintf('\n\nWARNING: Bad time resampling!\n');
end

fprintf('\n\nTime Remapping = %d s\n',to);

restimeVec = [0:to:timeVec(end)+to];
restimeVec = [restimeVec(1:size(resdata,1))]';

% Print other generic information
fprintf('\n\nNyquist Frequency = %.2f mHz\n',1000/(2*to));
fprintf('\n\nTime Vector Length = %d\n',length(restimeVec));
fprintf('\n\nLowest Frequency = %.2f mHz\n',1000/restimeVec(end));
fprintf('\n\nOctaves Number = %d\n',floor(log2(length(restimeVec))));

% Write new data file
dlmwrite('newsampledata.csv',[restimeVec,resdata])

% Plot original traces of calcium recording
if (sub == 1)
	
	for k = 1:size(data,2)
		
		y = data(:,k);
		resy = resdata(:,k);
		
		figure
		plot(timeVec,y,'-b','LineWidth',1), hold on
		plot(restimeVec,resy,'-r','LineWidth',1)
		
		xlim([restimeVec(1),restimeVec(end)])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(restimeVec(end))]));
		xlabel('Time (s)','FontSize',s2)
		ylabel('Ratio','FontSize',s2)
		title('Trace Resampling','FontSize',s3)
		text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ',num2str(roi(k))],'FontSize',s2,'Color','k')
		
		% Print output graph
		if (output)
			print(['traceout',num2str(k),'.eps'],'-depsc') %Progressive enumeration of output files
		end
	
	end
	
else
	
	for k = 1:size(data,2)
		
		y = data(:,k);
		resy = resdata(:,k);
		
		if (mod(k,sub) == 1)
			figure
			count = 1;
		end
		
		subplot(sub,1,count), hold on
		plot(timeVec,y,'-b','LineWidth',1)
		plot(restimeVec,resy,'-r','LineWidth',1)
		
		xlim([restimeVec(1),restimeVec(end)])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(restimeVec(end))]));
		ylabel('Ratio','FontSize',s2)
		text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ',num2str(roi(k))],'FontSize',s2,'Color','k')
		
		if (count == sub || k == size(data,2))
			xlabel('Time (s)','FontSize',s2)
			subplot(sub,1,1)
			title('Traces Resampling','FontSize',s3)
			
			% Print output graph
			if (output)
				print(['traceout',num2str(floor(k/sub)+(mod(k,sub)~=0)),'.eps'],'-depsc') %Progressive enumeration of output files
			end
		end
		
		count = count+1;
		
	end
	
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