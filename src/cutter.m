function cutter(filename,bounds,unit,roi,output)

%
%--------------------------------------------------------------------------------
% Trace Edge Cutter
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% cutter(filename,bounds,unit,roi,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% filename -> string  -> 'filename.csv' -> File Name
% bounds   -> array   -> [100,3000]     -> Lower and Upper Bound, respectively
% unit     -> string  -> 's'            -> Bound Time Unit: 's' or 'ms'
% roi      -> array   -> [2,3,9]        -> ROIs TO SHOW ([x:y] = from x to y) - NOT IMPLEMENTED!
% output   -> boolean -> true           -> Print .eps Output Graphs
%
% OUTPUT      TYPE                         MEANING
% ------      ----                         -------
% -none-   -> plot                      -> 1 Plot - Sampling Quality Control
%

% Graphic parameters
s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

% Output print control - Default value = false = Print nothing
if (nargin < 5)
	output = false;
end
% ROI visualization control - Default value = Show all traces
if (nargin < 4)
	roi = [];
end
% Time unit control - Default value = 's'
if (nargin < 3)
	unit = 's';
end
if ~(strcmp(unit,'s') | strcmp(unit,'ms'))
	fprintf('\n\nWARNING: Invalid Time Unit\n');
	fprintf('\n\n');
	return
end
% Cutting bound control - Default value = No bounds
if (nargin < 2)
	bounds = [];
end

% Open data file
data = dlmread(filename);

% Extract time vector
timeVec = data(:,1);
timeVec = timeVec - timeVec(1); % Start from t=0s
if (strcmp(unit,'ms')) % Measured in s
	timeVec = timeVec/1000;
end
dt = mode(diff(timeVec)); % Sampling time mode

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

% Select data to show, keeping into account ROI 1 is the background
if (length(roi) == 0)
	roi = [2:size(data,2)+1];
end

% Check user input
if (size(data,2) < 1)
	fprintf('\n\nWARNING: This file appears not to have any calcium recordings\n');
	fprintf('\n\n');
	return; %Function breaks
end

% Cut according to bounds
low = bounds(1);
high = bounds(2);

[useless1 lowindex] = min(abs(timeVec-low));
[useless2 highindex] = min(abs(timeVec-high));

cutdata = data(lowindex:highindex,:);
cuttimeVec = timeVec(lowindex:highindex,1);
cuttimeVec = cuttimeVec - cuttimeVec(1); % Restart from t=0s

% Print other generic information
fprintf('\n\n');
fprintf('\n\nFinal Parameters\n');
fprintf('================\n');
fprintf('\n\nSampling Time (Mode) = %.2f s\n',mode(diff(cuttimeVec)));
fprintf('\n\nNyquist Frequency = %.2f mHz\n',1000/(2*mode(diff(cuttimeVec))));
fprintf('\n\nTime Vector Length = %d\n',length(cuttimeVec));
fprintf('\n\nLowest Frequency = %.2f mHz\n',1000/cuttimeVec(end));
fprintf('\n\nOctaves Number = %d\n',floor(log2(length(cuttimeVec))));

% Write new data file
dlmwrite('cutdata.csv',[cuttimeVec,cutdata])

% Plot original traces of calcium recording
figure

subplot(2,1,1), hold on
	
	for k = 1:size(data,2)
		
		y = data(:,k);
		
		p = plot(timeVec,y,'-b','LineWidth',1);
		
		% Blue->Red color gradient
		set(p,'Color',[(k-1)/(size(data,2)-1),0,(size(data,2)-k)/(size(data,2)-1)])
	
	end
	
	xlim([timeVec(1),timeVec(end)])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Ratio','FontSize',s2)
	title('Trace Edge Cutting','FontSize',s3)
	text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI from ',num2str(roi(1)),' to ',num2str(roi(end))],'FontSize',s2,'Color','k')
	
	plot([low,low],[min(ylim),max(ylim)],'-r','LineWidth',1)
	plot([high,high],[min(ylim),max(ylim)],'-r','LineWidth',1)

subplot(2,1,2), hold on
	
	for k = 1:size(cutdata,2)
		
		cuty = cutdata(:,k);
		
		p = plot(cuttimeVec,cuty,'-b','LineWidth',1);
		
		% Blue->Red color gradient
		set(p,'Color',[(k-1)/(size(cutdata,2)-1),0,(size(cutdata,2)-k)/(size(cutdata,2)-1)])
	
	end
	
	xlim([cuttimeVec(1),cuttimeVec(end)])
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(cuttimeVec(end))]));
	xlabel('Time (s)','FontSize',s2)
	ylabel('Ratio','FontSize',s2)

% Print output graph
if (output)
	print -depsc traceout.eps
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