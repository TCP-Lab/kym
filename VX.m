function VX(filename,sub,unit,mark,roi,output)

%
%--------------------------------------------------------------------------------
% Trace Visualization
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% VX(filename,sub,unit,mark,roi,output)
%
% INPUT       TYPE       EXAMPLE           MEANING
% -----       ----       -------           -------
% filename -> string  -> 'filename.csv' -> File Name
% sub      -> scalar  -> 3              -> Traces per Plot (Integer,'S','M','C')
% unit     -> string  -> 's'            -> Time Unit: 's' or 'ms'
% mark     -> array   -> [300,450,500]  -> Time Marker Set (Step Number)
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
% Marker array control - Default value = [] = No markers
if (nargin < 4)
	mark = [];
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
end
dt = mode(diff(timeVec)); % Sampling time mode

% Accepted time error in seconds for each sample
timeError = 0.5;

% Check if time steps are equal (evenly spaced)
if ~(isempty(find(abs(diff(timeVec) - dt) > timeError)))
	fprintf('\n\nWARNING: Not equal time steps\nArtifacts may be introduced in Wavelet Transform computation!\n');
end

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

for k = 1:length(mark) % If length(mark)==0 this does nothing
	plot([mark(k),mark(k)],[min(ylim),max(ylim)],'-r','LineWidth',1)
end

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

% Plot original traces of calcium recording
if (ischar(sub) == 0)
	
	if (sub == 1)
		
		for k = 1:size(data,2)
			
			y = data(:,k);
			
			figure
			plot(timeVec,y,'-b','LineWidth',1), hold on
			
			xlim([timeVec(1),timeVec(end)])
			set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
			xlabel('Time (s)','FontSize',s2)
			ylabel('Ratio','FontSize',s2)
			title('Trace Visualization','FontSize',s3)
			text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ',num2str(roi(k))],'FontSize',s2,'Color','k')
			
			for l = 1:length(mark)
				plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			end
			
			% Print output graph
			if (output)
				print(['traceout',num2str(k),'.eps'],'-depsc') %Progressive enumeration of output files
			end
		
		end
		
	else
		
		for k = 1:size(data,2)
			
			y = data(:,k);
			
			if (mod(k,sub) == 1)
				figure
				count = 1;
			end
			
			subplot(sub,1,count), hold on
			plot(timeVec,y,'-b','LineWidth',1)
			
			xlim([timeVec(1),timeVec(end)])
			set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
			ylabel('Ratio','FontSize',s2)
			text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI ',num2str(roi(k))],'FontSize',s2,'Color','k')
			
			for l = 1:length(mark)
				plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			end
			
			if (count == sub || k == size(data,2))
				xlabel('Time (s)','FontSize',s2)
				subplot(sub,1,1)
				title('Traces Visualization','FontSize',s3)
				
				% Print output graph
				if (output)
					print(['traceout',num2str(floor(k/sub)+(mod(k,sub)~=0)),'.eps'],'-depsc') %Progressive enumeration of output files
				end
			end
			
			count = count+1;
			
		end
		
	end

else
	
	if (size(data,2) < 2)
		fprintf('\n\nWARNING: This visualization option requires more than one calcium trace\n');
		fprintf('\n\n');
		return; %Function breaks
	else
		
		switch (sub)
			
			case 'S' % Superimposition: Plot all traces together with a Blue->Red color gradient
				
				figure
				
				for k = 1:size(data,2)
				
					y = data(:,k);
					
					hold on
					p = plot(timeVec,y,'-b','LineWidth',1);
					
					% Blue->Red color gradient
					set(p,'Color',[(k-1)/(size(data,2)-1),0,(size(data,2)-k)/(size(data,2)-1)])
					
				end
				
				xlim([timeVec(1),timeVec(end)])
				set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
				xlabel('Time (s)','FontSize',s2)
				ylabel('Ratio','FontSize',s2)
				title('Traces Visualization','FontSize',s3)
				text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI from ',num2str(roi(1)),' to ',num2str(roi(end))],'FontSize',s2,'Color','k')
				
				for l = 1:length(mark)
					plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
				end
				
				% Print output graph
				if (output)
					print -depsc traceout.eps
				end
			
			case 'M' % Mean: Plot a single trace representing the mean value +/- SEM at each point in time
				
				figure
				
				y = mean(data,2);
				sem = std(data,0,2)/sqrt(size(data,2));
				
				hold on
				plot(timeVec,y,'-r','LineWidth',2);
				plot(timeVec,y+sem,'-b','LineWidth',1);
				plot(timeVec,y-sem,'-b','LineWidth',1);
				
				xlim([timeVec(1),timeVec(end)])
				set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
				xlabel('Time (s)','FontSize',s2)
				ylabel('Ratio','FontSize',s2)
				title('Traces Visualization','FontSize',s3)
				text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI from ',num2str(roi(1)),' to ',num2str(roi(end))],'FontSize',s2,'Color','k')
				
				for l = 1:length(mark)
					plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
				end
				
				% Print output graph
				if (output)
					print -depsc traceout.eps
				end
			
			case 'C' % Colormap: Plot all traces together with a Blue->Red color gradient and a colormap (x,y,z)=(time,roi,ratio)
				
				figure
				subplot(2,2,1), hold on
					
					for k = 1:size(data,2)
				
						y = data(:,k);
					
						hold on
						p = plot(timeVec,y,'-b','LineWidth',1);
					
						% Blue->Red color gradient
						set(p,'Color',[(k-1)/(size(data,2)-1),0,(size(data,2)-k)/(size(data,2)-1)])
					
					end
					
					xlim([timeVec(1),timeVec(end)])
					set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
					ylabel('Ratio','FontSize',s2)
					title('Traces Visualization','FontSize',s3)
					text(max(xlim)*(3/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),['ROI from ',num2str(roi(1)),' to ',num2str(roi(end))],'FontSize',s2,'Color','k')
					
					for l = 1:length(mark)
						plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
					end
				
				subplot(2,2,3), hold on
					
					imagesc(timeVec(1):timeVec(end),[1:size(data,2)],data')
					xlim([timeVec(1),timeVec(end)])
					ylim([0.5,size(data,2)+0.5])
					set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
					if (size(data,2) >= 5)
						set(gca,'YDir','normal','YTick',unique([1:floor(size(data,2)/5):size(data,2),size(data,2)]));
						set(gca,'YTickLabel',num2str([roi(1:floor(size(data,2)/5):end),roi(end)]'));
					else
						set(gca,'YDir','normal','YTick',[1:1:size(data,2)]);
						set(gca,'YTickLabel',num2str(roi'));
					end
					ylabel('ROI number','FontSize',s2)
					xlabel('Time (s)','FontSize',s2)
				
				% Add color bar
				subplot(2,2,[2 4]), hold on
					
					t = [0:1:100];
					imagesc([1:2],[1:101],t') % X-support [1:2] prevents a "division by zero" in OCTAVE
					xlim([0.5,2.5])
					ylim([0.5,101.5])
					set(gca,'FontSize',s1,'XTick',[0]);
					set(gca,'XTickLabel','');
					set(gca,'YDir','normal','YTick',[1:10:101]);
					massimo = max(max(data));
					minimo = min(min(data));
					magnVec = floor([minimo:(massimo-minimo)/10:massimo]*100)/100; % In order to get only 2 decimal digits
					set(gca,'YTickLabel',num2str(magnVec'));
					xlabel('Ratio','FontSize',s2)
				
				% Resize -> Sintax template: set(gca,'Position',[left bottom width height])
				set(gca,'Position',[0.96,0.11,0.02,0.815])
				
				subplot(2,2,1)
				set(gca,'Position',[0.12,0.58384,0.77,0.34116])
				
				subplot(2,2,3)
				set(gca,'Position',[0.12,0.11,0.77,0.34116])
				
				% Print output graph
				if (output)
					print -depsc traceout.eps
				end
			
			otherwise
				
				fprintf('\n\nWARNING: Invalid subplot parameter\n');
			
		end
		
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