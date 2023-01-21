function FEAT(wt,par,ratiomat,output)

%
%--------------------------------------------------------------------------------
% Multitrace Features Extraction
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% FEAT(wt,par,ratiomat,output)
%
% INPUT       TYPE         MEANING
% -----       ----         -------
% wt       -> matrix    -> 1st WT Output - Continuous Wavelet Transform
% par      -> structure -> 2nd WT Output - Parameters
% ratiomat -> matrix    -> [PDOutput1;PDOutput2;PDOutput3;...]
% output   -> boolean   -> Print .eps Output Graphs
%
% OUTPUT      TYPE         MEANING
% ------      ----         -------
% -none-   -> plot      -> 1 Plot Resulting from Analysis
%
% WARNING: This code represents a temporary implementation of FEAT function
% ======== because it does not allow comparison among traces acquired under
%          different recording parameters (sampling time, length,...).
%

%---------------------------------------------------------------------------
% Graphic Parameters
%---------------------------------------------------------------------------

s1 = 16; % X-Y TickLabel size
s2 = 19; % X-Y Label and text size
s3 = 24; % Title size

%---------------------------------------------------------------------------
% Variable Assignment
%---------------------------------------------------------------------------

timeVec = par.x;
freqVec = par.y;
periodVec = par.z;
lowest = freqVec(end)/(2^(length(freqVec)-1));
mark = par.m;

%x = abs(wt);
%y = abs(real(wt));

nscale = size(wt,1);
nvoice = size(wt,1)/(length(freqVec)-1);
n = size(wt,2);

e = exp(1); % Essential just for MATLAB

%---------------------------------------------------------------------------
% Input Control
%---------------------------------------------------------------------------

% Print plots - Default value = false = Print nothing
if (nargin < 4)
	output = false;
end

%---------------------------------------------------------------------------
% Multitrace Features Extraction
%---------------------------------------------------------------------------

meanratio = mean(ratiomat);
stdevratio = std(ratiomat);
stderratio = stdevratio./sqrt(size(ratiomat,1));

% Upper and Lower 95% Confidence Interval bounds
% Under the hypothesis of log-normal data distribution
mu = log(meanratio./sqrt((stderratio.^2)./(meanratio.^2)+1));
sigma = sqrt(log((stderratio.^2)./(meanratio.^2)+1));

upci = e.^(mu + 1.96*sigma);
lowci = e.^(mu - 1.96*sigma);

% Under the hypothesis of normal data distribution
%upci = meanratio + 1.96*stderratio;
%lowci = meanratio - 1.96*stderratio;

%---------------------------------------------------------------------------
% Plot
%---------------------------------------------------------------------------

figure

plot([1:nscale],meanratio,'-b','LineWidth',2), hold on
plot([1:nscale],ones(1,nscale),'-g','LineWidth',2)
%plot([1:nscale],upci,'-r','LineWidth',1)
plot([1:nscale],lowci,'-r','LineWidth',1)

xlim([1,nscale])
q = ylim;
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(freqVec'));
xlabel('Frequency (mHz)','FontSize',s2)
ylabel('post / pre Ratio','FontSize',s2)

% Resize -> Syntax Template: set(gca,'Position',[left bottom width height])
set(gca,'Position',[0.12,0.11,0.80,0.68],'Color','none')
axes('Position',[0.12,0.11,0.80,0.68],'XAxisLocation','top','YAxisLocation','left','Color','none');
xlim([1,nscale])
ylim(q)
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(periodVec'));
set(gca,'YTick',[]);
xlabel('Period (s)','FontSize',s2)
title(['Multitrace Features Extraction'],'FontSize',s3)

% Print output graph
if (output)
	print -depsc featureout.eps
end


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