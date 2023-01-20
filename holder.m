function holder(x,sig,unit,US,timeVec,freqVec,periodVec,skex)

%
%--------------------------------------------------------------------------------
% Holder Exponent (alpha) Computation through WTMM
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% holder(x,sig,unit,US,timeVec,freqVec,periodVec,skex)
%
% INPUT        TYPE       MEANING
% -----        ----       -------
% x         -> matrix  -> Continuous Wavelet Modulus
% sig       -> array   -> 3rd WT Output - Original Signal
% unit      -> string  -> Time Unit: 's' or 'ms'
% US        -> scalar  -> US-Upsampling factor
% timeVec   -> array   -> Time Vector
% freqVec   -> array   -> Frequency Vector
% mark      -> array   -> Time Marker Set (Step Number)
% iter      -> scalar  -> Number of dyadic splitting
%
% OUTPUT       TYPE        MEANING
% ------       ----        -------
% maxima    -> matrix   -> Thresholded Maxima Map
% maximaNT  -> matrix   -> Non-Thresholded Maxima Map
% -none-    -> plot     -> 1 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
[nscale,n] = size(x);
nvoice = nscale/(length(freqVec)-1);

sea = nscale-log2(US)*nvoice; % WTMM are like rivers, while sea is the edge of the finest scale
dimfreqVec = freqVec(1,1:end-log2(US));

% Scale Vector
scaleVec = 2.^([1:(sea/nvoice)+1]);
scaleVec = fliplr(scaleVec);

% Complete Frequency Vector in Hz
freqVecComp = dimfreqVec(end)*(2.^(-[0:1:sea-1]/nvoice));
if (strcmp(unit,'s'))
	freqVecComp = freqVecComp/1000;
end

% log2(W) vs. log2(v)
logfreqVecComp = log2(freqVecComp);
logx = log2(x);
% To prevent log2(0) = -Inf
%logx = log2(x + min(x(find(x))));

%---------------------------------------------------------------------------
% Up the rivers!
%---------------------------------------------------------------------------

% Find the river mouths
outlet = find(skex(sea,:));

% River reconstruction for comparison and coerence test
reco = zeros(nscale,n);
reco(sea,outlet) = 1;

% River curve parameter - It must be odd
l = 3;
hl = (l-1)/2;

% Create empty cell arrays
c = {};
d = {};

for h = 1:length(outlet)
	
	c{1,h}(1)=logx(sea,outlet(h));
	
	river = outlet(h);
	upstep = 1;
	
	while (river > hl && river <= n-hl && upstep < sea)
	
		index = find(skex(sea-upstep,river-hl:river+hl));
		
		if (length(index)==1)
		
			river = river+index-(hl+1);
			c{1,h}(upstep+1) = logx(sea-upstep,river);
			reco(sea-upstep,river) = 1;
			
			upstep = upstep+1;
		
		else % length(index)>1->fork --- length(index)==0->source or too curvy
		
			river = 0; % To exit while-loop
		
		end
	
	end
	
	c{2,h} = upstep; % River length

end

% Comparison and coerence test
figure
imagesc(timeVec,1:nscale,skex)
xlim([timeVec(1),timeVec(end)])
ylim([0.5,nscale+0.5])
set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'YTickLabel',num2str(freqVec',4));
if (strcmp(unit,'ms'))
	ylabel('Frequency (Hz)','FontSize',s2)
else
	ylabel('Frequency (mHz)','FontSize',s2)
end
xlabel('Time (s)','FontSize',s2)

figure
imagesc(timeVec,1:nscale,reco)
xlim([timeVec(1),timeVec(end)])
ylim([0.5,nscale+0.5])
set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
set(gca,'YDir','normal','YTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'YTickLabel',num2str(freqVec',4));
if (strcmp(unit,'ms'))
	ylabel('Frequency (Hz)','FontSize',s2)
else
	ylabel('Frequency (mHz)','FontSize',s2)
end
xlabel('Time (s)','FontSize',s2)



massimo = -Inf;
minimo = +Inf;
coeff{1} = []; % Visible out of for-scope

% Use the highest frequency (the finest scale) for alpha estimate
k = 3; % Octaves neglected (starting from the sea == Nyquist frequency)
interval = [k*nvoice+1:(k+1)*nvoice];
%interval = [1:3];

count = 1;
for h = 1:length(outlet)
	if (c{2,h} >= interval(end))
		d{count} = c{1,h};
		count = count+1;
	end
end

%size(c)
%size(c{1,1})
%size(c{1,2})
%size(c{1,3})
%size(c{1,4})
%size(d)

figure

set(gca,'Position',[0.12,0.11,0.80,0.68])
xlim([1,sea])
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:sea])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(dimfreqVec',4));
ax1 = gca;
hold on

for h=1:size(d,2)
	
	if (max(d{h})>massimo)
		massimo = max(d{h});
	end
	if (min(d{h})<minimo)
		minimo = min(d{h});
	end
	
	% Linear fit
	coeff{h} = polyfit(logfreqVecComp(interval),d{h}(interval),1);
	coeffgraph{h} = polyfit(interval,d{h}(interval),1);
	tang = polyval(coeffgraph{h},[1:sea]);
	
	hold on
	p1 = plot([sea-length(d{h})+1:sea],fliplr(d{h}),'-b','LineWidth',1);
	p2 = plot([1:sea],fliplr(tang),'-r','LineWidth',1);
	
	if (size(d,2)>1)
		% Blue->Red color gradient
		set(p1,'Color',[(h-1)/(size(d,2)-1),0,(size(d,2)-h)/(size(d,2)-1)])
		% Cyan->Green color gradient
		set(p2,'Color',[0,1,(size(d,2)-h)/(size(d,2)-1)])
	end

end

ylim([round(minimo-1),round(massimo+1)])
q = ylim;
if (strcmp(unit,'ms'))
	xlabel('Frequency (Hz)','FontSize',s2)
else
	xlabel('Frequency (mHz)','FontSize',s2)
end
ylabel('Log_2(W)','FontSize',s2)

ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','YTick',[],'YTickLabel',[],'Color','none','XColor','k','YColor','k');
linkaxes([ax1 ax2],'x');

xlim([1,sea])
ylim(q)
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:sea])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(scaleVec',4));
xlabel('Scale','FontSize',s2)

title(['Holder Exponent \alpha Computation'],'FontSize',s3)

for j=1:size(coeff,2)
-coeff{j}(1)
end




return



%---------------------------------------------------------------------------
% Compute Holder exponent (alpha) and plot results
%---------------------------------------------------------------------------

figure

subplot(3,3,[1 2 3]), hold on
	
	plot(timeVec,sig,'-b','LineWidth',1)
	
	xlim([timeVec(1),timeVec(end)])
	
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Ratio','FontSize',s2)
	%title(['Holder Exponents - ROI ',num2str(par.r)],'FontSize',s3)
	
	%for l = 1:length(mark)
	%	plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	%end

subplot(3,3,[4 5 6]), hold on
	
	ascissa1 = 0;
	ascissa2 = 0;
	localalpha = [];
	
	for j = 1:size(a,2)
		
		ascissa1 = ascissa2 + 1;
		ascissa2 = ascissa2 + size(a{j},2);
		
		localalpha = [localalpha;(timeVec(ascissa1)+timeVec(ascissa2))/2 -coeff{j}(1)];
		
	end
	
	plot(localalpha(:,1),localalpha(:,2),'-b','LineWidth',1)
	
	xlim([timeVec(1),timeVec(end)])
	
	set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
	ylabel('Holder \alpha','FontSize',s2)

	%for l = 1:length(mark)
	%	plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
	%end
	
subplot(3,3,8), hold on

	%hist(localalpha(:,2),[-1:0.1:1],1)
	hist(localalpha(:,2),10)

mean(localalpha(:,2))
std(localalpha(:,2))

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