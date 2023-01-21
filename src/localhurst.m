function localhurst(x,sig,unit,US,timeVec,freqVec,mark,iter);

%
%--------------------------------------------------------------------------------
% Hurst Exponent (H) Computation
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% hurst(x,sig,unit,US,timeVec,freqVec,mark,iter)
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
% OUTPUT       TYPE       MEANING
% ------       ----       -------
% -none-    -> plot    -> 2 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
[nscale,n] = size(x);
nvoice = nscale/(length(freqVec)-1);

% Scale Vector
scaleVec = 2.^([1:(nscale/nvoice)+1]-log2(US));
scaleVec = fliplr(scaleVec);

% Complete Frequency Vector in Hz
freqVecComp = freqVec(end)*(2.^(-[0:1:nscale-1]/nvoice));
freqVecComp = flipud(freqVecComp');
if (strcmp(unit,'s'))
	freqVecComp = freqVecComp/1000;
end

% log2(W) vs. log2(v)
logfreqVecComp = log2(freqVecComp);
logx = log2(x);
% To prevent log2(0) = -Inf
%logx = log2(x + min(x(find(x))));

%---------------------------------------------------------------------------
% Scalogram Dyadic Split
%---------------------------------------------------------------------------

% Cell Array
a{1} = logx;

h = 1;
while h <= iter
	for j = 1:size(a,2)
		bisec = round(size(a{j},2)/2);
		b{2*j-1} = a{j}(:,1:bisec);
		b{2*j} = a{j}(:,bisec+1:end);		
	end
	a = b;
	h = h+1;
end

% Coherence Test
%test=[];
%for j = 1:size(a,2)
%	size(a{j})
%	test=[test a{j}];
%end
%size(x)
%size(test)
%find(test ~= logx)

%---------------------------------------------------------------------------
% Averaged Wavelet Coefficient (AWC)
%---------------------------------------------------------------------------

figure

set(gca,'Position',[0.12,0.11,0.80,0.68])
xlim([1,nscale])
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(freqVec',4));
ax1 = gca;
hold on

massimo = -Inf;
minimo = +Inf;
coeff{1} = []; % Visible out of for-scope

% Use the highest frequency (the finest scale) for H estimate
k = 3; % Octaves neglected (starting from the Nyquist frequency)
interval = [nscale-round((log2(US)+k+1)*nvoice):nscale-(log2(US)+k)*nvoice];

for j = 1:size(a,2)
	
	AWC = mean(a{j},2);
	if (max(AWC)>massimo)
		massimo = max(AWC);
	end
	if (min(AWC)<minimo)
		minimo = min(AWC);
	end
	
	% Linear fit
	coeff{j} = polyfit(logfreqVecComp(interval),AWC(interval),1);
	coeffgraph{j} = polyfit(interval',AWC(interval),1);
	tang = polyval(coeffgraph{j},[1:nscale]);
	
	hold on
	p1 = plot([1:nscale],AWC,'-b','LineWidth',1);
	p2 = plot([1:nscale],tang,'-r','LineWidth',1);
	
	if (size(a,2)>1)
		% Blue->Red color gradient
		set(p1,'Color',[(j-1)/(size(a,2)-1),0,(size(a,2)-j)/(size(a,2)-1)])
		% Cyan->Green color gradient
		set(p2,'Color',[0,1,(size(a,2)-j)/(size(a,2)-1)])
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

xlim([1,nscale])
ylim(q)
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(scaleVec',4));
xlabel('Scale','FontSize',s2)

% Nyquist limit: the log2(US) highest octaves
hold on
plot([nscale-log2(US)*nvoice,nscale-log2(US)*nvoice],[min(ylim),max(ylim)],'-k','LineWidth',1)

title(['Averaged Wavelet Coefficient'],'FontSize',s3)

%---------------------------------------------------------------------------
% Local Hurst Exponent (H)
%---------------------------------------------------------------------------

ascissa1 = 0;
ascissa2 = 0;
localH = [];

for j = 1:size(a,2)
	
	ascissa1 = ascissa2 + 1;
	ascissa2 = ascissa2 + size(a{j},2);
	
	localH = [localH;(timeVec(ascissa1)+timeVec(ascissa2))/2 -coeff{j}(1)];

end

if (size(a,2)>1)
	
	figure
	
	subplot(2,1,1), hold on
		
		plot(timeVec,sig,'-b','LineWidth',1)
		
		xlim([timeVec(1),timeVec(end)])
		
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		ylabel('(nominal Volt)','FontSize',s2)
		title(['Local Hurst Exponent'],'FontSize',s3)
		
		for l = 1:length(mark)
			plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
		end
	
	subplot(2,1,2), hold on
		
		plot(localH(:,1),localH(:,2),'-b','LineWidth',1)
		
		xlim([timeVec(1),timeVec(end)])
		
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		ylabel('Local H','FontSize',s2)
		xlabel('Time (s)','FontSize',s2)
		
		for l = 1:length(mark)
			plot([timeVec(mark(l)),timeVec(mark(l))],[min(ylim),max(ylim)],'-r','LineWidth',1)
		end

end

% Information about Slices and Points per slice
if (size(a,2)>1)
	
	sli = 2^iter;
	pps = round(n/(2^iter));
	
	fprintf('\n\nLocal Hurst Exponent is computed over %d slices\n',sli);
	fprintf('\nEach slice is made up of about %d points\n',pps);
	fprintf('\n\n');

end

H = mean(localH(:,2));
devH = std(localH(:,2));
semH = devH/sqrt(2^iter);

fprintf('\n\nHurst Exponent\n');
fprintf('\n\n H = %.4f\n',H);
fprintf('\n\n');




%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
if 0
% Remove DC offset and filter
sig = sig - mean(sig);

%sig = diff(sig);
[autoc,lags] = xcorr(sig,'coeff');
hold on
plot(lags,autoc)
return

end
%%%%%%%%%%%%%%%%%%%


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