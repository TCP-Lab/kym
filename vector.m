function vratio = vector(x,coimask,freqVec,periodVec,mark,roi,output)

%
%--------------------------------------------------------------------------------
% CWT Vector Approach
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% vratio = vector(x,coimask,freqVec,periodVec,mark,roi,output)
%
% INPUT        TYPE        MEANING
% -----        ----        -------
% x         -> matrix   -> Continuous Wavelet Modulus
% coimask   -> matrix   -> COI Mask
% freqVec   -> array    -> Frequency Vector
% periodVec -> array    -> Period Vector
% mark      -> array    -> Time Marker Set (Step Number)
% roi       -> scalar   -> Current ROI
% output    -> boolean  -> Print .eps Output Graphs
%
% OUTPUT       TYPE        MEANING
% ------       ----        -------
% vratio    -> array    -> post/pre Spectrum Ratio Vector
% -none-    -> plot     -> 2 Plots Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

n = size(x,2);
nscale = size(x,1);
nvoice = (size(x,1))/(length(freqVec)-1);

x = x .* coimask;

v(1,:) = 1000*(sum(x(:,1:mark(1)),2)/(mark(1)))';
if (length(mark) > 1)
	for h = 2:length(mark)
		v(h,:) = 1000*(sum(x(:,mark(h-1):mark(h)),2)/(mark(h)-mark(h-1)+1))';
	end
end
v(length(mark)+1,:) = 1000*(sum(x(:,mark(end):end),2)/(n-mark(end)+1))';

vratio = v(2,:)./v(1,:);
vratio(find(isnan(vratio) | isinf(vratio))) = 0; % Prevent to NaN and Inf from being plotted

% Plotting
figure

for k = 1:(length(mark)+1)
	
	hold on
	p = plot([1:nscale],v(k,:),'-b','LineWidth',1);
	
	% Blue->Red color gradient
	col(k,:) = [(k-1)/length(mark),0,(length(mark)+1-k)/length(mark)]; % Color code
	set(p,'Color',col(k,:))
	lgn{k} = ['W_',num2str(k)]; % Cell Arrays
	
end

xlim([1,nscale])
q = ylim;
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(freqVec'));
xlabel('Frequency (mHz)','FontSize',s2)
ylabel('10^3 V Amplitude','FontSize',s2)

mylegend = legend(lgn,'Location','East');
set(mylegend,'FontSize',s1);

% Resize -> Syntax Template: set(gca,'Position',[left bottom width height])
set(gca,'Position',[0.12,0.11,0.80,0.68],'Color','none')
axes('Position',[0.12,0.11,0.80,0.68],'XAxisLocation','top','YAxisLocation','left','Color','none');
xlim([1,nscale])
ylim(q)
set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
set(gca,'XTickLabel',num2str(periodVec'));
set(gca,'YTick',[]);
xlabel('Period (s)','FontSize',s2)
title(['Continuous Wavelet Transform Vector Analysis - ROI ',num2str(roi)],'FontSize',s3)

if (length(mark) == 1)
	
	dist = sqrt(sum((v(2,:)-v(1,:)).^2)/(nvoice));
	
	normpre = sqrt(sum(v(1,:).^2)/(nvoice));
	normpost = sqrt(sum(v(2,:).^2)/(nvoice));
	delta = (normpost-normpre);
	
	theta = acos((sum(v(2,:).*v(1,:))/(nvoice))/(normpost*normpre));
	theta = (theta/(pi))*180; % Radians -> Degrees conversion
	
	dist = round(dist*100)/100; % In order to get only 2 decimal digits
	delta = round(delta*100)/100; % In order to get only 2 decimal digits
	theta = round(theta*100)/100; % In order to get only 2 decimal digits
	
	text(max(xlim)*(75/100),max(ylim)-(max(ylim)-min(ylim))*(6/100),['W_2 vs. W_1'],'FontSize',s2,'Color','r')
	text(max(xlim)*(75/100),max(ylim)-(max(ylim)-min(ylim))*(12/100),['D = ',num2str(dist)],'FontSize',s2,'Color','k')
	text(max(xlim)*(75/100),max(ylim)-(max(ylim)-min(ylim))*(18/100),['\Delta = ',num2str(delta)],'FontSize',s2,'Color','k')
	text(max(xlim)*(75/100),max(ylim)-(max(ylim)-min(ylim))*(24/100),['\theta = ',num2str(theta),'^\circ'],'FontSize',s2,'Color','k')
	
	% Print output graph
	if (output)
		print -depsc vectorout.eps
	end
	
else
	
	% Print output graph
	if (output)
		print -depsc vectorout.eps
	end
	
	figure
	
	for h = 0:(length(mark)-1)
		for k = 1:(length(mark)-h)
			
			dist = sqrt(sum((v(end-h,:)-v(end-k-h,:)).^2)/(nvoice));
			
			normpre = sqrt(sum(v(end-k-h,:).^2)/(nvoice));
			normpost = sqrt(sum(v(end-h,:).^2)/(nvoice));
			delta = (normpost-normpre);
			
			theta = acos((sum(v(end-h,:).*v(end-k-h,:))/(nvoice))/(normpost*normpre));
			theta = (theta/(pi))*180; % Radians -> Degrees conversion
			
			dist = round(dist*100)/100; % In order to get just 2 decimal digits
			delta = round(delta*100)/100; % In order to get just 2 decimal digits
			theta = round(theta*100)/100; % In order to get just 2 decimal digits
			
			subplot(length(mark),length(mark),length(mark)*(length(mark)-h)-k+1-h), hold on
			p1 = plot([1:nscale],v(end-h,:),'-r','LineWidth',1);
			set(p1,'Color',col(end-h,:))
			p2 = plot([1:nscale],v(end-k-h,:),'-b','LineWidth',1);
			set(p2,'Color',col(end-k-h,:))
			
			xlim([1,nscale])
			ylim(q) % Fix y-range to enable comparison among different plots
			set(gca,'FontSize',s1,'XTick',unique([1,nvoice:nvoice:nscale])); % Element 1 is needed to display freqVec(1) value
			set(gca,'XTickLabel',num2str(freqVec'));
			
			if (h == 0)
				xlabel('Frequency (mHz)','FontSize',s2)
			end
			if (k == (length(mark)-h))
				ylabel('10^3 V Amplitude','FontSize',s2)
			end
			
			if (h == length(mark)-1)
				title(['Paired Vector Analysis - ROI ',num2str(roi)],'FontSize',18)
			end
			
			% Text Size exception
			text(max(xlim)*(65/100),max(ylim)-(max(ylim)-min(ylim))*(10/100),['W_',num2str(length(mark)+1-h),' vs. W_',num2str(length(mark)+1-k-h)],'FontSize',s1,'Color','r')
			text(max(xlim)*(65/100),max(ylim)-(max(ylim)-min(ylim))*(20/100),['D = ',num2str(dist)],'FontSize',s1,'Color','k')
			text(max(xlim)*(65/100),max(ylim)-(max(ylim)-min(ylim))*(30/100),['\Delta = ',num2str(delta)],'FontSize',s1,'Color','k')
			text(max(xlim)*(65/100),max(ylim)-(max(ylim)-min(ylim))*(40/100),['\theta = ',num2str(theta),'^\circ'],'FontSize',s1,'Color','k')
		
		end
	end
	
	% Print output graph
	if (output)
		print -depsc pairvecout.eps
	end
	
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