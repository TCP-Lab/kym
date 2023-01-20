function [WAI,MEAN] = wai(x,coimask,thr,sig,timeVec,freqVec,mark,maxima,indextype,grafunc)

%
%--------------------------------------------------------------------------------
% CWT Wave-Activity Index - WAI
%--------------------------------------------------------------------------------
%
%
% Function Definition
%
% [WAI,MEAN] = wai(x,timeVec,freqVec,nvoice,mark,maxima,flag,grafunc)
%
% INPUT       TYPE        MEANING
% -----       ----        -------
% x        -> matrix   -> Continuous Wavelet Modulus
% timeVec  -> array    -> Time Vector
% freqVec  -> array    -> Frequency Vector
% nvoice   -> scalar   -> Lines of Pixel per Octave
% mark     -> array    -> Time Marker Set (Step Number)
% maxima   -> matrix   -> paths.m Output - Maxima Coordinates
% flag     -> scalar   -> 0 = Thresholded WAI / 1 = NO Threshold
% grafunc  -> boolean  -> Enable Graphical Functions
%
% OUTPUT      TYPE        MEANING
% ------      ----        -------
% WAI      -> matrix   -> Activity Index Values
% MEAN     -> matrix   -> Time Averaged Activity
% -none-   -> plot     -> 1/2 Plot Resulting from Analysis
%

% Graphic Parameters
s1 = 16; % X-Y TickLabel Size
s2 = 19; % X-Y Label and Text Size
s3 = 24; % Title Size

% Variables Assignment
nscale = size(x,1);
nvoice = (size(x,1))/(length(freqVec)-1);
n = size(x,2);

lowest = (freqVec(end)/1000)/(2^(length(freqVec)-1));

% Frequency matrix for the Complete WAI
freqVecComp = (freqVec(end)/1000).*2.^(-([0:1:nscale-1])./nvoice);
freqVecComp = flipud(freqVecComp');
freqVecComp = freqVecComp*ones(1,n);



% Thresholded index
xthrsld = (x-min(min(x)))/max(max(x-min(min(x))));
threshold = thr/100;
xthrsld = (xthrsld >= threshold);
xthrsld = x .* xthrsld .* coimask;
x = x .* coimask;

% l = Regularization window width - Mean convolutive filter
% l value must be odd! - l=1 means NO FILTER
l = 1;



switch (indextype)
	
	case 'I'
		
		WAI(1,:) = 2*sum(freqVecComp.*(xthrsld.^2),1)*(log(2)/nvoice);
		WAI(2,:) = 2*sum((freqVecComp.*((x.*maxima).^2)),1)*log(2); %??
		
	case 'Y'
		
		WAI(1,:) = 2*sum(log2(freqVecComp/lowest).*freqVecComp.*(xthrsld.^2))*(log(2)/nvoice);
		WAI(2,:) = 2*sum((log2(freqVecComp/lowest).*freqVecComp.*((x.*maxima).^2)))*log(2);
		
	case 'J'
		
		WAI(1,:) = 2*(1/lowest)*sum((freqVecComp.^2).*(xthrsld.^2))*(log(2)/nvoice);
		WAI(2,:) = 2*(1/lowest)*sum(((freqVecComp.^2).*((x.*maxima).^2)))*log(2); %??
		
end

% Smoothing convolution
WAItemp = WAI(2,:);
for h = 1:n-l+1
	WAI(2,h+(l-1)/2) = (1./l)*sum(WAItemp(:,h:h+l-1),2);
end

% Mean values
if (length(mark) > 0)
	MEAN(:,1) = (1/mark(1))*sum(WAI(:,1:mark(1)),2);
	if (length(mark) > 1)
		for h = 2:length(mark)
			MEAN(:,h) = (1/(mark(h)-mark(h-1)+1))*sum(WAI(:,mark(h-1):mark(h)),2);
		end
	end
	MEAN(:,length(mark)+1) = (1/(n-mark(end)+1))*sum(WAI(:,mark(end):end),2);
else
	MEAN(:,1) = (1/size(WAI,2))*sum(WAI,2);
end

% Enable graphical functions
if (grafunc)
	
	subplot(3,1,1), hold on
		
		plot(timeVec,sig,'-b','LineWidth',1)
		
		xlim([timeVec(1),timeVec(end)])
		set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
		ylabel('Ratio','FontSize',s2)
		
		for k = 1:length(mark) % If length(mark)==0 this does nothing
			plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
		end
		
	for g = 1:2
		
		subplot(3,1,1+g), hold on
			
			plot(timeVec,WAI(g,:),'-b','LineWidth',1)
			
			xlim([timeVec(1),timeVec(end)])
			set(gca,'FontSize',s1,'XTick',unique([get(gca,'XTick'),floor(timeVec(end))]));
			
			if (g == 1)	
				ylabel(['Index ' indextype ' Complete'],'FontSize',s2)
			else
				ylabel(['Index ' indextype ' Maxima'],'FontSize',s2)
				xlabel('Time (s)','FontSize',s2)
			end
			
			if (length(mark) > 0)
				text(timeVec(15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(g,1),5)],'FontSize',s2,'Color','k')
				if (length(mark) > 1)
					for h = 2:length(mark)
						text(timeVec(mark(h-1)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(g,h),5)],'FontSize',s2,'Color','k')
					end
				end
				text(timeVec(mark(end)+15),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(g,end),5)],'FontSize',s2,'Color','k')
				
				r = MEAN(g,2)/MEAN(g,1);
				text(max(xlim)+15,max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(r,2)],'FontSize',s2,'Color','r')
			else
				text(max(xlim)*(75/100),max(ylim)-(max(ylim)-min(ylim))*(15/100),[num2str(MEAN(g,1),5)],'FontSize',s2,'Color','k')
			end
			
			for k = 1:length(mark) % If length(mark)==0 this does nothing
				plot([timeVec(mark(k)),timeVec(mark(k))],[min(ylim),max(ylim)],'-r','LineWidth',1)
			end
			
	end
	
	subplot(3,1,1), hold on
	
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