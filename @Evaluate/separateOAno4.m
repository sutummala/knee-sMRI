function P = separateOAno4(ev, knees, vals, KL, caption, makeFigure, OAKL, doKL01)

	doSTD = false; % Alternative is to do SEM
	
	if isempty(KL)
		disp('Finding KL elsewhere...')
		cc = CCBR;
		[localKnee, localKLs, KL2006] = GetKL(cc);
		localKL = 1;
	else
		localKL = 0;
	end
	
	if ~exist('OAKL')
		OAKL = 1;
	end
	if ~exist('doKL01')
		doKL01 = 0;
	end
		
	Healthy = [];
	OA = [];
	OAKLhist = zeros(4,1);
	count = length(vals);
	KL0 = []; KL1 = []; KL2 = []; KL3 = [];
	for p = 1:count
		if localKL
			lk = knees{p};
			lk(end-1) = ' ';
			lk = strmatch(lk,localKnee,'exact');
			if isempty(lk)
				KLp = -1;
			else
				KLp = localKLs(lk);
			end
		else
			KLp = KL(p);
		end
		% Check if we have KL value for the knee
		val  = vals(p);
		if KLp<0
			disp(['No KL for ',knees{p}])
		else
			switch KLp
				case 0, KL0(end+1) = val;
				case 1, KL1(end+1) = val;
				case 2, KL2(end+1) = val;
				case 3, KL3(end+1) = val;
				case 4, KL3(end+1) = val;
			end
			if KLp >= OAKL
				OA(end+1) = val;
				OAKLhist(KLp) = OAKLhist(KLp)+1;
			else
				Healthy(end+1) = val;
			end
		end
	end
	[H, P, CI] = ttest2(Healthy, OA);
	disp(sprintf('SEP %35s: %3d healthy (%f) and %3d OA (%f)(%d,%d,%d,%d)-> p %.6f',...
		caption,length(Healthy),mean(Healthy),length(OA),mean(OA),OAKLhist,P));
	if doKL01
		[H, P01, CI] = ttest2(KL0, KL1);		
		disp(sprintf('%39s  %3d KL 0 (%f) and %3d KL 1 (%f) -> p %.5f',...
			' ',length(KL0),mean(KL0),length(KL1),mean(KL1),P01));
	end

	if makeFigure
		doKLhist = 1;
		figure
		fs = 16;
		eHe = std(Healthy);
		eOA = std(OA);
		if doKLhist
			eKL0 = std(KL0);
			eKL1 = std(KL1);
			eKL2 = std(KL2);
			eKL3 = std(KL3);
		end
		if ~doSTD
			eHe = eHe / sqrt(length(Healthy));
			eOA = eOA / sqrt(length(OA));
			if doKLhist
				eKL0 = eKL0 / sqrt(length(KL0));
				eKL1 = eKL1 / sqrt(length(KL1));
				eKL2 = eKL2 / sqrt(length(KL2));
				eKL3 = eKL3 / sqrt(length(KL3));
			end
		end
		if ~doKLhist
			h = errorbar(1:2,[mean(Healthy),mean(OA)], [eHe,eOA],'LineStyle','none','color','blue','marker','+');
			set(gca,'xtick',1:2,'xticklabel',...
			{'Healthy','OA (KL>=1)';sprintf('(N=%d)',length(Healthy)),sprintf('(N=%d)',length(OA))},'fontsize',fs)
				%	{sprintf('Healthy (N=%d)',length(Healthy)),sprintf('OA (KL>=1, N=%d)',length(OA))},'fontsize',fs)
		else
			idxs = [1:2,4:7];
			means = [mean(Healthy),mean(OA),mean(KL0),mean(KL1),mean(KL2),mean(KL3)];
			errors = [eHe,eOA,eKL0,eKL1,eKL2,eKL3];
			% disp(means)
			if OAKLhist(4)
				tag3 = 'KL=3/4';
			else
				tag3 = 'KL=3';
			end
%			h = errorbar(idxs,means, errors, 'LineStyle','none','color','blue','marker','+');
			h = errorbar(idxs,means, errors, 'LineStyle','none','color','black','marker','+');
			set(gca,'xtick',1:7,'xticklabel',{'Healthy','OA','','KL=0','KL=1','KL=2',tag3},'fontsize',fs)
			xlabel(' ') % make room for extra text
			y = min(ylim)-0.08*(max(ylim)-min(ylim));
			text(1,y,sprintf('N=%d',length(Healthy)),'horizontalalignment','center','fontsize',fs)
			text(2,y,sprintf('N=%d',length(OA)),'horizontalalignment','center','fontsize',fs)
			text(4,y,sprintf('N=%d',length(KL0)),'horizontalalignment','center','fontsize',fs)
			text(5,y,sprintf('N=%d',length(KL1)),'horizontalalignment','center','fontsize',fs)
			text(6,y,sprintf('N=%d',length(KL2)),'horizontalalignment','center','fontsize',fs)
			text(7,y,sprintf('N=%d',length(KL3)),'horizontalalignment','center','fontsize',fs)
			line([3 3],ylim,'color','black','linestyle',':')
		end
		set(h,'linewidth',2,'markersize',16);
		if doSTD
			ylabel({caption,'Mean and Std. Dev.'},'fontsize',fs)
		else
			ylabel({caption,'Mean and SEM'},'fontsize',fs)
		end
		set(gca,'xtickmode','manual','ytickmode','manual','xticklabelmode','manual','yticklabelmode','manual')
		if makeFigure>0
			title('{\bfSeparation of OA from Healthy}','fontsize',fs)
		end
		if ~doKLhist
			xlim([0.5,2.5])
		else
			xlim([0.5,7.5])
		end
		set(gca,'yticklabel',fixTickLabels(get(gca,'ytick')));
		xc = mean(xlim);
		yl = ylim;
		yc = yl(2)-0.06*(yl(2)-yl(1));
		if P<0.05
			verdict = 'Successful';
		else
			verdict = 'Failed';
		end
		% text(xc,yc,sprintf('%s t-test (p = %.5f)',verdict, P),'horizontalalignment','center','fontsize',fs)
		yc = yl(1)+0.06*(yl(2)-yl(1));
		text(1.75,yc,sprintf('p = %.5f', P),'horizontalalignment','center','fontsize',14)
		% axis square
	end
