function multiKLprogression(ev, vals, KL, caption, offset, isFinal, normalize)

	if ~exist('normalize')
		normalize = false;
	end
	
	OAKL = 1;
	skipKL4 = true;
   withLines = 1;
		
	% Normalizing the values - First make Std = 1
	if normalize
      okvals = find(~isnan(vals));
		m = mean(vals(okvals));
		s = std(vals(okvals));
      if s<1000*eps
         vals = 0;
      else
         vals = vals-m;
         vals = vals / s;
      end
	end
	
	count = length(vals);
	OA = [];
	Healthy = [];
	KL0 = []; KL1 = []; KL2 = []; KL3 = []; KL4 = [];
	for p = 1:count
		% Check if we have KL value for the knee
		val  = vals(p);
      if isnan(val)
         continue
      end
		if KL(p)<0
			error('No KL!')
		else
			switch KL(p)
				case 0, KL0(end+1) = val;
				case 1, KL1(end+1) = val;
				case 2, KL2(end+1) = val;
				case 3, KL3(end+1) = val;
				case 4
               if skipKL4
                  KL3(end+1) = val;
               else
                  KL4(end+1) = val;
               end                  
			end
			if KL(p) >= OAKL
				OA(end+1) = val;
			else
				Healthy(end+1) = val;
			end
		end
   end
   if ~std(Healthy)
      disp(sprintf(' %27s: All same value',caption))
   else
      [H, P, CI] = ttest2(Healthy, OA);
      [H, P01, CI] = ttest2(KL0, KL1);
      if P01<0.1
         sample01 = SampleSize(KL0, KL1);
         odds01 = OddsRatio(KL0,KL1);
      else
         sample01 = NaN;
         odds01 = NaN;
      end
      disp(sprintf(' %27s: %3d healthy and %3d OA -> pHO %.5f, p01 %.5f (%3.0f, %3.1f)',...
         caption,length(Healthy),length(OA),P,P01,sample01,odds01));
   end
   
	% Normalizing the values - Then make mean of KL0 to be zero
	if normalize
		m = mean(KL0);
		KL0 = KL0-m;
		KL1 = KL1-m;
		KL2 = KL2-m;
		KL3 = KL3-m;
		KL4 = KL4-m;
	end

	if ~length(KL4)
		skipKL4 = true;
	end
	
	persistent idx
	persistent legtext
	if isempty(idx)
		figure
		hold on
		idx = 1;
		legtext = {};
	else
		idx = idx+1;
	end
	
	fs = 16;
   eKL0 = sem(KL0);
   eKL1 = sem(KL1);
   eKL2 = sem(KL2);
   eKL3 = sem(KL3);
   if ~length([KL2,KL3,KL4])
      idxs = 0:1;
      means = [mmean(KL0),mmean(KL1)];
      errors = [eKL0,eKL1];
   elseif skipKL4
      idxs = 0:3;
      means = [mmean(KL0),mmean(KL1),mmean(KL2),mmean(KL3)];
      errors = [eKL0,eKL1,eKL2,eKL3];
   else
      eKL4 = sem(KL4);
      idxs = 0:4;
      means = [mmean(KL0),mmean(KL1),mmean(KL2),mmean(KL3),mmean(KL4)];
      errors = [eKL0,eKL1,eKL2,eKL3,eKL4];
   end
	idxs = idxs + offset;
	% disp(means)
   cols = hsv(18);
   if withLines
      lineSty = '-';
   else
      lineSty = 'none';
   end
	errorbar(idxs,means, errors, 'LineStyle',lineSty,'color',cols(idx,:),'marker','diamond',...
      'linewidth',3,'markersize',6,'markerfacecolor',cols(idx,:));
	legtext{idx} = caption;
	if isFinal
		set(gca,'xtick',0:4,'xticklabel',{'KL=0','KL=1','KL=2','KL=3','KL=4'},'fontsize',fs)
		xlabel(sprintf('KL score (N = %d, %d, %d, %d, %d)',length(KL0),length(KL1),length(KL2),length(KL3),length(KL4)))
      if normalize
         ylabel({'Normalized (zero mean, std one)','Mean and SEM'},'fontsize',fs)
      else
         ylabel({'Y label?','Mean and SEM'},'fontsize',fs)
      end
		set(gca,'xtickmode','manual','ytickmode','manual','xticklabelmode','manual','yticklabelmode','manual')
		xlim([-0.5,4.5])
		set(gca,'yticklabel',fixTickLabels(get(gca,'ytick')));
		legend(legtext,'Location','EastOutside');
		idx = [];
	end
	
	return
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

function s = sem(vals)
if isempty(vals)
   s = 0;
else
   s = std(vals) / sqrt(length(vals));
end

function m = mmean(vals)
if isempty(vals)
   m = 0;
else
   m = mean(vals);
end
