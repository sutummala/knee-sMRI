function [R,meanRelDif,meanDif,repro,reproKnee] = reproPlot(ev, kneeIds, vals, caption, makeFigure, doWorst)
	if ~exist('doWorst','var')
		doWorst = 0;
	end

	if isempty(kneeIds)
		repro = vals;
	else
		% Split V* rescans from the rest
		count = length(vals);
		isV = zeros(count,1);
		for k = 1:count
			knee = kneeIds{k};
			isV(k) = knee(1) == 'V';
		end
		idxV = find(isV);
		if isempty(idxV)
			disp('*************** No V scans in reproPlot ***************')
			R=Inf;meanRelDif=Inf;meanDif=Inf; return
		end
		idx = find(~isV);
		kneeIdsV = kneeIds(idxV);
		kneeIds = kneeIds(idx);
		valsV = vals(idxV);
		vals = vals(idx);
		
		% Evaluate reproducibility
		countV = length(valsV);
		repro = zeros(countV,2);
		reproKnee = cell(countV,1);
		for s = 1:countV
			repro(s,2) = valsV(s);
			knee = kneeIdsV{s};
			knee = knee(2:end);
			nonV = strmatch(knee,kneeIds,'exact');
			if isempty(nonV)
				error(['Did not find nonV version of ', knee])
			end
			repro(s,1) = vals(nonV);
			reproKnee{s} = knee;
		end
	end

	% Find correlation
	[R,P]=corrcoef(repro);
	R = R(2); P = P(2);
	% Find relative error
	ratio = repro(:,1) ./ repro(:,2);
	ratio = max(ratio, 1./ratio);
	meanRelDif = 100*(mean(ratio)-1);
	% Do it again, different method
	ratio2 = 2*abs(repro(:,1)-repro(:,2)) ./ (repro(:,1)+repro(:,2));
	meanRel2 = mean(ratio2)*100;
	% Mean abs diff
	meanDif = mean(abs(repro(:,1) - repro(:,2)));
	disp(sprintf('REPRO %23s: corr. %.2f, mean dif %4.1f%% / %4.1f%% / %f ',caption,R,meanRelDif,meanRel2,meanDif))
	
	% Find worst
	if doWorst && ~isempty(kneeIds)
		relerror = repro(:,1)./repro(:,2);
		relerror = max(relerror,1./relerror);
		[worst,map] = sort(relerror,1,'descend');
		idx = find(worst>1.15)';
		for i = idx
			disp(sprintf('  Worst Reproducibility: %s - %f (%f vs %f)',kneeIdsV{map(i)}, worst(i), repro(map(i),1), repro(map(i),2)))
		end
	end
	% Find Mean KL for these
	if 0
		[kneeKL, Sex, KL] = getKneeSexKL();
		kls = zeros(5,1);
		for k = 1:length(knee)
			idx = strmatch(knee{k},kneeKL,'exact');
			kls(KL(idx)+1) = kls(KL(idx)+1)+1;
		end
		klmean = 0*kls(1) + 1*kls(2) + 2*kls(3) + 3*kls(4) + 4*kls(5);
		klmean = klmean / sum(kls)
		kls'
	end

	% Make Figure
	if makeFigure
		figure
		fs = 16;
%		scatter(repro(:,1),repro(:,2),'filled','markerfacecolor','black')
		scatter(repro(:,1),repro(:,2),'filled')
		xl = xlim;
		yl = ylim;
		lims = [min(xl(1),yl(1)),max(xl(2),yl(2))];
		xlim(lims); ylim(lims);
		if 0
			set(gca,'xtickmode','manual','ytickmode','manual','xticklabelmode','manual','yticklabelmode','manual')
			if xl(2)>yl(2)
				set(gca,'xtick',get(gca,'ytick'))
			else
				set(gca,'ytick',get(gca,'xtick'))
			end
			set(gca,'xticklabel',fixTickLabels(get(gca,'xtick')),'yticklabel',fixTickLabels(get(gca,'ytick')));
		end
		axis square
		xlabel([caption,' (First scan)'],'fontsize',fs)
		ylabel([caption,' (Second scan)'],'fontsize',fs)
		set(gca,'fontsize',fs)
		if makeFigure>0
			title('{\bfScan-Rescan Reproducibility}','fontsize',fs)
		end
		xc = mean(xlim);
		yl = ylim;
		yc = yl(2)-0.06*(yl(2)-yl(1));
		text(xc,yc,sprintf('Linear correlation %.2f, Mean difference %.1f%%',R,meanRelDif),'horizontalalignment','center','fontsize',fs)
		set(gca,'box','on')
	end

