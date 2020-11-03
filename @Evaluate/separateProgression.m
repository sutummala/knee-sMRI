function separateProgression(ev, KL1, KL2, Vals, splitNonProg, yTag, doLog)

	if length(find(KL1<0))
		error('I did not expect this!')
	end
	if splitNonProg
		separateProgression3(KL1, KL2, Vals, yTag, doLog)
	else
		separateProgression2(KL1, KL2, Vals, yTag, doLog)
	end
	
function separateProgression3(KL04, KL06, CTX, yTag, doLog)
	ok = find(KL06>=0);
	KL04 = KL04(ok); KL06 = KL06(ok); CTX = CTX(ok);
	health = find(KL04==0 & KL06==KL04);
	oastat = find(KL04> 0 & KL06==KL04);
	oaprog = find(KL06>KL04);
	if ~isempty(find(KL06<KL04))
		error('This is weird')
	end
	if doLog
		CTX = log(CTX);
	end
	h = CTX(health);
	s = CTX(oastat);
	p = CTX(oaprog);
	means = [mean(h) mean(s) mean(p)];
	sems = [sem(h) sem(s) sem(p)];
	if doLog
		sems = exp(means+sems)-exp(means);
		means = exp(means);
	end
	figure
	bar(means)
	hold on
	eb = errorbar(means, sems);
	set(eb,'linewidth',2,'markersize',16,'LineStyle','none','color','red','marker','.');
	set(gca,'xtick',1:3,'xticklabel',{'Healthy','OA, Static','OA, Progressive'})
	[hy, Ps, ci] = ttest2(h, s);
	[hy, Pp, ci] = ttest2(h, p);
	yl = ylim;
	d = 0.05*(yl(2)-yl(1));
	c1 = means(2)+sems(2)+d;
	c2 = means(3)+sems(3)+d;
	line([1 2],[c1 c1]), text(1.5, c1+d, sprintf('p = %.3f',Ps),'horizontalalignment','center')
	line([1 3],[c2 c2]), text(2, c2+d, sprintf('p = %.3f',Pp),'horizontalalignment','center')
	ylim([yl(1), c2+4*d])
	ylabel(yTag)
	y = min(ylim)-0.08*(max(ylim)-min(ylim));
	text(1,y,sprintf('N=%d',length(h)),'horizontalalignment','center')
	text(2,y,sprintf('N=%d',length(s)),'horizontalalignment','center')
	text(3,y,sprintf('N=%d',length(p)),'horizontalalignment','center')

function separateProgression2(KL04, KL06, CTX, yTag, doLog)
	ok = find(KL06>=0);
	KL04 = KL04(ok); KL06 = KL06(ok); CTX = CTX(ok);
	stat = find(KL06==KL04);
	prog = find(KL06>KL04);
	if ~isempty(find(KL06<KL04))
		error('This is weird')
	end
	if doLog
		CTX = log(CTX);
	end
	s = CTX(stat);
	p = CTX(prog);
	means = [mean(s) mean(p)];
	sems = [sem(s) sem(p)];
	if doLog
		sems = exp(means+sems)-exp(means);
		means = exp(means);
	end
	figure
	bar(means)
	xlim([0.5 2.5])
	hold on
	eb = errorbar(means, sems);
	set(eb,'linewidth',2,'markersize',16,'LineStyle','none','color','red','marker','.');
	set(gca,'xtick',1:2,'xticklabel',{'Non-Progressive','Progressive OA'})
	[hy, Pp, ci] = ttest2(s, p);
	yl = ylim;
	d = 0.05*(yl(2)-yl(1));
	c1 = double(max(means(2)+sems(2),means(1)+sems(1))+d); % text chokes on single 
	line([1 2],[c1 c1]) 
	text(1.5, c1+d, sprintf('p = %.3f',Pp),'horizontalalignment','center')
	ylim([yl(1), c1+4*d])
	ylabel(yTag)
	y = min(ylim)-0.08*(max(ylim)-min(ylim));
	text(1,y,sprintf('N=%d',length(s)),'horizontalalignment','center')
	text(2,y,sprintf('N=%d',length(p)),'horizontalalignment','center')

function s = sem(v)
	s = std(v) / sqrt(length(v));
	
