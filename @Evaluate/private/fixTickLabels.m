function labels = fixTickLabels(ticks)
	labels = cell(length(ticks));
   dif = max(ticks)-min(ticks);
   decs = -round(log10(dif));
   decs = max(0,decs+1);
   formatstr = sprintf('%%.%.0ff',decs);
   for l = 1:length(ticks)
      labels{l} = sprintf(formatstr,ticks(l));
   end
