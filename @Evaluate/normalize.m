function [val, has] = normalize(ev, val, norma, power)
	has = find(~isnan(norma));
   if length(find(has)) / length(val) < 0.9
      error('Too few normalization values?')
   end
	meanNorma = mean(norma(has));
	norma = (norma ./ meanNorma) .^ power;
   for c = 1:size(val,2)
      val(has,c) = val(has,c) ./ norma(has);
   end
