function DoDemographics(ev, quanKnees, Quan, Knees, subset)
	
	len = length(Knees);
	if exist('subset','var')
		insubset = zeros(len,1);
		for k = 1:len
			kn = Knees{k};
			if ~isempty(strmatch(kn,subset,'exact'));
				insubset(k) = 1;
			end
		end
		Knees = Knees(find(insubset));
		len = length(Knees);
	end
	
	age = -ones(len,1);
	bmi = age;
	sex = age;
	for k = 1:len
		% Age
		knee = Knees{k};
		age(k) = approxAge(knee);
		q = strmatch(knee,quanKnees);
		if ~isempty(q)
			V1 = Quan(q).V1;
			% BMI
			if isfield(V1,'Height') && isfield(V1,'Weight')
				bmi(k) = V1.Weight / (V1.Height/100)^2;
			end
			% Sex
			if isfield(V1,'Sex')
				if V1.Sex == 'M'
					sex(k) = 0;
				elseif V1.Sex == 'F'
					sex(k) = 1;
				end
			end
		end
	end
	disp(sprintf('Knee count: %d',len))
	doStats(age,'Age')
	doStats(bmi,'BMI')
	doStats(sex,'Sex')

function doStats(val,tag)
	idx = find(val>=0);
	val = val(idx);
	disp(sprintf('%7s (%d): %.1f - %.1f, mean %.2f, std %.1f',tag,length(val),min(val),max(val),mean(val),std(val)))
	
function age = approxAge(scan)
  % Approximating baseline age - assuming that they were all scanned on September 1, 2004
  birth = scan(1:6);
	year = str2num(birth(5:6));
	if year < 10
		year = 2000 + year;
	else
		year = 1900 + year;
	end
	n = datenum(year,str2num(birth(3:4)),str2num(birth(1:2)));
	scanday = datenum(2004,09,01); % assuming they were all scanned that day
	days = scanday-n;
	age = days / 365;
