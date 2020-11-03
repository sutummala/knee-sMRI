function [congruity] = correctSomething(Val, KL, Age, BMI, sex)

hel = KL(:,1) == 0; % healthy

allmean = nanmean(Val(hel,1));

% Initialization of Residuals 
residuals = Val;
residuals(:) = NaN;

men = sex(:,1) == 0; % Male 
wmen = sex(:,1) == 1; % Female

%% Correct for Age
[Polm, r ,pm] = fitdata(Val, Age, hel, men);
[Polw, r1 ,pw] = fitdata(Val, Age, hel, wmen);

if pm(2) < 0.05
    residuals(men,1) = Val(men,1) - (Polm(2) + Polm(1).*Age(men,1));
    cong(men,1) = residuals(men,1) + allmean;
else
    cong(men,1) = Val(men,1);
end

if pw(2) < 0.05
    residuals(wmen,1) = Val(wmen,1) - (Polw(2) + Polw(1).*Age(wmen,1));
    cong(wmen,1) = residuals(wmen,1) + allmean;
else
    cong(wmen,1) = Val(wmen,1);
end

%% Correct for BMI
clear Polm Polw pm pw p
[Polm, r ,pm] = fitdata(cong, BMI, hel, men);
[Polw, r ,pw] = fitdata(cong, BMI, hel, wmen);

if pm(2) < 0.05
    residuals(men,1) = cong(men,1) - (Polm(2) + Polm(1).*BMI(men,1));
    congruity(men,1) = residuals(men,1) + allmean;
else
    congruity(men,1) = cong(men,1);
end

if pw(2) < 0.05
    residuals(wmen,1) = cong(wmen,1) - (Polw(2) + Polw(1).*BMI(wmen,1));
    congruity(wmen,1) = residuals(wmen,1) + allmean;
else
    congruity(wmen,1) = cong(wmen,1);
end
 

function [pol, r, p] = fitdata(S, A, hel, gen)

S = S(hel & gen, 1); % Only healthy and specified gender
A = A(hel & gen, 1);
ok = ~isnan(S) & ~isnan(A);
factor = A(ok);
sm = S(ok);
[r, p] = corrcoef(factor, sm);
pol = polyfit(factor,sm,1);

