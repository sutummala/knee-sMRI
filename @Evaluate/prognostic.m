function [p,N,OR,auc] = prognostic(ev,val,bKL,fKL,caption,OAKL)




if ~exist('OAKL')
    OAKL = 1;
end
if ~exist('makeFigure')
    makeFigure = 0;
end

l = length(val);
P = [];
NP = [];
base = val(:,1);
bKL0 = []; bKL1 = []; bKL2 = []; bKL3 = []; bKL4 = [];

fKL0 = []; fKL1 = []; fKL2 = []; fKL3 = []; fKL4 = [];

for i = 1:l
    bas = base(i);
    if isnan(bKL(i)) || isnan(fKL(i)) || isnan(bas);
       continue
    end
    if bKL(i)<0 || fKL(i)<0 % Check if we have KL value for the knee
      disp('No KL!')
      continue
    end
   switch bKL(i)
      case 0, bKL0(end+1) = bas;
      case 1, bKL1(end+1) = bas;
      case 2, bKL2(end+1) = bas;
      case 3, bKL3(end+1) = bas;
      case 4, bKL4(end+1) = bas;
   end
%    switch fKL(i)
%       case 0, fKL0(end+1) = fol;
%       case 1, fKL1(end+1) = fol;
%       case 2, fKL2(end+1) = fol;
%       case 3, fKL3(end+1) = fol;
%       case 4, fKL4(end+1) = fol;
%    end
%    
   if bKL(i) < OAKL && fKL(i) < OAKL
       NP(end+1) = bas;
   else if bKL(i) < OAKL && fKL(i) >= OAKL
           P(end+1) = bas;
       end
   end
end     

[h,p] = ttest2(NP,P);
N = SampleSize(NP,P);
[OR, stars] = RatioVsLow(NP,P,0.50);
[dA, sddA, pauc] = DeLongTest(NP,P,1);
auc = AUC([NP, P],[zeros(1,length(NP)),ones(1,length(P))], 0 * makeFigure);

disp(sprintf('%20s: %3d NonProgressors (%7.6f) and %3d Progressors   (%7.6f)(%d,%d,%d,%d)-> p %.20f (ss %4.0f or %3.1f, %3.1f auc %3.3f, %1.5f)',...
   caption,length(NP),mean(NP),length(P),mean(P),(length(fKL1)-length(bKL1)),(length(fKL2)-length(bKL2)),(length(fKL3)-length(bKL3)),...
   (length(fKL4)-length(bKL4)),p,SampleSize(NP,P),OR(2), stars(2),...
   AUC([NP, P],[zeros(1,length(NP)),ones(1,length(P))], 0 * makeFigure),pauc));

if makeFigure
   figure
   fs = 20;
   ebKL0 = std(bKL0);
   ebKL1 = std(bKL1);
   ebKL2 = std(bKL2);
   ebKL3 = std(bKL3);
   ebKL4 = std(bKL4);
   efKL0 = std(fKL0);
   efKL1 = std(fKL1);
   efKL2 = std(fKL2);
   efKL3 = std(fKL3);
   efKL4 = std(fKL4);
   
      ebKL0 = ebKL0 / sqrt(length(bKL0));
      ebKL1 = ebKL1 / sqrt(length(bKL1));
      ebKL2 = ebKL2 / sqrt(length(bKL2));
      ebKL3 = ebKL3 / sqrt(length(bKL3));
      ebKL4 = ebKL4 / sqrt(length(bKL4));
      efKL0 = efKL0 / sqrt(length(fKL0));
      efKL1 = efKL1 / sqrt(length(fKL1));
      efKL2 = efKL2 / sqrt(length(fKL2));
      efKL3 = efKL3 / sqrt(length(fKL3));
      efKL4 = efKL4 / sqrt(length(fKL4));
    
      
      hold on
      hb = bar(1:4,([(mean(bKL0)-mean(fKL0))/mean(bKL0),(mean(bKL1)-mean(fKL1))/mean(bKL1), (mean(bKL2)-mean(fKL2))/mean(bKL2), ...
          ((mean(bKL3)+mean(bKL4))-(mean(bKL3)+mean(fKL4)))/(mean(bKL3)+mean(bKL4))])*100,'w');
      he = errorbar(1:4,([(mean(bKL0)-mean(fKL0))/mean(bKL0),(mean(bKL1)-mean(fKL1))/mean(bKL1), ...
          (mean(bKL2)-mean(fKL2))/mean(bKL2), ((mean(bKL3)+mean(bKL4))-(mean(bKL3)+mean(fKL4)))/(mean(bKL3)+mean(bKL4))])*100, ...
          ([(ebKL0-efKL0)/ebKL0,(ebKL1-efKL1)/ebKL1, (ebKL2-efKL2)/ebKL2, ((ebKL3+ebKL4)-(efKL3+efKL4))/(ebKL3+ebKL4)])*100);
      ymin = min(([(mean(bKL0)-mean(fKL0))/mean(bKL0),(mean(bKL1)-mean(fKL1))/mean(bKL1), (mean(bKL2)-mean(fKL2))/mean(bKL2), ...
          ((mean(bKL3)+mean(bKL4))-(mean(bKL3)+mean(fKL4)))/(mean(bKL3)+mean(bKL4))])*100 - ([(ebKL0-efKL0)/ebKL0,(ebKL1-efKL1)/ebKL1,...
          (ebKL2-efKL2)/ebKL2, ((ebKL3+ebKL4)-(efKL3+efKL4))/(ebKL3+ebKL4)])*100);
      ymax = max(([(mean(bKL0)-mean(fKL0))/mean(bKL0),(mean(bKL1)-mean(fKL1))/mean(bKL1), (mean(bKL2)-mean(fKL2))/mean(bKL2), ...
          ((mean(bKL3)+mean(bKL4))-(mean(bKL3)+mean(fKL4)))/(mean(bKL3)+mean(bKL4))])*100 + ([(ebKL0-efKL0)/ebKL0,(ebKL1-efKL1)/ebKL1,...
          (ebKL2-efKL2)/ebKL2, ((ebKL3+ebKL4)-(efKL3+efKL4))/(ebKL3+ebKL4)])*100);
      
      set(hb,'FaceColor',[1,1,0])
      set(hb,'FaceColor',[.5,.5,.5])
      set(he,'linewidth',2,'markersize',16,'LineStyle','none','color','black','marker','.');
      xlabel('Kellgren & Lawrence Index','fontsize',fs)
      ylabel({caption,'% Fractional Longitudinal Change'},'fontsize',fs)
end

% if OAKL == 1
%          set(gca,'xtick',1:4,'xticklabel',{'0','1','2','3/4'},'fontsize',fs)
%       else         
%          set(gca,'xtick',1:4,'xticklabel',{['0-',num2str(OAKL-1)],['>',num2str(OAKL-1)],'0','1','2','3/4'},'fontsize',fs)
% end
% 
% ysz = ymax-ymin;
% %  if 1
% %      [H,P] = ttest(KL0);
% %      makestars(P, mean(KL0), eKL0, mean(KL0), eKL0, 5, 5, ysz)
% %        [H,P] = ttest(KL1);
% %             makestars(P, mean(KL1), eKL1, mean(KL1), eKL1, 6, 6, ysz)
% %             [H,P] = ttest(KL2);
% %             makestars(P, mean(KL2), eKL2, mean(KL2), eKL2, 7, 7, ysz)
% %             [H,P] = ttest(KL3);
% %             makestars(P, mean(KL3), eKL3, mean(KL3), eKL3, 8, 8, ysz)
% %  else
% starpair(bKL0, fKL0, ebKL0, efKL0, 1, 2, ysz)
% starpair(bKL1, fKL1, ebKL1, efKL1, 2, 3, ysz)
% starpair(bKL2, fKL2, ebKL2, efKL2, 3, 4, ysz)
% starpair(bKL3, fKL3, ebKL3, efKL3, 4, 4, ysz)
% %  end
% function starpair(KL0, KL1, eKL0, eKL1, i1, i2, ysz)
%   if numel(KL0) && numel(KL1)
%      [Hx, Px, CIx] = ttest2(KL0, KL1);
%      makestars(Px, mean(KL0), eKL0, mean(KL1), eKL1, i1, i2, ysz)
%   end
% 
% function makestars(P, a, ea, b, eb, from, to, ysize)
% if P < 0.05
%    star = '*';
%    if P < 0.01, star = '**'; end
%    if P < 0.001, star = '***'; end
%    if P < 0.0001, star = '****'; end
%    if max(a,b)>0
%       y = max([a+ea, b+eb]) + 0.03*ysize;
%       sy = y + 0.03*ysize;
%    else
%       y = min([a-ea, b-eb]) - 0.03*ysize;
%       sy = y - 0.08*ysize;
%    end
%    if from~=to
%       line([from+.1,to-.1],[y,y],'color','k','linewidth',4)
%    end
%    h = text((from+to)/2,sy,star);
%    set(h,'horizontalalignment','center','verticalalignment','middle','FontSize',16)
% end


                