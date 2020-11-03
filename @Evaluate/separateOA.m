function [Healthy,OA] = separateOA(ev, vals, KL, caption, makeFigure, OAKL, doKL01, addBarStar)

doSTD = false; % Alternative is to do SEM

if ~exist('makeFigure')
   makeFigure = 1;
end
if ~exist('OAKL')
   OAKL = 1;
end
if ~exist('doKL01')
   doKL01 = 1;
end
doKLhist = 1; % Choose 2 to visualise only cross-sectional seperation between different KL.
if ~exist('addBarStar')
   addBarStar = 1;
end

Healthy = [];
OA = [];
% OAKLhist = zeros(4,1);
count = length(vals);
KL0 = []; KL1 = []; KL2 = []; KL3 = []; KL4 = [];
for p = 1:count
   val = vals(p);
   if isnan(val) || isnan(KL(p))
      continue
   end
   if KL(p)<0 % Check if we have KL value for the knee
      disp('No KL!')
      continue
   end
   switch KL(p)
      case 0, KL0(end+1) = val;
      case 1, KL1(end+1) = val;
      case 2, KL2(end+1) = val;
      case 3, KL3(end+1) = val;
      case 4, KL4(end+1) = val;
   end
   if KL(p) >= OAKL
        OA(end+1) = val;
      % OAKLhist(KL(p)) = OAKLhist(KL(p))+1;
   else 
        Healthy(end+1) = val;
   end
end
if length(KL2)+length(KL3)+length(KL4)==0
   doKLhist = 0;
   doKL01 = 0;
end

% Do statistics on Helthy and OA & by KL scores
P = seperateAll(Healthy, OA, KL1, KL2, KL3, KL4, caption, makeFigure);
seperateKL(KL0, KL1, makeFigure, 0, 1);
seperateKL(KL1, KL2, makeFigure, 1, 2);
seperateKL(KL2, KL3, makeFigure, 2, 3);

if makeFigure
   figure
   %subplot(2,2,4)
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
   hold on
   if ~doKLhist
      hb = bar(1:2,[mean(Healthy),mean(OA)],'w');
      he = errorbar(1:2,[mean(Healthy),mean(OA)], [eHe,eOA]);
      ymin = min([mean(Healthy),mean(OA)] - [eHe,eOA]);
      ymax = max([mean(Healthy),mean(OA)] + [eHe,eOA]);
   else if doKLhist == 2
      hb = bar(1:4,[mean(KL0),mean(KL1),mean(KL2),mean(KL3)],'w');
      he = errorbar(1:4,[mean(KL0),mean(KL1),mean(KL2),mean(KL3)], ...
         [eKL0,eKL1,eKL2,eKL3]);
      ymin = min([mean(KL0),mean(KL1),mean(KL2),mean(KL3)]-[eKL0,eKL1,eKL2,eKL3]);
      ymax = max([mean(KL0),mean(KL1),mean(KL2),mean(KL3)]+[eKL0,eKL1,eKL2,eKL3]);
       else
      hb = bar([1:2,5:8],[mean(Healthy),mean(OA),mean(KL0),mean(KL1),mean(KL2),mean(KL3)],'w');
      he = errorbar([1:2,5:8],[mean(Healthy),mean(OA),mean(KL0),mean(KL1),mean(KL2),mean(KL3)], ...
         [eHe,eOA,eKL0,eKL1,eKL2,eKL3]);
      ymin = min([mean(Healthy),mean(OA),mean(KL0),mean(KL1),mean(KL2),mean(KL3)]-[eHe,eOA,eKL0,eKL1,eKL2,eKL3]);
      ymax = max([mean(Healthy),mean(OA),mean(KL0),mean(KL1),mean(KL2),mean(KL3)]+[eHe,eOA,eKL0,eKL1,eKL2,eKL3]);
%         ymin = 400;
%         ymax = 900;
       end
   end
   set(hb,'FaceColor',[1,1,0])
   set(hb,'FaceColor',[.5,.5,.5])
   set(he,'linewidth',2,'markersize',16,'LineStyle','none','color','black','marker','.');
   if addBarStar
      ysz = ymax-ymin;
      if addBarStar < 0
         [H,P] = ttest(Healthy);
         makestars(P, mean(Healthy), eHe, mean(Healthy), eHe, 1, 1, ysz)
         [H,P] = ttest(OA);
         makestars(P, mean(OA), eOA, mean(OA), eOA, 2, 2, ysz)
         if doKLhist
            [H,P] = ttest(KL0);
            makestars(P, mean(KL0), eKL0, mean(KL0), eKL0, 5, 5, ysz)
            [H,P] = ttest(KL1);
            makestars(P, mean(KL1), eKL1, mean(KL1), eKL1, 6, 6, ysz)
            [H,P] = ttest(KL2);
            makestars(P, mean(KL2), eKL2, mean(KL2), eKL2, 7, 7, ysz)
            [H,P] = ttest(KL3);
            makestars(P, mean(KL3), eKL3, mean(KL3), eKL3, 8, 8, ysz)
         end
         ymin = ymin - 0.1*(ymax-ymin);
      else
          if doKLhist == 1
         makestars(P, mean(Healthy), eHe, mean(OA), eOA, 1, 2, ysz)
          end
         if doKLhist == 2
            starpair(KL0, KL1, eKL0, eKL1, 1, 2, ysz)
            starpair(KL1, KL2, eKL1, eKL2, 2, 3, ysz)
            starpair(KL2, KL3, eKL2, eKL3, 3, 4, ysz)
         else 
            starpair(KL0, KL1, eKL0, eKL1, 5, 6, ysz)
            starpair(KL1, KL2, eKL1, eKL2, 6, 7, ysz)
            starpair(KL2, KL3, eKL2, eKL3, 7, 8, ysz)
         end
%             [Hx, Px, CIx] = ttest2(KL0, KL1);            
%             makestars(Px, mean(KL0), eKL0, mean(KL1), eKL1, 5, 6, ysz)
%             [Hx, Px, CIx] = ttest2(KL1, KL2);
%             makestars(Px, mean(KL1), eKL1, mean(KL2), eKL2, 6, 7, ysz)
%             [Hx, Px, CIx] = ttest2(KL2, KL3);
%             makestars(Px, mean(KL2), eKL2, mean(KL3), eKL3, 7, 8, ysz)
       end
         ymax = ymax + 0.1*(ymax-ymin);
    end
 end
   ymin = ymin - 0.2*ysz;
   ymax = ymax + 0.2*ysz;
   ylim([ymin,ymax])
   if doKLhist
      if OAKL == 1
          if doKLhist == 2
         set(gca,'xtick',1:4,'xticklabel',{'0','1','2','3/4'},'fontsize',fs)
          else
         set(gca,'xtick',[1:2,5:8],'xticklabel',{'0','>0','0','1','2','3/4'},'fontsize',fs)
          end
      else         
         set(gca,'xtick',[1:2,5:8],'xticklabel',{['0-',num2str(OAKL-1)],['>',num2str(OAKL-1)],'0','1','2','3/4'},'fontsize',fs)
      end
      if 1
         xlabel('Kellgren & Lawrence Index')
      else
         xl = xlabel(''); % find place for extra text
         pos = get(xl,'Position');
         y = pos(2);
         text(1,y,sprintf('N=%d',length(Healthy)),'horizontalalignment','center','fontsize',fs)
         text(2,y,sprintf('N=%d',length(OA)),'horizontalalignment','center','fontsize',fs)
         text(5,y,sprintf('N=%d',length(KL0)),'horizontalalignment','center','fontsize',fs)
         text(6,y,sprintf('N=%d',length(KL1)),'horizontalalignment','center','fontsize',fs)
         text(7,y,sprintf('N=%d',length(KL2)),'horizontalalignment','center','fontsize',fs)
         text(8,y,sprintf('N=%d',length(KL3)),'horizontalalignment','center','fontsize',fs)
         text(9,y,sprintf('N=%d',length(KL4)),'horizontalalignment','center','fontsize',fs)
      end
      if doKLhist == 1
      line([3.5 3.5],[ymin,ymax],'color','black','linestyle',':')
      end
   else
      set(gca,'xtick',1:2,'xticklabel',{'Healthy','OA'},'fontsize',fs)
   end
   if doSTD
      ylabel({caption,'Mean and Std. Dev.'},'fontsize',fs)
   else
      ylabel({caption,'Mean and SEM'},'fontsize',fs)
   end
   set(gca,'xtickmode','manual','ytickmode','manual','xticklabelmode','manual','yticklabelmode','manual')
   % title('{\bfSeparation of OA from Healthy}','fontsize',fs)
   if ~doKLhist
      xlim([0.5,2.5])
   else if doKLhist == 2
            xlim([0.25,4.75])
       else
            xlim([0.25,8.75]) 
       end
   end
   set(gca,'yticklabel',fixTickLabels(get(gca,'ytick')));
   if 0 && ~doKLhist
      xc = mean(xlim);
      yl = ylim;
      yc = yl(2)-0.06*(yl(2)-yl(1));
      if P < 0.05
         txt = sprintf('Successful t-test (p = %.5f)',P);
      else
         txt = sprintf('Failed t-test (p = %.5f)',P);
      end
      text(xc,yc,txt,'horizontalalignment','center','fontsize',fs)
   end
   axis square
end

function starpair(KL0, KL1, eKL0, eKL1, i1, i2, ysz)
  if numel(KL0) && numel(KL1)
     [Hx, Px, CIx] = ttest2(KL0, KL1);
     makestars(Px, mean(KL0), eKL0, mean(KL1), eKL1, i1, i2, ysz)
  end
end

function makestars(P, a, ea, b, eb, from, to, ysize)
if P < 0.05
   star = '*';
   if P < 0.01, star = '**'; end
   if P < 0.001, star = '***'; end
   if P < 0.0001, star = '****'; end
   if max(a,b)>0
      y = max([a+ea, b+eb]) + 0.03*ysize;
      sy = y + 0.03*ysize;
   else
      y = min([a-ea, b-eb]) - 0.03*ysize;
      sy = y - 0.08*ysize;
   end
   if from~=to
      line([from+.1,to-.1],[y,y],'color','k','linewidth',4)
   end
   h = text((from+to)/2,sy,star);
   set(h,'horizontalalignment','center','verticalalignment','middle','FontSize',16)
end
end

function [P] = seperateAll(Healthy, OA, KL1, KL2, KL3, KL4, caption, makefig)
[H, P, CI] = ttest2(Healthy, OA);
[P1] = ranksum(Healthy, OA);
[N] = SampleSize(Healthy,OA);
[OR,stars] = RatioVsLow(Healthy,OA,0.50);
[dA, sddA, p] = DeLongTest(Healthy,OA,1);
SRM = (mean(OA)-mean(Healthy))/std(horzcat(Healthy,OA));
[x, y, t, auc] = perfcurve([zeros(1,length(Healthy)),ones(1,length(OA))],[Healthy, OA], 0, 'NBoot', 100);
disp(sprintf('%20s: %3d Healthy (%7.6f) and %3d OA (%7.6f)(%d,%d,%d,%d)-> p Ttest%.25f/Ranksum%.25f (ss %4.0f or %3.1f,%3.1f auc %3.3f,CI(%3.3f-%3.3f),%0.15f, srm %3.3f)',...
   caption,length(Healthy),mean(Healthy),length(OA),mean(OA),length(KL1),length(KL2),length(KL3),length(KL4),P, P1,...
   N,OR(2),stars(2),...
   AUC([Healthy, OA],[zeros(1,length(Healthy)),ones(1,length(OA))], 0 * makefig),auc(2),auc(3),p,SRM));
end


function seperateKL(KL0, KL1, makefig, fir, las)
    [H, P01, CI] = ttest2(KL0, KL1);
    [P11] = ranksum(KL0, KL1);
    [N1] = SampleSize(KL0, KL1);
    [OR1,stars1] = RatioVsLow(KL0, KL1,0.50);
    [dA, sddA, p1] = DeLongTest(KL0, KL1,1);
    SRM1 = (mean(KL1)-mean(KL0))/std(horzcat(KL0, KL1));
    [x, y, t, auc1] = perfcurve([zeros(1,length(KL0)),ones(1,length(KL1))],[KL0, KL1], 0, 'NBoot', 100);
   disp(sprintf('%20s  %3d KL%d (%9.4f) and %3d KL%d (%9.4f) -> p Ttest%.25f/Ranksum%.25f (ss %4.0f or %3.1f,%3.1f auc %3.3f,CI(%3.3f-%3.3f),%0.15f, srm %3.3f)',...
      ' ',length(KL0),fir,mean(KL0),length(KL1),las,mean(KL1),P01,P11,N1,OR1(2),stars1(2),...
      AUC([KL0, KL1],[zeros(1,length(KL0)),ones(1,length(KL1))], 0 * makefig),auc1(2),auc1(3),p1,SRM1));
end
