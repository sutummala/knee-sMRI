function figs

  go('C:\Documents and Settings\stu\Desktop\Congruity\Matlab\ContactR\SCA.fig',...
     'Estimated Sample Size','SCA.jpg')

%   go('C:\Documents and Settings\stu\Desktop\Congruity\Matlab\ContactR\CA.fig',...
%      'RMS CV [%]','CA.jpg')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function go(figfile,ylab,jpeg)
  uiopen(figfile,1)
  set(gca,'fontsize',12)
%   legend('Tibial','Femoral')
%   xlabel('Number of Iterations')
%   ylabel(ylab)
%   title('Medial Knee Compartments')
  % print(gcf,'-depsc','FigSS.eps') % looks crap in Word
  print(gcf,'-djpeg99',jpeg)
  close
