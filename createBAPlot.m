function createBAPlot(X1, Y1)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 09-Dec-2011 10:56:07

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[197.695852534562 1197.69585253456]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-272.374532833762 227.625467166238]);
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot1 = plot((X1+Y1)/2,X1-Y1,'Marker','o','LineStyle','none','DisplayName','data 1');

% Create xlabel
xlabel({'Average of Scan-Rescan CA values'},'FontWeight','bold',...
    'FontSize',14);

% Create ylabel
ylabel({'Difference of Scan-Rescan CA values'},'FontWeight','bold',...
    'FontSize',14);

% Get xdata from plot
xdata1 = get(plot1, 'xdata');
% Get ydata from plot
ydata1 = get(plot1, 'ydata');
% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);

% Get axes xlim
axXLim1 = get(axes1, 'xlim');

% Find the std
ystd1 = std(ydata1);

% Prepare values to plot std; first find the mean
ymean1 = mean(ydata1);
% Compute bounds as mean +/- std
lowerBound1 = ymean1 - 1.96 * ystd1;
upperBound1 = ymean1 + 1.96 * ystd1;
% Get coordinates for the std bounds line
stdValue1 = [lowerBound1 lowerBound1 NaN upperBound1 upperBound1 NaN];
axXStdLim1 = [axXLim1 NaN axXLim1 NaN];

% Plot the bounds
statLine1 = plot(axXStdLim1,stdValue1,'DisplayName','   y std',...
    'Parent',axes1,...
    'Tag','std y',...
    'LineStyle','-.',...
    'Color',[0.75 0 0.75]);

% Set new line in proper position
setLineOrder(axes1, statLine1, plot1);

% Find the mean
ymean2 = mean(ydata1);
% Get coordinates for the mean line
meanValue1 = [ymean2 ymean2];
% Plot the mean
statLine2 = plot(axXLim1,meanValue1,'DisplayName','   y mean',...
    'Parent',axes1,...
    'Tag','mean y',...
    'LineStyle','-.',...
    'Color',[0 0.5 0]);

% Set new line in proper position
setLineOrder(axes1, statLine2, plot1);

% Get xdata from plot
xdata1 = get(plot1, 'xdata');
% Get ydata from plot
ydata1 = get(plot1, 'ydata');
% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);

% Get axes xlim
axXLim1 = get(axes1, 'xlim');

% Find the std
ystd1 = std(ydata1);

% Prepare values to plot std; first find the mean
ymean1 = mean(ydata1);
% Compute bounds as mean +/- std
lowerBound1 = ymean1 - ystd1;
upperBound1 = ymean1 + ystd1;
% Get coordinates for the std bounds line
stdValue2 = [lowerBound1 lowerBound1 NaN upperBound1 upperBound1 NaN];
axXStdLim2 = [axXLim1 NaN axXLim1 NaN];

% Plot the bounds
statLine3 = plot(axXStdLim2,stdValue2,'DisplayName','   y std',...
    'Parent',axes1,...
    'Tag','std y',...
    'LineStyle','-.',...
    'Color',[0.75 0 0.75]);

% Set new line in proper position
setLineOrder(axes1, statLine3, plot1);

%-------------------------------------------------------------------------%
function setLineOrder(axesh1, newLine1, associatedLine1)
%SETLINEORDER(AXESH1,NEWLINE1,ASSOCIATEDLINE1)
%  Set line order
%  AXESH1:  axes
%  NEWLINE1:  new line
%  ASSOCIATEDLINE1:  associated line

% Get the axes children
hChildren = get(axesh1,'Children');
% Remove the new line
hChildren(hChildren==newLine1) = [];
% Get the index to the associatedLine
lineIndex = find(hChildren==associatedLine1);
% Reorder lines so the new line appears with associated data
hNewChildren = [hChildren(1:lineIndex-1);newLine1;hChildren(lineIndex:end)];
% Set the children:
set(axesh1,'Children',hNewChildren);
