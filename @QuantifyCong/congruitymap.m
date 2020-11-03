function congruitymap(qc, cong, Tsubs, Fsubs, Tibiofem, smooth)

cong = (1./cong)/10;

lent = length(Tsubs);
lenf = length(Fsubs);
data = Tibiofem;
Tibiofem(Tibiofem == 2) = 0;
Tibiofem(Tibiofem == 1) = 0;

if lent < lenf
    
    for i = 1:lent
        Tibiofem(Tsubs(i,1), Tsubs(i,2), Tsubs(i,3)) = cong(i);
    end
    Tibiofem = double(Tibiofem);
else
    for i = 1:lenf
        Tibiofem(Fsubs(i,1), Fsubs(i,2), Fsubs(i,3)) = cong(i);
    end
    Tibiofem = double(Tibiofem);
end
    
    figure
    if smooth % Not Good
        Tibiofem = smooth3(Tibiofem, 'gaussian',3,1);
    end
    
    % Renders Congruity over Contact AREA
    cdata = smooth3(Tibiofem, 'gaussian',3,1);
    fv = isosurface(Tibiofem, min(Tibiofem(:)));
    fv = smoothpatch(fv, 1,1,1,1); % Written by D.Kroon, see the function for details on how smoothing is done. Also OK without this.
    p = patch(fv);
    isonormals(Tibiofem,p);
    isocolors(cdata,p);
    caxis([0 5]);
    set(p, 'FaceColor', 'interp', 'EdgeColor', 'interp');
    hold on
    
    F = smooth3(data == 1,'gaussian',3,1); % Tibial
    p = patch(isosurface(F,0.1));
    set(p,'FaceColor','cyan','EdgeColor','none','AmbientStrength',.3);
    alpha(p, 0.1)
    hold on
    F = smooth3(data == 2,'gaussian',3,1); % Femoral 
    p1 = patch(isosurface(F,0.1));
    set(p1,'FaceColor',[1 0.6 0.7],'EdgeColor','none','AmbientStrength',.3);
    alpha(p1,0.1);
    l = light;
    light('Position', -get(l,'Position'))
    lighting gouraud
    axis off tight
    view(178,6)
        
% Create colorbar
colorbar([0.918802083333333 0.086708314335721 0.0208333333333333 0.815],...
    'FontSize',25);

% Create textbox
annotation('textbox',...
    [0.0400311192468617 0.31161891943857 0.105694560669456 0.0584415584415584],...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation('textbox',...
    [0.645187369246864 0.0332070755623919 0.105694560669456 0.0584415584415584],...
    'String',{'Posterior'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textarrow
annotation('textarrow',[0.559623430962342 0.553347280334727],...
    [0.459415584415584 0.327922077922078],'TextEdgeColor','none','FontSize',30,...
    'String',{'Contact Area'});

% Create textbox
annotation('textbox',...
    [0.0103436192468618 0.347958084983657 0.242781380753138 0.0584415584415584],...
    'String',{'Center of knee'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create arrow
annotation('arrow',[0.0453125 0.16171875],...
    [0.210305518169583 0.493943472409152]);

% Create doublearrow
annotation('doublearrow',[0.259375 0.7671875],...
    [0.0918667563930013 0.0915208613728129]);

% Create textbox
annotation('textbox',...
    [0.443468619246862 0.892857142857143 0.105694560669456 0.0584415584415584],...
    'String',{'KL 0'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation('textbox',...
    [0.256437369246862 0.0316727552393772 0.105694560669456 0.0584415584415584],...
    'String',{'Anterior'},...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textarrow
annotation('textarrow',[0.33046875 0.38297071129707],...
    [0.183465458663647 0.22025974025974],'TextEdgeColor','none','FontSize',30,...
    'String',{'Tibial cartilage'});

% Create textarrow
annotation('textarrow',[0.365625 0.306485355648536],...
    [0.65118912797282 0.508116883116883],'TextEdgeColor','none','FontSize',30,...
    'String',{'Femoral cartilage'});


    