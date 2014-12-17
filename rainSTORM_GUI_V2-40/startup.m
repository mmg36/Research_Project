function startup(gui)
    addpath(genpath(pwd));
    
    if nargin > 1
        error('Program usage: startup (gui default) or startup(''gui'') or startup(''nogui'')');
    elseif nargin == 0
        gui='gui';
    end
    
    params = rainSTORM_params_struct();
    switch gui;
        case 'gui'
            rainSTORM(params);
        case 'nogui'
            rainSTORM_NOGUI();
        otherwise
            error('Possible parameters: gui or nogui');
    end
end