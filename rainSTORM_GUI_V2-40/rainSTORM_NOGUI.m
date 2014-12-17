
function varargout = rainSTORM_NOGUI(varargin)

if nargin == 1
    params=varargin{1};
else
    params = readparam('parameters.txt');
    [pathstr, name, ext] = fileparts(params.rawdata_mgr.filename);

    params.rawdata_mgr.filename=fullfile(pathstr,name);
    params.rawdata_mgr.ext=ext;
end

display(params.filename);


params = rainSTORM_main(params);

assignin('base','params',params);
end

function output = readparam(file)
    fid=fopen(file);
    
    f1 = textscan(fid, '%s = %q',1);  %filename
    f2 = textscan(fid, '%s = %s',1);  %algorithm
    f3 = textscan(fid, '%s = %s',1);  %pixelWidth
    f4 = textscan(fid, '%s = %s',1);  %initSig
    f5 = textscan(fid, '%s = %s',1);  %rad
    f6 = textscan(fid, '%s = %s',1);  %tol
    f7 = textscan(fid, '%s = %s',1);  %Thresh
    f8 = textscan(fid, '%s = %s',1);  %maxIts
    f9 = textscan(fid, '%s = %s',1);  %scalebar
    f10 = textscan(fid, '%s = %s',1); %sumimage
    
    
    output = rainSTORM_params_struct();
    
    output.rawdata_mgr.filename = f1{2};
    output.rawdata_mgr.myImInfo.pixelWidth= f3{2};
    
    output.localization.settings.Thresh = f7{2};
    output.localization.algo_id = f2{2};
    output.localization.settings.initSig = f4{2};
    output.localization.settings.maxIts = f8{2};
    output.localization.settings.rad = f5{2};
    output.localization.settings.tol = f6{2};
    
    output.flags.SB = f9{2};
    output.flags.Sum = f10{2};
    
    fclose(fid);
end

