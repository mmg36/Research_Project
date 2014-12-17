function st = rainSTORM_params_struct()

% % % myFits(lpPx,:)=[fitRowPos,fitColPos];

% % % myParams(lpPx,1)=myPixels(lpPx,3); % Averaged magnitude of this signal

% % % myParams(lpPx,2)=(residueRows+residueCols)/2; % Mean critical tol for fit

% % % myParams(lpPx,3)=sum(yCols); % Sum of signal (counts) for this fit

% % % myParams(lpPx,4)=sigX;  % X-width (sigma, rows, fitted) of this Gaussian

% % % myParams(lpPx,5)=sigY;  % Y-width (sigma, cols, fitted) of this Gaussian

% % % myParams(lpPx,6)=bkgdSig;  % Background for each ROI

% % %

% % % SupResParams(lpPx).frame_idx=frame_index;

% % % SupResParams(lpPx).x=fitRowPos;

% % % SupResParams(lpPx).y=fitColPos;

% % % SupResParams(lpPx).z=0;

% % % SupResParams(lpPx).I=myPixels(lpPx,3); % Averaged magnitude of this signal

% % % SupResParams(lpPx).sig_x=sigX;  % X-width (sigma, rows, fitted) of this Gaussian

% % % SupResParams(lpPx).sig_y=sigY;  % Y-width (sigma, cols, fitted) of this Gaussian

% % % SupResParams(lpPx).avg_brigthness=bkgdSig;  % Background for each ROI

% % % SupResParams(lpPx).res=(residueRows+residueCols)/2; % Mean critical tol for fit

% % % SupResParams(lpPx).res_Row=residueRows;

% % % SupResParams(lpPx).res_Col=residueCols;

% % % SupResParams(lpPx).Sum_signal=sum(yCols); % Sum of signal (counts) for this fit



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global parameters
global parCompTool imgProcTool

    v = ver;
    parCompTool = 0;
    imgProcTool = 0;
    for k=1:length(v)
        if strfind(v(k).Name, 'Parallel Computing Toolbox')
            parCompTool = str2num((sprintf( v(k).Version)));
        end
        if strfind(v(k).Name, 'Image Processing Toolbox')
            imgProcTool = str2num((sprintf( v(k).Version)));
        end
    end
    if imgProcTool <= 0
        error('There is no Image Processing Toolbox installed! ');
    end
    
    if (parCompTool > 0) && (matlabpool('size') == 0)
        matlabpool open;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rawdata_mgr

rawdata_mgr=struct(...
    'ext',[],...
    'filename',[],...
    'countsPerPhoton',[],...
    'myImInfo',[]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% localization



% SupResIm: delete, this shall be generated by reviewer, calls of rainSTORM_recon and imshow has to be

% removed from rainSTORM_main_GUI, actually rainSTORM_main_GUI shall be renamed as rainSTORM_main

% SupResNumber: -> numberOfLocs

% SupResParams: redefine as planed

% SupResPosits: delete, delegate into SupResParams

% Thresh: settings

% alg: -> algo_id (string)

% allowSig: settings

% allowX: settings

% estNum: settings

% ext: -> rawdata_mgr

% filename: -> rawdata_mgr

% initSig: settings

% initX0: settings

% maxIts: settings

% myFrame: delete, this just the last localized frame, in rainSTORM_reviewer.m replace size(myFrame) with

% params.rawdata_mgr.myImInfo.frame_size, see later at sizeOfCCDFrame

% numberOfFiles: -> numberOfFrames

% prevSF: settings

% rad: settings

% sizeOfCCDFrame: -> rawdata_mgr.myImInfo.frame_size or similar

% sumFrame: -> results

% tol: -> settings





settings=struct(...
    'Thresh',[],...
    'allowSig',[],...
    'allowX',[],...
    'estNum',[],...
    'initSig',[],...
    'initX0',[],...
    'maxIts',[],...
    'prevSF',[],...
    'rad',[],...
    'tol',[]...
    );



results=struct(...
    'SupResParams',[],...
    'numberOfFrames',[],...
    'numberOfLocs',[],...
    'sumFrame',[]...
    );



localization=struct(...
    'algo_id',[],...
    'settings',settings,...
    'results',results);

%localization = rainSTORM_createLocStruct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reviewer



% SparrowThompsonLimit: results

% SupResDeltaX: -> localization.results.SupResParams

% SupResIm: results

% SupResParams: -> localization.results.SupResParams

% SupResPosits: -> localization.results.SupResParams

% SavedSupResParams:

% SavedSupResPosits:

% alg: -> localization.algo_id

% algVisual: settings, replace to a string as like localization.algo_id

% boxCols: settings

% boxRows: settings

% boxedColour: settings

% boxedDeltaX: settings

% boxedParams: settings

% boxedPosits: settings

% countsPerPhoton: -> rawdata_mgr.countsPerPhoton

% densestPoint: results

% figHistsImage: results

% figNewReconHandle: results

% figOffsetPairs: results

% filename: -> rawdata_mgr.filename

% jhImage: results

% jhParams: results

% jhLinMag: results

% linMag: settings

% markerAnchor: results

% markerFrames: results

% markerParams: results

% markerPosits: results

% meanRevDeltaX: results

% myError: rename and store in params.error

% myError2: rename and store in params.error

% myError3: rename and store in params.error

% MError: rename and store in params.error% myFrame: delete, this just the last localized frame, in rainSTORM_reviewer.m replace size(myFrame) with

% newFrames: -> settings.filter_settings

% newPrecision: -> settings.filter_settings

% newSigma: -> settings.filter_settings

% newTol: -> settings.filter_settings

% numberAcceptedPerFrame: results

% numberRejectedPerFrame: results

% numberOfFiles: -> localization.numberOfFrames

% numberOfFits: -> localization.numberOfLocs

% offPairsCH1: params.opoffCh1

% offPairsCH2: params.opoffCh2

% optOffTFORM: results

% pixelWidth: rawdata_mgr.myImInfo.pixelWidth

% prevSF: settings

% qualityApprovedRows: rename to qualityApprovedLocs and store in results

% reconstructionScaleFactor: delete, if used replace with linMag



% reviewedDeltaX: results.reviewedSupResParams accordingly to localization.results.SupResParams

% reviewedParams: results.reviewedSupResParams accordingly to localization.results.SupResParams

% reviewedPhotonNums: results.reviewedSupResParams accordingly to localization.results.SupResParams

% reviewedPosits: results.reviewedSupResParams accordingly to localization.results.SupResParams

% results.reviewedSupResParams shall store the data of the selected localizations

% either by filter_settings of boxtracking



% scaleBarLn: results

% selectedRows: -> boxtrack_params.selectedLocs

% sizeOfCCDFrame: -> rawdata_mgr.myImInfo.frame_size or similar

% stdRevDeltaX: results

% sumFrame: localization.results.SupResParams

% flags: see params.flags

% opoffCh1: see params.opoffCh1

% opoffCh2: see params.opoffCh2



boxtrack_params=struct(...
    'boxCols',[],...
    'boxRows',[],...
    'boxedColour',[],...
    'boxedParams',[],...
    'selectedLocs',[]...
    );



filter_settings=struct(...
    'newFrames',[],...
    'newPrecision',[],...
    'newSigma',[],...
    'newTol',[],...
    'newThresh',[]...
    );



settings=struct(...
    'algVisual',[],...
    'boxtrack_params',boxtrack_params,...
    'linMag',[],...
    'filter_settings',filter_settings,...
    'prevSF',[],...
    'opoff_Nmarkers',[]...
    );



results=struct(...
    'SparrowThompsonLimit',[],...
    'SupResIm',[],...
    'stdRevDeltaX',[],...
    'densestPoint',[],...
    'figHistsImage',[],...
    'figNewReconHandle',[],...
    'figOffsetPairs',[],...
    'jhImage',[],...
    'jhParams',[],...
    'jhLinMag',[],...
    'markerAnchor',[],...
    'markerFrames',[],...
    'markerParams',[],...
    'markerPosits',[],...
    'meanRevDeltaX',[],...
    'numberAcceptedPerFrame',[],...
    'numberRejectedPerFrame',[],...
    'optOffTFORM',[],...
    'qualityApprovedLocs',[],...
    'reviewedSupResParams',[],...
    'scaleBarLn',[]...
    );



others=struct(...
    );



reviewer=struct(...
    'settings',settings,...
    'results',results,...
    'error',[],...
    'others',[]...
);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flags and others

flags=struct('Boxed',false,...
    'CalculatedJH',false,...
    'HistsPlotted',false,...
    'LocsPerFrame',false,...    %other variable
    'MarkAnchor',false,...
    'OpOffFound',false,...      %other variable
    'OpOffCorrected',false,...  %other variable
    'SB',false,...
    'Saved',false,...
    'SavedSupResData',false,...
    'SelectedAll',false,...
    'Sum',false,...
    'SomeCuts',false,...        %other variable
    'Undrifted',false);



opoffCh1=struct(...
    'reviewedSupResParams',[],...
    'offPairs',[],...
    'SupResIm',[]);



opoffCh2=struct(...
    'reviewedSupResParams',[],...
    'offPairs',[],...
    'SupResIm',[]);



st = struct(...
    'rawdata_mgr',rawdata_mgr,...
    'localization',localization,...
    'reviewer',reviewer,...
    'error',[],...
    'flags',flags,...
    'opoffCh1',opoffCh1,...
    'opoffCh2',opoffCh2...
    );

end

















