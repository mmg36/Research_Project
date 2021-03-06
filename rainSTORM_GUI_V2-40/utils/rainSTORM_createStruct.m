function st = rainSTORM_createStruct()

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

    st = struct('parCompTool', parCompTool,...
                    'imgProcTool', imgProcTool, ...
                    'SparrowThompsonLimit',[],...
                    'SupResDeltaX',[],...
                    'SupResConf',[],...                     %other variable
                    'SupResIm',[],...
                    'SupResImHist',[],...                   %other variable
                    'SupResImPrec',[],...                   %other variable
                    'SupResNumber',[],...
                    'SupResParams',[],...
                    'SupResPosits',[],...
                    'SavedSupResParams',[],...              %other variable
                    'SavedSupResPosits',[],...              %other variable
                    'Thresh',[],...
                    'alg',[],...
                    'algVisual',[],...
                    'allowSig',[],...
                    'allowX',[],...
                    'boxCols',[],...                        %other variable
                    'boxRows',[],...                        %other variable
                    'boxedColour',[],...                    %other variable
                    'boxedDeltaX',[],...                    %other variable
                    'boxedParams',[],...                    %other variable
                    'boxedPosits',[],...                    %other variable
                    'countsPerPhoton',[],...
                    'densestPoint',[],...
                    'estNum',[],...
                    'ext',[],...
                    'figHistsHandle',[],...
                    'figHistsImage',[],...                  %other variable
                    'figNewReconHandle',[],...
                    'figOffsetPairs',[],...                 %other variable
                    'filename',[],...
                    'initSig',[],...
                    'initX0',0,...
                    'jhImage',[],...                        %other variable
                    'jhParams',[],...                       %other variable
                    'jhLinMag',[],...                       %other variable
                    'linMag',[],...
                    'maxIts',[],...
                    'markerAnchor',[],...                   %other variable
                    'markerFrames',[],...                   %other variable
                    'markerParams',[],...                   %other variable
                    'markerPosits',[],...                   %other variable
                    'meanRevDeltaX',[],...
                    'myError',[],...
                    'myError2',[],...                       %other variable
                    'myError3',[],...                       %other variable
                    'MError',[],...                         %other variable
                    'myFrame',[],...
                    'myImInfo',[],...
                    'newFrames',[],...
                    'newPrecision',[],...
                    'newSigma',[],...
                    'newTresh',[],...
                    'newTol',[],...
                    'numberAcceptedPerFrame',[],...         %other variable
                    'numberRejectedPerFrame',[],...         %other variable
                    'numberOfFiles',[],...
                    'numberOfFits',[],...
                    'offPairsCH1',[],...                    %other variable
                    'offPairsCH2',[],...                    %other variable
                    'opoff_Nmarkers',[],...                 %other variable
                    'optOffTFORM',[],...                    %other variable
                    'pixelWidth',[],...
                    'prevSF',[],...
                    'qualityApprovedRows',[],...
                    'rad',[],...
                    'reconstructionScaleFactor',[],...
                    'reviewedDeltaX',[],...
                    'reviewedParams',[],...
                    'reviewedPhotonNums',[],...
                    'reviewedPosits',[],...
                    'scaleBarLn',[],...
                    'selectedRows',[],...                   %other variable
                    'sizeOfCCDFrame',[],...
                    'stdRevDeltaX',[],...
                    'sumFrame',[],...
                    'tol',[],...
                    'flags',...
                        struct('Boxed',false,...
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
                                'Undrifted',false),...
                     'opoffCh1',...
                        struct('reviewedDeltaX',[],...
                                'reviewedParams',[],...
                                'reviewedPosits',[],...
                                'SupResIm',[]),...
                     'opoffCh2',...
                        struct('reviewedDeltaX',[],...
                                'reviewedParams',[],...
                                'reviewedPosits',[],...
                                'SupResIm',[])...
                   );
end