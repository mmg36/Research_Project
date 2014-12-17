% Gausian Fit. Local Minimum Background Subtraction - not perfect! 
% Halt After 3 Misfits (Wolter). 
% Copyright 2012. Refer to 00_license.txt for details.
% See, e.g. http://mathworld.wolfram.com/
function [SupResParams] = rainSTORM_fitLocGF3(frameIdx, myFrame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts)
%function [myFits,myParams] = rainSTORM_fitLocGF3(frameIdx, myFrame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts)
persistent Nfails 

% Determine if main localisation algorith will use old method (subtract
% minimum of ROI) of the new method (subtract average of NON-ROI).
% 0 = old method (Min(ROI)), 1 = new method (Mean(non-ROI))
% fitting (fitting DC term to Gaussian)
FlagState= 1;


% Fit [x0,C,sigX] to make f(x)=C*exp( -(x-x0)^2/(2*sigX^2) ) 
% Work on rows then cols. Reject fits with far-out x0, sigX, or residual.

% myFits = -ones(size(myPixels,1),2);   % [row col] matrix for each candidate
% myParams = -ones(size(myPixels,1),7); % Parameters for accepted fits
SupResParams = [];
index = 1;
Nfails = 0; % Reset the number of rejected fits for this new frame

% New method 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FlagState == 1;
  % ResBackground =(sum(sum(myFrame))-sum(sum(myROI)))/((size(myFrame,1)^2)-(size(myROI,1)^2)); % Background is average count outside ROI.
   InitialBackground = mean (mean(double(myFrame)));
   stdev = std(double(myFrame(:)));
   area = (size(myFrame,1))^2;
   SumBackground = 0; 
   
   for rowNum = 1: size(myFrame,1)
       for colNum = 1 : size(myFrame, 1)
           if myFrame(rowNum, colNum) < 3*stdev+InitialBackground 
               SumBackground = SumBackground + double(myFrame(rowNum, colNum));
           else 
               area = area - 1;
           end;
       end;
   end;
   
   ResBackground = SumBackground / area;
   ResBackground 
   % 'new method'
end            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for lpPx = 1:size(myPixels,1); % For local maxima in descending order
myRow = myPixels(lpPx,1);
myCol = myPixels(lpPx,2);
myROI = myFrame(myRow-rad:myRow+rad,myCol-rad:myCol+rad); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Old method"

if FlagState == 0 ;
    ResBackground = min(myROI(:));  % square region to fit. Subtract minimum.
   % 'old method'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myROI = myROI - ResBackground;
flagRowFits = false;   % Begin by noting the centre-position is not fitted
flagColFits = false;

xx = (-rad:rad)';      % x-positions (rows) (pixel widths) as column vector
yy = (-rad:rad)';      % y-positions (cols) (pixel widths) as column vector
yRows = sum(myROI,2);  % Sum all columns (dim=2) to get row intensities
yCols = sum(myROI,1)'; % Sum all rows (dim=1) in each column.

% Fit a Gaussian to the observed intensities, binned by row:
x0 = initX0;
sigX = initSig;
C  = yRows(rad+1); % Guess height of f(x). Centre value is a good guess.
  fofX = C*exp(-(xx-x0).^2/(2*sigX^2)); %+ D ; % Initial guess of f(x)
  Beta = yRows - fofX; % Change needed in f(x)
  for lpLSF = 1:maxIts
  %if FlagState==2;
      %d = ones  (size (xx,1),1); 
      %A = [fofX/C,fofX.*(xx-x0)/sigX^2,fofX.*(xx-x0).^2/sigX^3,d]; % Jacobian
  %else
   A = [fofX/C,fofX.*(xx-x0)/sigX^2,fofX.*(xx-x0).^2/sigX^3]; % Jacobian
  %end;
  b = A'*Beta;
  a = A'*A;
  dL= a\b;
  C = C + dL(1);
  x0 = x0 + dL(2);
  sigX = sigX + dL(3);
   %if FlagState==2;
     % D = D + 0.01*dL(4); % Include damping term
    %  fofX = C*exp(-(xx-x0).^2/(2*sigX^2)) + D; 
  %else
     fofX = C*exp(-(xx-x0).^2/(2*sigX^2));
  %end;
  Beta = yRows - fofX;
  
 % store iteration variables
 Dfactorx(lpLSF,1,lpPx) = C;
 Dfactorx(lpLSF,2,lpPx) = x0;
 Dfactorx(lpLSF,3,lpPx) = sigX;
 Dfactorx(lpLSF,4,lpPx) = ResBackground;
 Dfactorx(lpLSF,5,lpPx) = sum(Beta.^2); 

    if(abs(x0)>allowX || (sigX < -allowSig(1)) || (sigX > allowSig(2)) )
     break; % Stop iterating if solution drifts too far
    end
  end
  % Judge the fit. Accept if residue is a small proportion of |y^2|, etc.
  residueRows = sum(Beta.^2)/sum(yRows.^2);
  if (residueRows<tol && abs(x0)<allowX && sigX>allowSig(1) && sigX<allowSig(2))
  fitRowPos = double(myRow)+x0-0.5; % Note (-0.5) for image registration
  flagRowFits = true; % Flag the row-direction fit as acceptable
  end
 
  % Fit a Gaussian to the observed intensities, binned by Col
  if(flagRowFits) % Don't fit col-direction if the row-axis fit was rejected
  y0 = initX0;
  sigY = sigX;% Keep sigX from row-fit. It should match the column-fit.
  C = yCols(rad+1);
 % if FlagState==2;
 %   fofX = C*exp(-(yy-y0).^2/(2*sigY^2)) + D; % Initial guess of f(x)
 %   d = ones  (size (xx,1),1);
 % else 
    fofX = C*exp(-(yy-y0).^2/(2*sigY^2));
 % end
  Beta = yCols - fofX; % Change needed in f(x)
  
  for lpLSF = 1:maxIts
      
% if FlagState==2;
%    A = [fofX/C,fofX.*(yy-y0)/sigY^2,fofX.*(yy-y0).^2/sigY^3,d]; % Jacobian
% else
     A = [fofX/C,fofX.*(yy-y0)/sigY^2,fofX.*(yy-y0).^2/sigY^3];
% end
  b = A'*Beta;
  a = (A'*A);
  dL= a\b;
  C = C+dL(1);
  y0 = y0 + dL(2);
   sigY = sigY + dL(3);
 % if FlagState==2;
 %   D = D + dL(4)*0.001;
 %   fofX = C*exp(-(yy-y0).^2/(2*sigY^2))+ D; 
 % else 
    fofX = C*exp(-(yy-y0).^2/(2*sigY^2));
 % end
  Beta = yCols - fofX;
 % Damping factor optimisation, store iteration variables
  Dfactory(lpLSF,1,lpPx) = C;
  Dfactory(lpLSF,2,lpPx) = y0;
  Dfactory(lpLSF,3,lpPx) = sigY;
  Dfactory(lpLSF,4,lpPx) = ResBackground;
  Dfactory(lpLSF,5,lpPx) = sum(Beta.^2); 
 
  if(abs(y0)>allowX || (sigY < allowSig(1)) || (sigY > allowSig(2)) )
     break; % Stop iterating if solution drifts too far
   end
  end
  % Judge the column direction fit:
  residueCols = sum(Beta.^2)/sum(yCols.^2);
   if (residueCols<tol && abs(y0)<allowX && sigY>allowSig(1) && sigY<allowSig(2))
    fitColPos = double(myCol)+y0-0.5;
    flagColFits = true; % Flag the column-direction fit as acceptable
   end
    end % End if, which only tries fitting Col-wise if Row-wise fitted OK
 
 if(flagRowFits && flagColFits )     % Accept iff Row and Col fits are OK
%   myFits(lpPx,:)=[fitRowPos,fitColPos];
%   myParams(lpPx,1)=myPixels(lpPx,3); % Averaged magnitude of this signal
%   myParams(lpPx,2)=(residueRows+residueCols)/2; % Mean critical tol for fit
%   myParams(lpPx,3)=sum(yCols); % Sum of signal (counts) for this fit
%   myParams(lpPx,4)=sigX;  % X-width (sigma, rows, fitted) of this Gaussian
%   myParams(lpPx,5)=sigY;  % Y-width (sigma, cols, fitted) of this Gaussian
%   myParams(lpPx,6)=bkgdSig;  % Background for each ROI

  
  SupResParams(index).frame_idx = frameIdx;
  SupResParams(index).x = fitRowPos;
  SupResParams(index).y = fitColPos;
  SupResParams(index).z = 0;
  SupResParams(index).I = myPixels(lpPx,3); % Averaged magnitude of this signal
  SupResParams(index).I = myPixels(lpPx,3); % Averaged magnitude of this signal    
  SupResParams(index).sig_x = sigX;  % X-width (sigma, rows, fitted) of this Gaussian
  SupResParams(index).sig_y = sigY;  % Y-width (sigma, cols, fitted) of this Gaussian
  SupResParams(index).avg_brigthness = bkgdSig; % Background for each ROI
  SupResParams(index).res = (residueRows + residueCols) / 2; % Mean critical tol for fit
  SupResParams(index).res_Row = residueRows;
  SupResParams(index).res_Col = residueCols;
  SupResParams(index).Sum_signal = sum(yCols); % Sum of signal (counts) for this fit
% if FlagState==2
%     SupResParams(index).Sum_signal = sum(yCols)-7*D; % Sum of signal (counts) for this fit
% else
    SupResParams(index).Sum_signal = sum(yCols); % Sum of signal (counts) for this fit
%  end  
  % Print D
  %D
  
  index=index+1;
 else 
  Nfails = Nfails+1;
 end
 if( Nfails>2 )
     break % After 3 failed fits, finish this frame (as Wolter-2010).
 end
 
end  % Loop to the next local maximum

% myParams = myParams(myFits(:,1)~=-1,:); %Params of accepted fits only
% myFits = myFits(myFits(:,1)~=-1,:); % Return accepted fits only. (-1)s are rejected fits.

end