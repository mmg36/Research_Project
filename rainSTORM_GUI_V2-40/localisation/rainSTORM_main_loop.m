% function params = loop(Frame,params)
function SupResParams = rainSTORM_main_GUI_loop(frameIdx,Frame,initX0, initSig, rad, tol, Thresh, maxIts, alg, allowSig, allowX)
%function [myFits, myParams] = rainSTORM_main_GUI_loop(frameIdx,Frame,initX0, initSig, rad, tol, Thresh, maxIts, alg, allowSig, allowX)

%         myFits = []; 

        % 1. Identify local maxima above threshold, and return [row,col,magnitude]
        myPixels = rainSTORM_avNMS(Frame,rad);

        % 2. Thresholding. To remove noise maxima and weak signals.
        myPixels((myPixels(:,3))<Thresh,:)=[]; % Apply threshold
        myPixels = flipud(sortrows(myPixels,3)); % Sort for human-readability

        bkgdSig = std(double(Frame(Frame < mean(Frame(:))))); % Avoids signal

        % 3. Localise centre of each pixellated blur (and reject if not Gaussian)
        %Implement selected image processing algorithm.
        if alg==1
            SupResParams = rainSTORM_fitLocGF3(frameIdx, Frame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts);
%           [myFits,myParams] = rainSTORM_fitLocGF3(frameIdx, Frame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts);

        elseif alg==2
            % Least squares Gaussian Fitting without halting
            SupResParams = rainSTORM_fitLocGF(frameIdx, Frame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts);
%            [myFits,myParams] = rainSTORM_fitLocGF(frameIdx, Frame,myPixels,initX0,initSig,allowSig,rad,tol,allowX,bkgdSig,maxIts);

        elseif alg==3 % Or fit by Centre of Mass: find 1st+2nd moment of intensity
            SupResParams = rainSTORM_fitCoM(frameIdx, Frame,myPixels,rad,bkgdSig); 
%            [myFits,myParams] = rainSTORM_fitCoM(frameIdx, Frame,myPixels,rad,bkgdSig); 
        end
end
