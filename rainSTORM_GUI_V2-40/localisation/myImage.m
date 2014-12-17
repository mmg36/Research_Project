classdef myImage 
    %myImage   Image reading class
    % obj = myImage(params) or
    % obj = myImage(params, profile) creates a myImage object associated with the
    % image filename (and ext)

    properties %(SetAccess = private, GetAccess = private)
        file = '';  %Filename
        fullfile = ''; %Full file path
        ext = '';
        type = '';  %File extension
        myImInfo;   %Information about the file
        width = 0;  %The width of the image
        height = 0; %The height of the image
        frameSize = [];
        numofim = 0;%The number of the images in the file
        Frames = [];%The image container
        Profile = [];
    end
    
    %Constructor
     methods
%      function obj = myImage(filename, ext, par)
      function obj = myImage(varargin)
         par = 0; 
         
         if nargin >= 1
             var = varargin{1};
             obj.file = var.filename;
             obj.ext = var.ext;
             obj.fullfile = [var.filename var.ext];
             obj.Profile.transpose = true;
             if nargin == 2
                par = varargin{2};
             end
         else
             error('Not enough parameters!')
         end
          
          
         if exist(obj.fullfile, 'file') == 2
            
  
            switch obj.ext
                case {'.tif','.tiff'}
                    obj.myImInfo = imfinfo(obj.fullfile);
                    obj.numofim= numel(obj.myImInfo);
                    obj.width=obj.myImInfo(1).Width; 
                    obj.height=obj.myImInfo(1).Height; 
                    obj.type = 'tif';
                    obj.frameSize= [obj.height,obj.width];
                    
                    TifLink = Tiff(obj.fullfile, 'r');
                    obj.Frames = zeros(obj.width,obj.height,obj.numofim, 'uint16');
% TODO: Parallel Tif read?
%                     if (par)
%                       
%                     else
                        for lpIm = 1:obj.numofim
                            TifLink.setDirectory(lpIm);
                            obj.Frames(:,:,lpIm) = TifLink.read();
                        end
%                     end
                    TifLink.close();
                case '.raw'
                    obj.type='raw';
                    [obj.height,obj.width,obj.numofim] = rainSTORM_xmlread(obj.file);
                    obj.frameSize = [obj.height,obj.width];
                    
                    
                    if (par)
                        myFrameSize = obj.frameSize;
                        myHeight = obj.height;
                        myWidth = obj.width;
                        myFile = obj.file;
                        myExt = obj.ext;
                        myFrames = [];
                        parfor lpIm = 1:obj.numofim; 

                             frameOffset = (lpIm - 1) * myHeight * myWidth * 2;

                             fileID = fopen([myFile myExt],'r'); % Open the file for reading
                             fseek( fileID, frameOffset, 'bof' );
                             frame = fread(fileID,myFrameSize,'uint16'); % X,Y order - NEEDS CHECK
                             if obj.Profile.transpose
                                 tframe = frame'; % Image orientation ???????
                             else
                                 tframe = frame;
                             end
                             myFrames(:,:,lpIm) = tframe;
                            fclose(fileID);
                        end
                        obj.Frames = myFrames;
                    else
                        fileID = fopen([obj.file obj.ext],'r'); % Open the file for reading

                        for lpIm = 1:obj.numofim; 

                            frameOffset = (lpIm - 1) * obj.height * obj.width * 2;

                            fseek( fileID, frameOffset, 'bof' );
                            frame = fread(fileID,obj.frameSize,'uint16'); % X,Y order - NEEDS CHECK
                            obj.Frames(:,:,lpIm)=frame'; % Image orientation ???????
                        end
                        fclose(fileID);
                    end
                otherwise
                    error('Not supported file type:\n%s', obj.type);
            end
         else
            % File does not exist.
            error('File does not exist:\n%s', obj.file);
         end
      end % myImage
     end
   
    methods
        function filename = getFilename(obj)
         filename = obj.file;
        end % Filename get function
        
        function width = getWidth(obj)
         width = obj.width;
        end % Filename get function
        
        function height = getHeight(obj)
         height = obj.height;
        end % Filename get function
        
        function type = getType(obj)
         type = obj.type;
        end % Filename get function
        
        function info = getInfo(obj)
         info = obj.myImInfo;
        end % Filename get function
         
        function numofimg = getNumberOfImages(obj)
         numofimg = obj.numofim;
        end % Filename get function
        
        function frameSize = getFrameSize(obj)
            frameSize = obj.frameSize;
        end
        
        function frame = getFrame(obj, index)
            if (size(obj.Frames,3) < index) || (index <= 0)
                error('Index out of bounds:\n%d\nNumber of images: %d\n', index, obj.numofim);
            else
                frame = obj.Frames(:,:,index);
            end
        end % Frame get function
      
        function disp(obj)
            fprintf(1,'Image name: %s\nImage type: %s\nImage Number: %d\nHeight: %d\nWidth: %d\n',...
            obj.file,obj.type,obj.numofim,obj.height,obj.width);
        end % disp
    end
    
    methods (Access = 'private')
    end
end

