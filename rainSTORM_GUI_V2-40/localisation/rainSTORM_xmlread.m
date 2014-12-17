% rainSTORM_xmlread
% by Mark Deimund
%
% This routine determines the image size and number of frames for .raw data
% from xml data.
%
%
function [ImSzX,ImSzY,numberOfFrames, xmlstruct] = rainSTORM_xmlread(filename)

readXml  = [filename '.xml'];

xmlstruct = xml2struct(readXml);
xmax = xmlstruct.acquisition_metadata.Camera_settings.Subimage__binning.x_max;
xmin = xmlstruct.acquisition_metadata.Camera_settings.Subimage__binning.x_min;
ymax = xmlstruct.acquisition_metadata.Camera_settings.Subimage__binning.y_max;
ymin = xmlstruct.acquisition_metadata.Camera_settings.Subimage__binning.y_min;
numberOfFrames = xmlstruct.acquisition_metadata.Camera_settings.Acquisition_progress.Saved;

% Calculate image size and number of frames.
ImSzX = xmax - xmin + 1;
ImSzY = ymax - ymin + 1;

end
