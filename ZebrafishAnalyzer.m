classdef ZebrafishAnalyzer
    
    properties
        currentDirectory;
        metaDataContainer;
        imageSeries;
        imageMetaData;
        currentFrame;
        %colorMap;
    end
    
    methods
        
        function obj = ZebrafishAnalyzer(dir, filename)
            currentDirectory = dir;
            fid = fopen(strcat(dir,'positionData.txt'),'r');
            ifile = 0;
            while feof(fid)==0
                metaDataContainer(ifile,1) = fgets(fid);
                tmp = fscanf(fid,'p=(%f,%f,%f)\n');
                metaDataContainer(ifile,2) = tmp(1);
                metaDataContainer(ifile,3) = tmp(2);
                metaDataContainer(ifile,4) = tmp(3);
                ifile = ifile + 1;
            end
            imageSeries = Tiff(strcat(dir,filename),'r');
            for i=1:size(metaDataContainer(:,1)
                metaDataFound = 0;
                if strcmp(filename,metaDataContainer(i,1))==1
                    imageMetaData(1) = metaDataContainer(2);
                    imageMetaData(2) = metaDataContainer(3);
                    imageMetaData(3) = metaDataContainer(4);
                    metaDataFound = 1;
                    break;
                end
            end
            currentFrame = imageSeries.read();
        end
        
        function obj = loadImage(filename)
            imageSeries = Tiff(strcat(currentDirectory,filename),'r');
            for i=1:size(metaDataContainer(:,1)
                metaDataFound = 0;
                if strcmp(filename,metaDataContainer(i,1))==1
                    imageMetaData(1) = metaDataContainer(2);
                    imageMetaData(2) = metaDataContainer(3);
                    imageMetaData(3) = metaDataContainer(4);
                    metaDataFound = 1;
                    break;
                end
            end
            currentFrame = imageSeries.read();
        end
        
        function obj = loadFrame(index)
            imageSeries.setDirectory(index);
            currentFrame = imageSeries.read();
        end
        
    end
    
end