classdef ZebrafishAnalyzer
    
    properties
        directoryList;
        currentDirectory;
        fileMetaDataIndex;
        posDataContainer;
        imageSeries;
        imagePosData;
        currentFrame;
        %colorMap;
    end
    
    methods
        
        function obj = ZebrafishAnalyzer(dir, filename)
            obj.directoryList(1,:) = char(dir);
            obj.currentDirectory = dir;
            fid = fopen(strcat(dir,'positionData.txt'),'r');
            ifile = 1;
            while feof(fid)==0
                tmp = fgets(fid);
                obj.fileMetaDataIndex(ifile,:) = tmp(1:end-1);
                tmp = fscanf(fid,'p=(%f,%f,%f)\n');
                obj.posDataContainer(ifile,1) = tmp(1);
                obj.posDataContainer(ifile,2) = tmp(2);
                obj.posDataContainer(ifile,3) = tmp(3);
                ifile = ifile + 1;
            end
            obj.imageSeries = Tiff(strcat(dir,filename),'r');
            for i=1:length(obj.fileMetaDataIndex(:,1))
                posDataFound = 0;
                if strcmp(filename,char(obj.fileMetaDataIndex(i,:)))==1
                    char(obj.fileMetaDataIndex(i,:))
                    obj.imagePosData(1) = obj.posDataContainer(i,1);
                    obj.imagePosData(2) = obj.posDataContainer(i,2);
                    obj.imagePosData(3) = obj.posDataContainer(i,3);
                    posDataFound = 1;
                    break;
                end
            end
            obj.currentFrame = obj.imageSeries.read();
        end
        
        function obj = loadImage(filename)
            obj.imageSeries = Tiff(strcat(obj.currentDirectory,filename),'r');
            for i=1:length(obj.fileMetaDataIndex(:,1))
                posDataFound = 0;
                if strcmp(filename,char(obj.fileMetaDataIndex(i,:)))==1
                    obj.imagePosData(1) = obj.posDataContainer(i,1);
                    obj.imagePosData(2) = obj.posDataContainer(i,2);
                    obj.imagePosData(3) = obj.posDataContainer(i,3);
                    posDataFound = 1;
                    break;
                end
            end
            obj.currentFrame = obj.imageSeries.read();
        end
        
        function obj = changeDirectory(dir)
            obj.currentDirectory = dir;
            newDir = 1;
            numDir = length(obj.directoryList(:,1));
            for i=1:numDir
                if strcmp(dir,char(obj.directoryList(i,:)))==1
                    newDir = 0;
                    break;
                end
            end
            if newDir==1
                obj.directoryList(numDir+1,:) = dir;
                ifile = length(obj.fileMetaDataIndex(:,1)) + 1;
                while feof(fid)==0
                    tmp = fgets(fid);
                    obj.fileMetaDataIndex(ifile,:) = tmp(1:end-1);
                    tmp = fscanf(fid,'p=(%f,%f,%f)\n');
                    obj.posDataContainer(ifile,1) = tmp(1);
                    obj.posDataContainer(ifile,2) = tmp(2);
                    obj.posDataContainer(ifile,3) = tmp(3);
                    ifile = ifile + 1;
                end
            end
        end
        
        function obj = loadFrame(index)
            obj.imageSeries.setDirectory(index);
            obj.currentFrame = obj.imageSeries.read();
        end
        
    end
    
end