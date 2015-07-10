# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:11:09 2015

@author: asood
"""

import sys
sys.path.append('/Users/asood/Documents/neuroscience/nand/CalBlitz/')
import numpy as np
from XMovie import XMovie
import time
import pylab as plt
plt.ion()
import copy
import exifread
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import label
from math import sqrt

stimulusConditions = range(21)

stimSequence = np.zeros(90,dtype=np.uint32)

for i in range(8):
    stimSequence = np.concatenate((stimSequence,2**(i+1)*np.ones(30,dtype=np.uint32)))

stimSequence = np.concatenate((stimSequence,np.zeros(90,dtype=np.uint32)))

for i in range(8):
    stimSequence = np.concatenate((stimSequence,2**(i+9)*np.ones(30,dtype=np.uint32)))
    
stimSequence = np.concatenate((stimSequence,np.zeros(90,dtype=np.uint32)))

stimSequence = np.concatenate((stimSequence,2**17*np.ones(6,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**18*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**17*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**18*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**17*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**18*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**17*np.ones(9,dtype=np.uint32)))

stimSequence = np.concatenate((stimSequence,np.zeros(90,dtype=np.uint32)))

stimSequence = np.concatenate((stimSequence,2**20*np.ones(6,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**19*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**20*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**19*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**20*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**19*np.ones(15,dtype=np.uint32)))
stimSequence = np.concatenate((stimSequence,2**20*np.ones(9,dtype=np.uint32)))

stimSequence = np.concatenate((stimSequence,np.zeros(90,dtype=np.uint32)))

class ZebraFishTools(object):

    def __init__(self,files=None,nchan=1):
        #For holding raw data
        self.filenames = []
        self.data = []
        self.stimData = []
        self.scopePos = []
        #To be filled when analyzing a particular data file
        self.movieCenter = [0,0,0]
        #self.stimLabels = None
        #self.motionCorrectedMovie = None
        #self.components = None
        self.noise = 0
        self.neuronMasks = []
        self.neuronFs = []
        self.neuronDFFs = []
        self.neuronCentroids = []
        self.neuronSizes = []
        if not files == None:
            for phil in files:
                self.filenames.append(phil)
                channels = self.splitChannels(XMovie(phil,frameRate=0.33),nchan)
                self.data.append(channels[0])
                if nchan > 1: self.stimData.append(channels[1])
                self.scopePos.append(self.getPositionData(phil))
        

    def addDataFile(self,phil=None,nchan=1):
        if phil == None:
            raise Exception('file name not provided')
        self.filenames.append(phil)
        channels = self.splitChannels(XMovie(phil,frameRate=0.33),nchan)
        self.data.append(channels[0])
        if nchan > 1: self.stimData.append(channels[1])
        self.scopePos.append(self.getPositionData(phil))

    def getPositionData(self,filename):
        f = open(filename,'rb')
        tags = exifread.process_file(f)
        pos = [0,0,0]
        for item in tags['IFD 2 ImageDescription'].values.split('\r'):
            if not 'state.motor.rel' in item:
                continue
            if 'state.motor.relXPosition' in item:
                pos[0] = item.split('=')[1]
            elif 'state.motor.relYPosition' in item:
                pos[1] = item.split('=')[1]
            elif 'state.motor.relZPosition' in item:
                pos[2] = item.split('=')[1]
        self.scopePos.append(pos)
    
    def splitChannels(self,compMov=None,nchan=1):
        if compMov == None:
            raise Exception('Need to provide movie to split')
        frameList = range(np.shape(compMov.mov)[0])
        channels = []
        for ich in range(nchan):
            channels.append(compMov.makeSubMov(frameList[ich::nchan]))
        return channels

    def labelStimulusConditions(self,stimMov=None,filename=None,index=-1):
        if not filename == None:
            stimMov = self.stimData[self.filenames.index(filename)]
        elif not index < 0:
            stimMov = self.stimData[index]
        if stimMov == None:
            raise Exception('Must provide stimulus movie, filename, or index')
        totalPhotons = stimMov.mov.sum(axis=1).sum(axis=1)
        meanPhotons = totalPhotons.mean()
        runningProduct = 1
        framesRunning = 0
        alignToFrame = -1
        for i in range(len(totalPhotons)):
            runningProduct = runningProduct*(totalPhotons[i]>meanPhotons)
            if runningProduct == 1:
                framesRunning = framesRunning + 1
                if framesRunning > 29:
                    alignToFrame = i - 29
                    break
            else:
                framesRunning = 0
                runningProduct = 1
        startingStim = 330 - alignToFrame
        startingStim = startingStim - len(stimSequence)*int(startingStim/len(stimSequence))
        alignedStimulus = stimSequence[startingStim:]
        alignedStimulus = np.concatenate((alignedStimulus,stimSequence[0:startingStim]))
        nframes = np.shape(stimMov.mov)[0]
        while len(alignedStimulus) < nframes:
            alignedStimulus = np.concatenate((alignedStimulus,stimSequence[startingStim:]))
            alignedStimulus = np.concatenate((alignedStimulus,stimSequence[0:startingStim]))
        self.stimLabels = alignedStimulus[0:nframes]
        
    def motionCorrect(self,filename=None,index=-1):
        if not filename == None:
            index = self.filenames.index(filename)
        elif index < 0:
            raise Exception('Must provide filename or index')
        m = copy.copy(self.data[index])
        templates=[];
        shifts=[];
        max_shift=5;
        num_iter=3;
        for j in range(0,num_iter):
            template,shift=m.motion_correct(max_shift=max_shift,template=None,show_movie=False);
            templates.append(template)
            shift=np.asarray(shift)
            shifts.append(shift)
            #plt.plot(np.asarray(shifts).reshape((j+1)*shift.shape[0],shift.shape[1]))
            #plt.show()
        m.crop(max_shift,max_shift,max_shift,max_shift)
        self.motionCorrectedMovie = m
        self.movieCenter = self.scopePos[index]
        
    def doPCAICA(self,ncomp=50):
        initTime=time.time()
        mdff=copy.copy(self.motionCorrectedMovie)
        mdff.computeDFF(secsWindow=5,quantilMin=20,subtract_minimum=True)
        print 'elapsed time:' + str(time.time()-initTime)
        initTime=time.time()
        self.components=mdff.IPCA_stICA(components=50);
        print 'elapsed time:' + str(time.time()-initTime)
        
    def findNeuronsFromICs(self):
        masks=self.motionCorrectedMovie.extractROIsFromPCAICA(self.components, numSTD=8, gaussiansigmax=2 , gaussiansigmay=2)
        nframes,h,w = np.shape(self.motionCorrectedMovie.mov)
        flatMovie = np.reshape(self.motionCorrectedMovie.mov,(nframes,h*w))
        #calculate noise using pixels not in ROIs
        noiseMask = np.asarray(masks[0])
        for i in xrange(1,len(masks)):
            noiseMask = noiseMask + np.asarray(masks[i])
        noiseMask = (noiseMask < 1)
        noisePix = np.sum(noiseMask)
        if noisePix < 10000:
            print 'Warning: Small number of pixels used for noise calculation'
        flatNoiseMask = np.reshape(noiseMask, (1,h*w))
        noiseFTrace = np.dot(flatNoiseMask,np.transpose(flatMovie)) / float(noisePix)
        self.noise = np.std(noiseFTrace) / np.mean(noiseFTrace)
        #cycle through each potential ROI
        for i in range(np.shape(masks)[0]):
            for j in range(np.max(masks[i])):
                tmpMask = (np.asarray(masks[i]) == j+1)
                tmpSize = float(np.sum(tmpMask))
                #reject if too small or too large
                if tmpSize < 30 or tmpSize > 1000:
                    continue;
                #require range of fluorescence signal to be greater than noise
                flatMask = np.reshape(tmpMask,(1,h*w))
                tmpFTrace = np.dot(flatMask,np.transpose(flatMovie)) / tmpSize
                tmpFTrace = tmpFTrace[0]
                #if np.max(tmpFTrace) - np.min(tmpFTrace) < 3*self.noise:
                if np.std(tmpFTrace) / np.mean(tmpFTrace) < 2*self.noise:
                    continue
                #keeping this ROI, calculate other quanties of interest
                #first, ROI center and convert position/size from pixels to microns
                pixels = np.where(np.asarray(tmpMask) == 1)
                tmpCentroid = [np.sum(pixels[1])/tmpSize,np.sum(pixels[0])/tmpSize,float(self.movieCenter[2])]
                tmpCentroid[0] = (tmpCentroid[0] - w/2)*0.5 + float(self.movieCenter[0])
                tmpCentroid[1] = (tmpCentroid[1] - h/2)*0.5 + float(self.movieCenter[1])
                tmpSize = tmpSize * 0.25
                #finally, dF/F
                window=int(10/self.motionCorrectedMovie.frameRate)
                minQuantile=20
                traceBL=[np.percentile(tmpFTrace[k:k+window],minQuantile) for k in xrange(1,len(tmpFTrace)-window)]
                missing=np.percentile(tmpFTrace[-window:],minQuantile);
                missing=np.repeat(missing,window+1)
                traceBL=np.concatenate((traceBL,missing))
                tmpDFFtrace = (tmpFTrace-traceBL)/traceBL
                self.neuronMasks.append(tmpMask)
                self.neuronSizes.append(tmpSize)
                self.neuronCentroids.append(tmpCentroid)
                self.neuronFs.append(tmpFTrace)
                self.neuronDFFs.append(tmpDFFtrace)
    
    def findNeuronsFromStDev(self,applyFilter=True,sigmax=3,sigmay=3,tRange=None,sizeLimit=40,distanceLimit=4.5):
        stdFrame = np.zeros(np.shape(self.motionCorrectedMovie.mov[0]))
        h,w = np.shape(stdFrame)
        for r in range(h):
            for c in range(w):
                stdFrame[r][c] = np.std(self.motionCorrectedMovie.mov[:,r,c])
        if applyFilter:
            stdFrame = gaussian_filter(stdFrame,[sigmax,sigmay])
        if tRange == None:
            med = np.median(stdFrame)
            stdStd = np.std(stdFrame)
            tRange = [med + 5*stdStd, med + 4*stdStd, med + 3*stdStd, med + 2*stdStd, med + stdStd, med]
        conmat = np.ones((3,3)) #np.array([[0,1,0],[1,1,1],[0,1,0]])
        ntot = 0
        for threshold in tRange:
            stdFrameFree = stdFrame
            for m in self.neuronMasks:
                stdFrameFree = stdFrameFree * (1-(m>0))
            mask = stdFrameFree*(stdFrameFree > threshold)
            mask, n = label(mask > 0, conmat)
            nRejected = 0
            nPruned = 0
            avgPixelsPruned = 0
            for i in range(1,n+1):
                iROI = (mask == i)
                iSize = np.sum(iROI)
                stdROI = iROI * stdFrameFree
                [[rmax],[cmax]] = np.where(stdROI == np.max(stdROI))
                r,c = np.where(iROI == 1)
                for j in range(len(r)):
                    if ((r[j]-rmax)**2 + (c[j]-cmax)**2)**0.5 < distanceLimit:
                        iROI[r[j],c[j]] = 0
                mask = mask * (1-iROI)
                diff = np.sum(iROI)
                iSizePruned = iSize - diff
                if iSizePruned < sizeLimit:
                    iROI = (mask == i)
                    mask = mask * (1-iROI)
                    nRejected = nRejected + 1
                elif diff > 0:
                    nPruned = nPruned + 1
                    avgPixelsPruned = avgPixelsPruned + diff
            avgPixelsPruned = avgPixelsPruned / float(nPruned)
            mask, nfinal = label(mask > 0, conmat)
            self.neuronMasks.append(mask)
            print '%i potential neurons found above threshold %f' % (n,threshold)
            print '%i rejected after pruning for having size < %i' % (nRejected,sizeLimit)
            print '%i kept neurons pruned, average %i pixels removed' % (nPruned,avgPixelsPruned)
            ntot = ntot + n - nRejected
        return ntot
    
    def resolveOverlaps(self):
        nover = 0
        for i in range(len(self.neuronMasks)):
            sumi = np.sum(self.neuronMasks[i])
            if sumi == 0:
                continue
            for j in xrange(i+1,len(self.neuronMasks)):
                overlap = (self.neuronMasks[i] + self.neuronMasks[j] == 2)
                cov = np.mean(np.dot(self.neuronDFFs[i],np.transpose(self.neuronDFFs[j]))) - np.mean(self.neuronDFFs[i])*np.mean(self.neuronDFFs[j])
                cov = cov / (len(self.neuronDFFs[i]) - 1)
                if np.sum(overlap) > 0.5 * np.min([sumi,np.sum(self.neuronMasks[j])]):
                    nover = nover + 1
                    self.neuronMasks[i] = (self.neuronMasks[i] + self.neuronMasks[j] > 0)
                    self.neuronMasks[j] = self.neuronMasks[j] * 0
        print '%i overlaps found' % nover
        self.computeNeuronVars()
        for i in reversed(range(len(self.neuronMasks))):
            if np.sum(self.neuronMasks[i]) == 0:
                del self.neuronMasks[i]
                del self.neuronSizes[i]
                del self.neuronCentroids[i]
                del self.neuronFs[i]
                del self.neuronDFFs[i]
    
    def computeNeuronVars(self):
        nframes,h,w = np.shape(self.motionCorrectedMovie.mov)
        flatMovie = np.reshape(self.motionCorrectedMovie.mov,(nframes,h*w))
        for mask in self.neuronMasks:
            for i in range(1,np.max(mask)+1):
                tmpMask = (np.asarray(mask) == i)
                tmpSize = float(np.sum(tmpMask))
                flatMask = np.reshape(tmpMask,(1,h*w))
                tmpFTrace = np.dot(flatMask,np.transpose(flatMovie)) / tmpSize
                tmpFTrace = tmpFTrace[0]
                pixels = np.where(np.asarray(tmpMask) == 1)
                tmpCentroid = [np.sum(pixels[1])/tmpSize,np.sum(pixels[0])/tmpSize,self.movieCenter[2]]
                #tmpCentroid[0] = (tmpCentroid[0] - w/2)*0.5 + self.movieCenter[0]
                #tmpCentroid[1] = (tmpCentroid[1] - h/2)*0.5 + self.movieCenter[1]
                #tmpSize = tmpSize * 0.25
                #finally, dF/F
                window=int(10/self.motionCorrectedMovie.frameRate);
                minQuantile=20;
                traceBL=[np.percentile(tmpFTrace[k:k+window],minQuantile) for k in xrange(1,len(tmpFTrace)-window)]
                missing=np.percentile(tmpFTrace[-window:],minQuantile);
                missing=np.repeat(missing,window+1)
                traceBL=np.concatenate((traceBL,missing))
                tmpDFFtrace = (tmpFTrace-traceBL)/traceBL
                self.neuronSizes.append(tmpSize)
                self.neuronCentroids.append(tmpCentroid)
                self.neuronFs.append(tmpFTrace)
                self.neuronDFFs.append(tmpDFFtrace)
                
    def clearNeuronVars(self):
        self.neuronSizes = []
        self.neuronCentroids = []
        self.neuronFs = []
        self.neuronDFFs = []
            
    def analyzeData(self,filename=None,index=-1):
        if not filename == None:
            index = self.filenames.index(filename)
        if index < 0:
            raise Exception('Must provided either filename or index')
        print 'Labeling stimuli...'
        self.labelStimulusConditions(index=index)
        print 'Motion correcting...'
        self.motionCorrect(index=index)
        print 'PCA/ICA...'
        self.doPCAICA()
        print 'Finding neurons...'
        self.findNeurons()
        print 'Looking for redundancies...'
        self.resolveOverlaps()
        print 'Done. %i neurons found.' % len(self.neuronMasks)

    
    def makeCompositeMask(self):
        compMask = np.zeros(np.shape(self.neuronMasks[0]))
        for i in xrange(len(self.neuronMasks)):
            addMask = (self.neuronMasks[i] + np.max(compMask)) * (self.neuronMasks[i] > 0)
            compMask = np.add(compMask,addMask)
        return compMask
        
    def plotROIsOverlay(self,vmin=-1,vmax=-1,saveAs=None,customMask=None):
        blankFrame = np.zeros(np.shape(self.motionCorrectedMovie.mov[0]))
        if not customMask == None:
            compMask = copy.copy(customMask)
        else:
            compMask = self.makeCompositeMask()
        h,w = np.shape(blankFrame)
        for r in range(h):
            contourFound = False
            for c in range(w-3):
                blankFrame[r][c] = np.std(self.motionCorrectedMovie.mov[:,r,c])
                if not contourFound:
                    if compMask[r][c]*compMask[r][c+1]>0:
                        contourFound = True
                else:
                    if compMask[r][c+1] == 0:
                        contourFound = False
                    elif compMask[r][c+2]*compMask[r][c+3]>0:
                        compMask[r][c+1] = 0
            for c in range(w-3,w):
                blankFrame[r][c] = np.std(self.motionCorrectedMovie.mov[:,r,c])
        
        compMask = np.ma.masked_where(compMask==0,compMask)

        fig, ax = plt.subplots()
        if vmin > 0 and vmax > 0:
            ax.imshow(blankFrame, cmap=plt.cm.Greys_r,vmin=vmin,vmax=vmax)
        else:
            ax.imshow(blankFrame, cmap=plt.cm.Greys_r)
        ax.imshow(compMask, interpolation='none')
        #mplt.savefig('overlay.png')
        plt.show()
        if not saveAs == None:
            plt.savefig(saveAs)
        
    def plotStDev(self,vmin=-1,vmax=-1):
        blankFrame = np.zeros(np.shape(self.motionCorrectedMovie.mov[0]))
        h,w = np.shape(blankFrame)
        for r in range(h):
            for c in range(w):
                blankFrame[r][c] = np.std(self.motionCorrectedMovie.mov[:,r,c])
        if vmin > 0 and vmax > 0:
            plt.imshow(blankFrame, cmap=plt.cm.Greys_r,vmin=vmin,vmax=vmax)
        else:
            plt.imshow(blankFrame, cmap=plt.cm.Greys_r)
    
    def showNeuronCorrelations(self):
        l = len(self.neuronDFFs)
        corrMat = np.zeros((l,l))
        for i in range(l):
            for j in range(l):
                mx = np.mean(self.neuronDFFs[i])
                my = np.mean(self.neuronDFFs[j])
                cov = np.sum((self.neuronDFFs[i] - mx)*(self.neuronDFFs[j] - my))
                corr = cov/(sqrt(np.sum((self.neuronDFFs[i] - mx)**2))*sqrt(np.sum((self.neuronDFFs[j] - my)**2)))
                corrMat[i,j] = corr
        plt.imshow(corrMat,interpolation='none')
        plt.jet()
        plt.colorbar()
        plt.show()