clear;
clc;

Name = ['bn-dmpc-top-YH'];
outputFolder = ['output-',Name];
VideoName = ['MOV_',Name];

imageNames = dir(fullfile(outputFolder,'*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter([VideoName,'.mp4'],'MPEG-4');
outputVideo.FrameRate = 60;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(outputFolder,imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)