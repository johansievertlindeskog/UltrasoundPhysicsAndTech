
% Generate data
load('PreRF_ImageA.mat')

y = preBeamformed.Signal;
vel = preBeamformed.SoundVel;
ft = preBeamformed.TransmittFreq;
fs = preBeamformed.SampleFreq;

%% Distance calculation
format long
pathlength = vel/fs/2; %distance from central element to line
width = preBeamformed.ElementWidth; %width of element
fullwidth = width*31; %width of element 1:31
timevec = zeros(32,2048,128);
for j=1:128
for k=1:2048
    for i=0:31
        diffpathlength = sqrt((fullwidth-width*i)^2 + (pathlength*k)^2);
        timevec(i+1,k,j) = (pathlength*k + diffpathlength)/vel;
        
    end
end
end
%%
ref = timevec(32,:,:); %The time it takes for the echo to reach the innermost element
ref = ones(32,1,1).*ref;
timedelay1 = (timevec - ref); %The time-delay for each element, using the outer element as ref
timedelay2 = timedelay1(end:-1:1,:,:);
timedelaytot = [timedelay1; timedelay2];
figure
plot(timedelaytot(:,:,1))
%% sample delay matrix with deadzone compensation
deadzonesamples = ceil(3*10^(-3)*fs/vel); %distance / speed * samplefrequency
sampledelay = round(timedelaytot.*fs + deadzonesamples);

%sampledelay(:,1) = sampledelay(:,1) + deadzonesamples;
%% Beamforming of first line
y1 = y(:,:,:);
realY = zeros(2048,64,128);
for j = 1:128
    for k = 1:64
        for i = 1:2048
            yindex = sampledelay(k,i,j);
            if(yindex+i < 2049)
                realY(i,k,j) = y1(yindex+i,k,j);
            end
        end
    end
end
%%
threetotwo = squeeze(sum(realY,2));
threetotwo1 = squeeze(sum(y1,2));
Image2 = abs(hilbert(threetotwo));
Image = abs(hilbert(threetotwo1));
figure
imagesc(Image);
colormap(gray);
figure
imagesc(Image2);
colormap(gray);


