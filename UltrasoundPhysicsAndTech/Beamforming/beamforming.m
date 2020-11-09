% %% Beamforming assignment
% 
% load 'PostRF_Phantom.mat'; % stored in PostRF
% postbeamformed = PostRF.Signal;
% image = abs(hilbert(postbeamformed));
% 
% % figure()
% % plot(image)
% 
% % picture of 1 colonn
% figure()
% plot(postbeamformed(:,1));
% 
% figure()
% imagesc(image);
% colormap(gray);
% 
% % periodogram of 1 colonn of demeaned postbeamformed data
% meanPostbeamformed = postbeamformed - mean(postbeamformed);
% P = 2^12;
% ff = (0:P-1)/P-0.5;
% per = fftshift(abs(fft(meanPostbeamformed(:,1), P)).^2);
% 
% % HP filtering
% % % HPfilter = designfilt('highpassfir', 'StopbandFrequency', 1, 'PassbandFrequency', 5, 'StopbandAttenuation', 60, 'PassbandRipple', 1, 'SampleRate', 25);
% % % filteredMeanPostbeamformed = filter(HPfilter, meanPostbeamformed);
% % % imageFiltered = abs(hilbert(filteredMeanPostbeamformed));
% 
% filteredMeanPostbeamformed = highpass(meanPostbeamformed,0.3);
% 
% 
% % After filtering, periodogram of 1 colonn of demeaned postbeamformed data
% perFiltered = fftshift(abs(fft(filteredMeanPostbeamformed(:,1), P)).^2);
% 
% 
% % HP filtering the signal
% % % [B,A] = butter(N,Wn,'high');
% % % postbeamformedFiltered = filtfilt(HPfilter, postbeamformed);
% 
% % periodogram before and after filering
% figure()
% subplot(2,1,1)
% semilogy(ff, per)
% legend('Periodogram of the demeaned postbeamformed data')
% subplot(2,1,2)
% semilogy(ff, perFiltered)
% legend('Periodogram of the filtered demeaned postbeamformed data')
% 
% 
% % figure()
% % plot(image)
% 
% % figure()
% % imagesc(imageFiltered);
% % colormap(gray);

%% Examining a pre-beamformed B-mode image, the 2:nd dimension
%load 'PreRF_ImageB.mat' % stored in preBeamformed
% prebeamformed = preBeamformed.Signal;
% threeToTwo = squeeze(sum(prebeamformed,2));
% image2 = abs(hilbert(threeToTwo));
% 
% figure()
% imagesc(image2);
% colormap(gray)

%% Creating dynamic focus
clear all
clc

% loading a pre beamformed picture
load('PreRF_ImageA.mat') % struct stored in preBeamformed
y = preBeamformed.Signal; % the data
Fs = preBeamformed.SampleFreq; % sampling frequency for the data
speedOfSound = preBeamformed.SoundVel; % 1540 m/s (not km/h!)

%% Distance calculating
format long;
vectical = speedOfSound/Fs/2; % distance from central element to focal point 1, divide by 2 since it is measured for the sound back and forth
elementWidth = preBeamformed.ElementWidth; % the width of every element
fullWidth = 31*elementWidth; % total width of the elements, 1st one calculated as zero
deadzone = 3*10^(-3); % the deadzone in mm


timeVec = zeros(32,2048,128);
for m = 1:128
    for k = 1:2048
        for i = 0:31
            horizontal = fullWidth - i*elementWidth;
            hyp = sqrt((vectical*k + deadzone)^2+(horizontal)^2);
            travelTime = (hyp+vectical*k)/speedOfSound;
            timeVec(i+1,k,m) = travelTime;
        end
    end
end

c = 1.05; % tuning parameter, different for imageA and B
ref = timeVec(32,:,:); % time taken for the echo to reach the innermost element
ref = ones(32,1,1).*ref*c;

timedelay1 = (timeVec-ref); %The time-delay for each element, using the outer element as ref
timedelay2 = timedelay1(end:-1:1,:,:);
timedelayTot = [timedelay1; timedelay2]; % nbr of sample we delay every element with, left end point--> right end point

% plotting the timedelay
% figure()
% plot(timedelayTot(:,:,1))

%%
% compensating for the deadzone
delaySamples = round(timedelayTot*Fs);

% using the index in delaySamples for finding where to look in the data to be beamformed

resultPostBeamforming = zeros(2048,64, 128);
for m = 1:128
    for k = 1:64
        for i = 1:2048
            yindex = delaySamples(k,i,m);
            if (yindex + i < 2049)
                resultPostBeamforming(i,k,m) = y(yindex+i,k,m);
            end
        end
    end
end

%% plotting results

threetotwo = squeeze(sum(resultPostBeamforming,2));
threetotwoComp = squeeze(sum(y,2));
Image2 = abs(hilbert(threetotwo));
Image = abs(hilbert(threetotwoComp));

figure()
imagesc(Image);
colormap(gray);

figure()
imagesc(Image2);
colormap(gray);


