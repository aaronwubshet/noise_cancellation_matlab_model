% ============
% Example MATLAB script to extract frequency response data and mess with it
% with defined transfer functions
% ============

clear;

% ============
% User parameters
============
sweepfname = 'data.csv';  % Filename of frequency sweep dataset ***TODO*** You'll need to edit this

% ============
% Read CSV file and average sweeps together
% ============
M = csvread(sweepfname,1,0);

% Below is based on how your CSV is structured:
Frequency = M(:,1);
Magg = M(:,2);
Phasee = M(:,3);
Phasee(Phasee>20) = Phasee(Phasee>20) - 360;
Magc = M(:,4);
Phasec = M(:,5);

%below then makes sure there's no replicates of frequencies (originally
%designed for an input file that generated repeated frequency sweeps (like
%in Lab 04/05/06 when we were characterizing the coil)

freq = sort(unique(Frequency)); % Pull out one copy of each frequency
count = zeros(size(freq)); 
mag = zeros(size(freq));
phase = zeros(size(freq));
magcommand = zeros(size(freq));
phasecommand = zeros(size(freq));

freqspot = containers.Map;
for n = 1:size(freq)
    freqspot(num2str(freq(n)))=n;
end

for n = 1:size(Frequency)
    tempn = freqspot(num2str(Frequency(n)));
    count(tempn)=count(tempn)+1;
    mag(tempn)=mag(tempn)+ Magg(n);
    phase(tempn)=phase(tempn)+Phasee(n);
    magcommand(tempn) = magcommand(tempn)+Magc(n);
    phasecommand(tempn) = phasecommand(tempn)+Phasec(n);
end

for n = 1:size(freq)
    mag(n) = mag(n)/count(n);
    phase(n) = phase(n)/count(n);
    magcommand(n) = magcommand(n)/count(n);
    phasecommand(n) = phasecommand(n)/count(n);
end


% ============
% Convert to complex frequency response and create idfrd object
% ============
y_resp = mag.*exp(phase*pi/180*1i);   % Y is the output (angle sensor measurement)


%Below is the good stuff: This is where we create a control object:
%y_resp is basically your bode plot data expressed as single complex values
%rather than mag and phase separated out.

%SPEAKER AND MIC TF
MS = idfrd(y_resp,freq*2*pi,0);    % create transfer function called 

%Gain
K = 10;

%30Hz low pass
R = 5305.16;
C = 1e-6;
[mag, phase] = bode(idtf([0 1], [R*C 1]), freq*2*pi);
lpfresponse = mag.*exp(1j*phase*pi/180);
lpf = idfrd(lpfresponse,2*pi*freq,0);

%300 Hz low pass
R = 5305.16;
C = .1e-6;
[mag, phase] = bode(idtf([0 1], [R*C 1]), freq*2*pi);
lp1fresponse = mag.*exp(1j*phase*pi/180);
lpf1 = idfrd(lp1fresponse,2*pi*freq,0);

%Lead compensator with phase bump 460 Hz (2890 rad/s)(using 2570 after
%trail and error)
%zero at -1/R1*C and pole at -(R1+R2)/(R1*R2*C) and low frequency gain of R3/(R1+R2)
C = 1e-6;
R1 = 1000;
R2 = 1/(R1*(2570^2 * C^2) - (1/R1));
R3 = (R1+R2)*100;
[mag, phase] = bode(idtf([R1*C*R3 R3], [R1*R2*C R1+R2]), freq*2*pi);
lcresponse = mag.*exp(1j*phase*pi/180);
leadcomp = idfrd(lcresponse, 2*pi*freq, 0);


% subplot(2,3,1);
% margin(MS)
% subplot(2,3,2);
% margin(.05*K*MS*lpf*lpf1*leadcomp);
% subplot(2,3,3);
% margin(K*MS*lpf*lpf*lpf1*lpf1*leadcomp);
% subplot(2,3,4);
% bode(lpf);
% subplot(2,3,5);
% bode(lpf1);
% subplot(2,3,6);
% bode(leadcomp);
%margin(lpf*lpf1*leadcomp*K*MS);


%margin(K*MS*leadcomp*lpf1);
%bode(lpf);
%bode(lpf*leadcomp*lpf*K*MS);
%margin(lpf1*lpf1*leadcomp*K*MS*leadcomp);
shg;