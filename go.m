%% Initialization
clear;
init;

%% Phased array 
harray = phased.URA(...
                    'Element', hant, ...
                       'Size', [30 30], ...
             'ElementSpacing', [lambda/2, lambda/2]);
%figure;
pattern(harray, fc, 'PropagationSpeed', v, 'Type', 'powerdb');
%% Radiator and collector

hradiator.Sensor = harray;
hcollector.Sensor = harray;
hradiator.WeightsInputPort = true;

%% Recalculate transmit power

hag = phased.ArrayGain('SensorArray', harray,'PropagationSpeed', v);
ag = step(hag, fc, [0; 0]);

snr_min = albersheim(pd, pfa, int_pulsenum);
peak_power = radareqpow(lambda, max_distance_range, snr_min, hwav.PulseWidth, ...
    'RCS', tgt_rcs, 'Gain', htx.Gain + ag);

htx.PeakPower = peak_power;

%% Scangrid

%volumnAz = initialAz - endAz;
scangrid = initialAz+scanstep/2:scanstep:endAz;
numscans = length(scangrid);    
pulsenum = int_pulsenum*numscans;

%% Target Definition

targetDefinition;
numtargets = length(htarget);

%% Pulse Synthesis

hstvec = phased.SteeringVector(...
                       'SensorArray', harray, ...
                  'PropagationSpeed', v);
       
hbeamform = phased.PhaseShiftBeamformer(...
                       'SensorArray', harray,...
                'OperatingFrequency', fc, ...
                  'PropagationSpeed', v,...
                   'DirectionSource', 'Input port');

for n = numtargets:-1:1
    htargetchannel{n} = phased.FreeSpace(...
                        'SampleRate', fs,...
                 'TwoWayPropagation', true,...
                'OperatingFrequency', fc);
    htargetpath{n} = zeros(1,2,pulsenum);
end

fast_time_grid = unigrid(0, 1/fs, 1/prf, '[)');
rx_pulses = zeros(numel(fast_time_grid), pulsenum);  % Result matrix
tgt_ang = zeros(2,numtargets);                       % Target angle buffer

for m = 1:pulsenum

    x = step(hwav);                              % Generate pulse
    [s, tx_status] = step(htx,x);                % Transmit pulse
    [ant_pos,ant_vel] = step(hantplatform,1/prf);% Update antenna position
    
    % Calculate the steering vector
    scanid = floor((m-1)/int_pulsenum) + 1;      %повторяется 10 раз
    sv = step(hstvec,fc,scangrid(scanid));       %здесь можно поиграться
    w = conj(sv);
    
    rsig = zeros(length(s),numtargets);
    for n = numtargets:-1:1                      % For each target
        [tgt_pos,tgt_vel] = step(htargetplatform{n}, 1/prf);           % Update target position
        %htargetpath{n}(:,:,m) = [tgt_pos(1) tgt_pos(2)];
        [~,tgt_ang(:,n)] = rangeangle(tgt_pos,...% Calculate range/angle
            ant_pos);                               
        tsig = step(hradiator,s,tgt_ang(:,n),w); % Radiate toward target 
        tsig = step(htargetchannel{n},...        % Propagate pulse
            tsig,ant_pos,tgt_pos,ant_vel,tgt_vel);                  
        rsig(:,n) = step(htarget{n},tsig);       % Reflect off target
    end
    rsig = step(hcollector,rsig,tgt_ang);        % Collect all echoes
    rsig = step(hrx, rsig,~(tx_status>0));       % Receive signal
    rsig = step(hbeamform,rsig,[scangrid(scanid);0]);  % Beamforming
    rx_pulses(:, m) = rsig;                       % Form data matrix
end

figure('Name', 'rx_pulses');
image(abs(rx_pulses)*10e7);

%% Matched Filter

matchingcoeff = getMatchedFilter(hwav);
hmf = phased.MatchedFilter(...
                'Coefficients', matchingcoeff,...
              'GainOutputPort', true);
[mf_pulses, mfgain] = step(hmf, rx_pulses);
mf_pulses = reshape(mf_pulses,[], int_pulsenum, numscans);

matchingdelay = size(matchingcoeff, 1)-1;
sz_mfpulses = size(mf_pulses);
mf_pulses = [mf_pulses(matchingdelay+1:end) zeros(1, matchingdelay)];

mf_pulses = reshape(mf_pulses, sz_mfpulses);
% Pulse integration
int_pulses = pulsint(mf_pulses, 'noncoherent');
int_pulses = squeeze(int_pulses);

figure('Name', 'int_pulses');
image(abs(int_pulses)*10e6);

% Visualize
r = v*fast_time_grid/2;  %здесь вылезли 5000
X = r'*cosd(scangrid); Y = r'*sind(scangrid);
clf;
pcolor(X,Y,pow2db(abs(int_pulses).^2));
axis equal tight 
shading interp
axis off
text(-800,0,'Array');
text((max(r)+10)*cosd(initialAz),(max(r)+10)*sind(initialAz),...
    [num2str(initialAz) '^o']);
text((max(r)+10)*cosd(endAz),(max(r)+10)*sind(endAz),...
    [num2str(endAz) '^o']);
text((max(r)+10)*cosd(0),(max(r)+10)*sind(0),[num2str(0) '^o']);
colorbar;

%% Vizualization
% We now visualize the detection process. To better represent the data, we
% only plot range samples beyond 50.

N = 51;
clf;
surf(X(N:end,:),Y(N:end,:),...
            pow2db(abs(int_pulses(N:end,:)).^2));
hold on;
% mesh(X(N:end,:),Y(N:end,:),...
%     pow2db(threshold*ones(size(X(N:end,:)))),'FaceAlpha',0.8);
view(0,56);
axis off