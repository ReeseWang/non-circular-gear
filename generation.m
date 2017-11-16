addpath(genpath('./'))
a = 15;                                         % Center distance
pAngle = 22.5*pi/180;                           % Pressure Angle
module = 1;                                     % Module
pitch = pi * module;                            % Curve pitch
addDist = 0.8 * module;                         % Addendum distance
dedDist = 1 * module;                           % Dedendum distance
posLimiterLeng = 3;                             % Position limiter length
errTol = 1e-4;                                  % Tolerance
angTol = 0.5*pi/180;                            % Angular tolerance, rad
toolRadius = 0.2;                               % Forming tool radius

%% Read input, interpolate, and generate pitch curves.
readfp;
interpAngles = polData(1,1):angTol/pi*180:polData(end,1);
interpAngles = interpAngles';
polDataInterp = [interpAngles, spline(polData(:,1), polData(:,2), interpAngles)];
polDataInterp(:,3) = polDataInterp(:,1)*pi/180;

driverPitch = zeros(size(polDataInterp, 1), 2);
followerPitch = zeros(size(polDataInterp, 1), 2);
pitchLengths = zeros(size(polDataInterp, 1), 1);

driverPitch(1,1) = polDataInterp(1, 2) / (1 + polDataInterp(1, 2));
followerPitch(1,1) = - 1 / (1 + polDataInterp(1, 2));
followerAngle = zeros(size(polDataInterp, 1), 1);
followerRad = 1 / (1 + polDataInterp(1, 2));
for i = 2:size(polDataInterp, 1)
    driverRad = polDataInterp(i, 2) / (1 + polDataInterp(i, 2));
    followerRadPrev = followerRad;
    followerRad = 1 - driverRad;
    driverPitch(i,:) = driverRad * [cos(polDataInterp(i,3)) sin(polDataInterp(i,3))];
    driverArcSquared = sum((driverPitch(i,:) - driverPitch(i-1,:)).^2);
    pitchLengths(i) = pitchLengths(i-1) + sqrt(driverArcSquared);
    followerAngle(i) = followerAngle(i-1) - ...
        acos((followerRad^2 + followerRadPrev^2 - driverArcSquared) / ...
        (2 * followerRad * followerRadPrev));
    followerPitch(i,:) = followerRad * [-cos(followerAngle(i)) 
        -sin(followerAngle(i))];
end

driverPitch = driverPitch * a;
followerPitch = followerPitch * a;
pitchLengths = pitchLengths * a;

driverPitchCl = [
    driverPitch;
    driverPitch(end,:)*(eye(2) + [0 1; -1 0]/norm(driverPitch(end,:))*posLimiterLeng);
    driverPitch(1,:)*(eye(2) + [0 -1; 1 0]/norm(driverPitch(1,:))*posLimiterLeng);
    ];
followerPitchCl = [
    followerPitch;
    followerPitch(end,:)*(eye(2) + [0 -1; 1 0]/norm(followerPitch(end,:))*posLimiterLeng);
    followerPitch(1,:)*(eye(2) + [0 1; -1 0]/norm(followerPitch(1,:))*posLimiterLeng);
    ];
%plot(driverPitch(:,1), driverPitch(:,2), followerPitch(:,1) + a, followerPitch(:,2))
%plotc(followerPitchCl(:,1)+a, followerPitchCl(:,2));
%plotc(driverPitchCl(:,1), driverPitchCl(:,2));

%% Generating rack teeth
zigZagHalfHight = pitch/4/tan(pAngle);
temp = polyout(driverPitchCl(:,1), driverPitchCl(:,2), zigZagHalfHight, 'm');
driverInitShape = [temp{1}{1} temp{2}{1}];
%plotc(driverInitShape(:,1),driverInitShape(:,2));
temp = polyout(followerPitchCl(:,1), followerPitchCl(:,2), zigZagHalfHight, 'm');
followerInitShape = [temp{1}{1} temp{2}{1}];
%plotc(followerInitShape(:,1)+a,followerInitShape(:,2));

driverPitchLeng = sum(sqrt(sum((driverPitch(2:end,:) - driverPitch(1:end-1,:)).^2, 2)));
followerPitchLeng = sum(sqrt(sum((followerPitch(2:end,:) - followerPitch(1:end-1,:)).^2, 2)));
assert(abs(driverPitchLeng - followerPitchLeng) < errTol)
assert(abs(pitchLengths(end) - followerPitchLeng) < errTol)

rackZigZag = ((-pitch/4:pitch/2:driverPitchLeng+pitch/2)'*[1 1]);
for i = 1:size(rackZigZag, 1)
    rackZigZag(i,1) = zigZagHalfHight * (-1)^i;
end
if rackZigZag(end,1) < 0 % Driver should drive the follower at the end
    rackZigZag = [rackZigZag; zigZagHalfHight rackZigZag(end,2) + pitch/2];
end


%% Generating rack shapes
driverClip = [
    -dedDist, rackZigZag(2,2); 
    -dedDist-3*zigZagHalfHight, rackZigZag(2,2); 
    -dedDist-3*zigZagHalfHight, rackZigZag(end,2); 
    -dedDist, rackZigZag(end,2)];
followerClip = [
    dedDist, rackZigZag(1,2); 
    dedDist+3*zigZagHalfHight, rackZigZag(1,2); 
    dedDist+3*zigZagHalfHight, rackZigZag(end,2); 
    dedDist, rackZigZag(end,2)];
driverAdd = [
    addDist, rackZigZag(1,2); 
    zigZagHalfHight, rackZigZag(1,2); 
    zigZagHalfHight, rackZigZag(end,2); 
    addDist, rackZigZag(end,2)];
followerAdd = [
    -addDist, rackZigZag(2,2); 
    -zigZagHalfHight, rackZigZag(2,2); 
    -zigZagHalfHight, rackZigZag(end,2); 
    -addDist, rackZigZag(end,2)];
%if rackZigZag(end,1) < 0
%    driverClip(3:4,2) = rackZigZag(end-1,2);
%    followerAdd(3:4,2) = rackZigZag(end-1,2);
%else
%    followerClip(3:4,2) = rackZigZag(end-1,2);
%    driverAdd(3:4,2) = rackZigZag(end-1,2);
%end

rackZigZag = [
    rackZigZag(1,:) - [0 a+posLimiterLeng];
    rackZigZag;
    rackZigZag(end,:) + [0 a+posLimiterLeng]
    ];
%plot(rackZigZag(:,1)+driverPitch(1,1), rackZigZag(:,2))

driverRack = [
    rackZigZag(1,:) + [3*zigZagHalfHight 0];
    rackZigZag;
    [rackZigZag(1,1) + 3*zigZagHalfHight, rackZigZag(end,2)]
    ];
temp = polyclip(driverRack, driverClip, 'dif');
driverRack = [temp{1}{1} temp{2}{1}];
temp = polyclip(driverRack, driverAdd, 'uni');
driverRack = [temp{1}{1} temp{2}{1}];
%plotc(driverClip(:,1)+driverPitch(1,1), driverClip(:,2));
%plotc(driverRack(:,1)+driverPitch(1,1), driverRack(:,2));

followerRack = [
    rackZigZag(1,:) - [zigZagHalfHight 0];
    rackZigZag;
    [rackZigZag(1,1) - zigZagHalfHight, rackZigZag(end,2)]
    ];
temp = polyclip(followerRack, followerClip, 'dif');
followerRack = [temp{1}{1} temp{2}{1}];
temp = polyclip(followerRack, followerAdd, 'uni');
followerRack = [temp{1}{1} temp{2}{1}];
%plotc(followerClip(:,1)+a+followerPitch(1,1), followerClip(:,2));
%plotc(followerRack(:,1)+a+followerPitch(1,1), followerRack(:,2));

%% Begin cutting!
%driverRack = driverRack + driverPitch(1,:);
%followerRack = followerRack + followerPitch(1,:);
temp = polyclip(driverInitShape, driverRack + driverPitch(1,:), 'dif');
driverProfile = [temp{1}{1} temp{2}{1}];
temp = polyclip(followerInitShape, followerRack + followerPitch(1,:), 'dif');
followerProfile = [temp{1}{1} temp{2}{1}];

f = figure;
ax = axes(f);
set(ax, 'XLim', [-a 2*a], 'YLim', [-a a], 'YLimMode', 'manual', 'DataAspectRatio', [1 1 1])
hold all
dph = fill(driverProfile(:,1), driverProfile(:,2), 'y');
fph = fill(followerProfile(:,1)+a, followerProfile(:,2), 'r');
dpr = plotc(driverRack(:,1), driverRack(:,2));
fpr = plotc(followerRack(:,1)+a, followerRack(:,2));
dpc = plotc(driverPitchCl(:,1), driverPitchCl(:,2), '--');
fpc = plotc(followerPitchCl(:,1)+a, followerPitchCl(:,2), '--');
rectangle('Position', [-1 -1 2 2]*a/10, 'Curvature', [1 1]);
rectangle('Position', [-1 -1 2 2]*a/10 + [a 0 0 0], 'Curvature', [1 1]);
drawnow

dispAngle = 0;
for i = 2:length(driverPitch)
    if floor(interpAngles(i)) > dispAngle || i == length(driverPitch)
        drawnow
        dispAngle = floor(interpAngles(i));
    end
    driverPitchNow = rotPolygon(driverPitch, [0 0], driverPitch(i,1), -driverPitch(i,2));
    followerPitchNow = rotPolygon(followerPitch, [0 0], -followerPitch(i,1), followerPitch(i,2));
    driverProfile = rotPolygon(driverProfile, [0 0], driverPitchNow(i-1,1), driverPitchNow(i-1,2));
    followerProfile = rotPolygon(followerProfile, [0 0], -followerPitchNow(i-1,1), -followerPitchNow(i-1,2));
    driverPitchClNow = rotPolygon(driverPitchCl, [0 0], driverPitch(i,1), -driverPitch(i,2));
    followerPitchClNow = rotPolygon(followerPitchCl, [0 0], -followerPitch(i,1), followerPitch(i,2));

    %pitchDist = norm(driverPitch(i,:) - driverPitch(i-1,:));
    %driverRack = driverRack + [norm(driverPitch(i,:))-norm(driverPitch(i-1,:)) -pitchDist];
    %followerRack = followerRack + [norm(driverPitch(i,:))-norm(driverPitch(i-1,:)) -pitchDist];
    driverRackNow = driverRack + [driverPitchNow(i,1) -pitchLengths(i)];
    followerRackNow = followerRack + [followerPitchNow(i,1) -pitchLengths(i)];
    temp = polyclip(driverProfile, driverRackNow, 'dif');
    driverProfile = [temp{1}{1} temp{2}{1}];
    temp = polyclip(followerProfile, followerRackNow, 'dif');
    followerProfile = [temp{1}{1} temp{2}{1}];

    set(dph, 'XData', driverProfile(:,1), 'YData', driverProfile(:,2));
    set(fph, 'XData', followerProfile(:,1)+a, 'YData', followerProfile(:,2));
    set(dpr, 'XData', closeArray(driverRackNow(:,1)), 'YData', closeArray(driverRackNow(:,2)));
    set(fpr, 'XData', closeArray(followerRackNow(:,1))+a, 'YData', closeArray(followerRackNow(:,2)));
    set(dpc, 'XData', closeArray(driverPitchClNow(:,1)), 'YData', closeArray(driverPitchClNow(:,2)));
    set(fpc, 'XData', closeArray(followerPitchClNow(:,1)+a), 'YData', closeArray(followerPitchClNow(:,2)));
end


function res = closeArray ( in )
    res = in([1:end 1]);
end


function res = rotPolygon ( polygon, pivot, cosrot, sinrot)
    res = (polygon - pivot) * [cosrot, sinrot; -sinrot, cosrot]/norm([sinrot cosrot]) + pivot;
end


function h = plotc ( x,y,varargin )
    %PLOTC .
    %   Auto close the line plotted
    %
    %   Parameters: x - X coords. (array)
    %               y - Y coords. (array)
    %               varargin - Rest of args. (any)
    %   Returns:    h - Handle returned by plot(). (handle)
    h = plot(closeArray(x),closeArray(y),varargin{:});
end

