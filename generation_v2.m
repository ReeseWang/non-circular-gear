function generation_v2 ()
    addpath(genpath('./'))
    a = 15;                                     % Center distance
    pAngle = 30*pi/180;                         % Pressure Angle
    module = 1;                                 % Module
    pitch = pi * module;                        % Curve pitch
    addDist = 0.9 * module;                     % Addendum distance
    dedDist = 1.1 * module;                     % Dedendum distance
    posLimiterLeng = 3;                         % Position limiter length
    errTol = 1e-4;                              % Error checking tolerance
    angTol = 5*pi/180;                          % Angular tolerance, rad
    toolTipRadius = 0.2;                        % Forming tool radius
    toolFullAngle = 45*pi/180;
    toolFluteLength = 5;
    toolDiameter = 3.175;
    toolStickOut = 15;
    cutSteps = 2;                               % How many steps to simulate a rack using tip tool
    cutDepth = 0.2;                             % Z step
    rightIsDriver = true;                           % Reverse the role of driver and follower (pos don't change)
    machineRef = true;                          % Demo in machine reference frame to avoid audience confusion.
    leftRotateMargin = 1*pi/180;                % Rotate margin of driver
    blankDia = 30;                              % Blank material diameter
    roughToolDia = 10;
    roughToolFluteLength = 50;

    %% Read input, interpolate, and generate pitch curves.
    filename = 'fp.txt';
    delimiter = ',';
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    polData = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;

    polData = [polData(:,1)/180*pi polData(:,2)];   % Convert to rad
    polData = polData([1 1:end end],:);
    polData(1,:) = polData(2,:) - leftRotateMargin/(polData(3,1) - polData(2,1))*(polData(3,:) - polData(2,:));
    polData(end,:) = polData(end-1,:) - leftRotateMargin/(polData(end-1,1) - polData(end-2,1))*(polData(end-1,:) - polData(end-2,:));
    dfStruct = spline(polData(:,1), polData(:,2));   % driverAngle = f(followerAngle)
    ddfStruct = fnder(dfStruct);                % 2nd order derivative
    fStruct = fnint(dfStruct);                  % Original form
    offset = -ppval(fStruct, 0);
    fFunc = @(theta) ppval(fStruct, theta) + offset; % Ensure f(0) = 0
    dfFunc = @(theta) ppval(dfStruct, theta);
    ddfFunc = @(theta) ppval(ddfStruct, theta);
    leftRadiusFunc = @(theta) a * (1 - 1./(1+dfFunc(theta))); % r_d = a * f'(theta) / (1 + f'(theta))

    dDRFunc = @(theta) a * ddfFunc(theta)./(1 + dfFunc(theta)).^2; % r_d' = a * f''(theta) / (1 + f'(theta))^2
    dsFunc = @(theta) sqrt(leftRadiusFunc(theta).^2 + dDRFunc(theta).^2); % ds = sqrt(r^2 + (df/dtheta)^2) dtheta
    tangentFunc = @(theta) ddfFunc(theta)./(dfFunc(theta).*(1+dfFunc(theta))); % tan(alpha-theta)=f''[theta]/(f'[theta](1+f'[theta]))

    leftPolarAngles = polData(1,1):angTol:polData(end,1);
    if leftPolarAngles(end) < polData(end,1)    % Last value should reach the max in domain
        leftPolarAngles = [leftPolarAngles polData(end,1)];
    end
    leftPolarAngles = leftPolarAngles';

    rightPolarAngles = fFunc(leftPolarAngles);
    rightPitchPolarRadius = a * 1./(1+dfFunc(leftPolarAngles));
    leftPitchPolarRadius = a - rightPitchPolarRadius;
    leftTangentAngles = tangentFunc(leftPolarAngles);
    %rightTangentAngles = -tangentFunc(leftPolarAngles); % Not necessary

    %% XY coords of driver follower pitch line, extend by tangent
    [leftPitch(:,1), leftPitch(:,2)] = pol2cart(leftPolarAngles, leftPitchPolarRadius);
    [rightPitch(:,1), rightPitch(:,2)] = pol2cart(rightPolarAngles, rightPitchPolarRadius);
    tanAngle = leftPolarAngles(end) - leftTangentAngles(end);
    leftPitch(end+1,:) = ...
        leftPitch(end,:) + posLimiterLeng * ...
        [-sin(tanAngle), cos(tanAngle)];
    tanAngle = rightPolarAngles(end) + leftTangentAngles(end);
    rightPitch(end+1,:) = ...
        rightPitch(end,:) + posLimiterLeng * ...
        [-sin(tanAngle), cos(tanAngle)];
    tanAngle = leftPolarAngles(1) - leftTangentAngles(1);
    leftPitch = [ ...
        leftPitch(1,:) - posLimiterLeng * [-sin(tanAngle) cos(tanAngle)]; leftPitch];
    tanAngle = rightPolarAngles(1) + leftTangentAngles(1);
    rightPitch = [ ...
        rightPitch(1,:) - posLimiterLeng * [-sin(tanAngle) cos(tanAngle)]; rightPitch];

    %plot(leftPitch(:,1), leftPitch(:,2), - rightPitch(:,1) + a, rightPitch(:,2));
    %axis equal

    %% Integrate pitch arc lengths.
    pitchArcLengths = zeros(size(leftPolarAngles));
    for i=2:length(leftPolarAngles)
        pitchArcLengths(i) = pitchArcLengths(i-1) + ...
            integral(dsFunc, leftPolarAngles(i-1), leftPolarAngles(i), ...
            'RelTol', 0, ...
            'AbsTol', 1e-15);
    end
    temp = leftPolarAngles;
    temp(temp <= 0) = Inf;
    [closestToZero, idx] = min(temp);
    offset = integral(dsFunc, 0, closestToZero, 'RelTol', 0, 'AbsTol', 1e-15) - ...
        pitchArcLengths(idx);
    pitchArcLengths = pitchArcLengths + offset;

    %% Expand pitch curve by addDist to obtain initial profiles.
    temp = polyout(leftPitch(:,1), leftPitch(:,2), addDist, 'm');
    leftProfileTarget = [temp{1}{1} temp{2}{1}];
    temp = polyout(rightPitch(:,1), rightPitch(:,2), addDist, 'm');
    rightProfileTarget = [temp{1}{1} temp{2}{1}];

    %fill(leftProfile(:,1), leftProfile(:,2), 'y', ...
    %    -rightProfile(:,1)+a, rightProfile(:,2), 'r');
    %axis equal

    %% Generate zigzag shape of rack to simulate
    zigZagHalfHight = pitch/4/tan(pAngle);
    rackZigZag = ((-pitch/4:pitch/2:pitchArcLengths(end)+pitch/2)'*[1 1]);
    for i = 1:size(rackZigZag, 1)
        rackZigZag(i,1) = zigZagHalfHight * (-1)^i;
    end
    if rackZigZag(end,1) < 0                    % Driver should drive the follower at the end
        rackZigZag = [rackZigZag; zigZagHalfHight rackZigZag(end,2) + pitch/2];
    end
    if rightIsDriver
        rackZigZag(:,1) = - rackZigZag(:,1);
    end

    %% Check if tool radius satisfy the requirement.
    maxToolRadius = (dedDist - addDist) / (1 - sin(pAngle));
    if toolTipRadius > maxToolRadius
        error(['Tool fillet radius may lead to interference. '...
            'Its radius should be lower than %.3f'], maxToolRadius)
    end
    maxToolRadius = (zigZagHalfHight - dedDist) * tan(pAngle) / tan(pi/4 - pAngle/2);
    if toolTipRadius > maxToolRadius
        warning(['Tool fillet wiped top land out, decrease addendum hight or decrease tool radius. ' ...
        'Maximum fillet radius is %.3f'], maxToolRadius)
    end

    %% Generate tool shape
    filletAngles = (0:angTol:(pi - toolFullAngle)/2)';
    if filletAngles(end) < (pi - toolFullAngle)/2
        filletAngles(end+1) = (pi - toolFullAngle)/2;
    end
    toolFull = [ ...
        -toolDiameter/2, toolStickOut;
        -toolDiameter/2, toolDiameter/2/tan(toolFullAngle/2) - toolTipRadius/sin(toolFullAngle/2);
        toolTipRadius * ...
            [-sin(filletAngles(end:-1:2)), -cos(filletAngles(end:-1:2));
            sin(filletAngles), -cos(filletAngles)];
        toolDiameter/2, toolDiameter/2/tan(toolFullAngle/2) - toolTipRadius/sin(toolFullAngle/2);
        toolDiameter/2, toolStickOut
    ];

    %plot(toolFull(:,1), toolFull(:,2))
    %axis equal

    toolCutMask = [ ...
        -toolDiameter/2-1, toolFluteLength;
        -toolDiameter/2-1, -toolTipRadius-1;
        toolDiameter/2+1, -toolTipRadius-1;
        toolDiameter/2+1, toolFluteLength
    ];
    temp = polyclip(toolFull, toolCutMask, 'int');
    toolCut = [temp{1}{1} temp{2}{1}];
    temp = polyclip(toolFull, toolCutMask, 'dif');
    toolNonCut = [temp{1}{1} temp{2}{1}];

    %fill(toolCut(:,1), toolCut(:,2), 'g', toolNonCut(:,1), toolNonCut(:,2), 'r');
    %axis equal

    distanceRangeHalf = (zigZagHalfHight - dedDist) * tan(pAngle) - toolTipRadius * ...
        tan(pi/4 - pAngle/2);
    cutOffsets = [ ...
        linspace(-pAngle+toolFullAngle/2, pAngle-toolFullAngle/2, cutSteps);
        linspace(-distanceRangeHalf, distanceRangeHalf, cutSteps)
    ]';

    %% Cutting the gear on the left

    f = figure('OuterPosition', [0 0 1200 800]);
    ax = axes(f);
    set(ax, 'XLim', [-a 2*a], 'YLim', [-a a], 'YLimMode', 'manual', 'DataAspectRatio', [1 1 1])
    hold all
    lph = fill(leftProfile(:,1), leftProfile(:,2), 'y', 'Visible', 'off');
    rph = fill(-rightProfile(:,1)+a, leftProfile(:,2), 'g', 'Visible', 'off');
    inh = fill([0], [0], 'b');                  % Intersection, materials being cut
    zh = plot(rackZigZag(:,1), rackZigZag(:,2), '-.');                        % Zigzag, rack to simulate
    tch = plot(toolCut(:,1), toolCut(:,2));     % Cutting portion of the tool
    trh = plot([0 0], [0 1]);                   % Tool reference vector
    lc = rectangle('Position', [-1 -1 2 2]*a/10, 'Curvature', [1 1], 'Visible', 'off');
    rc = rectangle('Position', [-1 -1 2 2]*a/10 + [a 0 0 0], 'Curvature', [1 1], 'Visible', 'off');
    pause(1)

    show = @(x) set(x, 'Visible', 'on');
    hide = @(x) set(x, 'Visible', 'off');

    show([lph lc]);
    isRightGear = false;
    cutTeeth();
    hide([lph lc]);

    show([rph rc]);
    isRightGear = true;
    cutTeeth();

    show([lph lc]);
    hide([inh tch trh zh]);
    motsim();

    function motsim ()
        for i = 1:size(leftPolarAngles, 1)
            isRightGear = false;
            replot(lph, ...
                rotPolygon(leftProfile, ...
                    cos(leftPolarAngles(i)), ...
                    -sin(leftPolarAngles(i))));
            isRightGear = true;
            replot(rph, ...
                rotPolygon(rightProfile, ...
                    cos(rightPolarAngles(i)), ...
                    -sin(rightPolarAngles(i))));
            drawnow
        end
    end

    function cutTeeth()
        if isRightGear
            polarAngles = rightPolarAngles;
            tangentAngles = -leftTangentAngles;
            pitchPolarRadius = rightPitchPolarRadius;
            iRange = find(rackZigZag(:,1) > 0)';
        else
            polarAngles = leftPolarAngles;
            tangentAngles = leftTangentAngles;
            pitchPolarRadius = leftPitchPolarRadius;
            iRange = find(rackZigZag(:,1) < 0)';
        end
        for offsIdx = 1:size(cutOffsets, 1)     % For every cut offset
            toolCutThisPass = rotPolygon(toolCut, sin(cutOffsets(offsIdx,1)), ... % Rotate 90 degrees cw additionally
                -cos(cutOffsets(offsIdx,1)), [0 cutOffsets(offsIdx,2)]);
            toolNonCutThisPass = rotPolygon(toolNonCut, sin(cutOffsets(offsIdx,1)), ...
                -cos(cutOffsets(offsIdx,1)), [0 cutOffsets(offsIdx,2)]);
            toolRefVectorThisPass = rotPolygon([0 0; 0 1], sin(cutOffsets(offsIdx,1)), ...
                -cos(cutOffsets(offsIdx,1)), [0 cutOffsets(offsIdx,2)]);
            for i = iRange  % For every tooth
                toolCutThisTooth = toolCutThisPass + [-dedDist+toolTipRadius rackZigZag(i,2)];
                toolRefVectorThisTooth = toolRefVectorThisPass + [-dedDist+toolTipRadius rackZigZag(i,2)];
                for j = 1:size(polarAngles, 1)
                    cosrot = cos(polarAngles(j) - tangentAngles(j));
                    sinrot = sin(polarAngles(j) - tangentAngles(j));
                    anchor = [0 pitchArcLengths(j)];
                    move = pitchPolarRadius(j) * ...
                       [cos(polarAngles(j)) sin(polarAngles(j))] - anchor;
                    fun = @(x) rotPolygon(x, cosrot, sinrot, move, anchor);
                    toolCutNow = fun(toolCutThisTooth);
                    toolRefVectorNow = fun(toolRefVectorThisTooth);
                    if isRightGear
                        temp = polyclip(rightProfile, toolCutNow, 'int');
                        zigZagNow = fun([-rackZigZag(:,1), rackZigZag(:,2)]);
                    else
                        temp = polyclip(leftProfile, toolCutNow, 'int');
                        zigZagNow = fun(rackZigZag);
                    end
                    if ~isempty(temp{1})        % If has intersection
                        toolCutNow = {toolCutNow};
                        toolRefVectorNow = {toolRefVectorNow};
                        if j == 1               % First cut of the tooth, might need Z step cut
                            fun = @(x) x{end} + ...
                                cutDepth * normr(toolRefVectorNow{1}(2,:) - toolRefVectorNow{1}(1,:));
                            while true
                                toolCutNow{end+1} = fun(toolCutNow);
                                if isRightGear
                                    temp = polyclip(rightProfile, toolCutNow{end}, 'int');
                                else
                                    temp = polyclip(leftProfile, toolCutNow{end}, 'int');
                                end
                                if isempty(temp{1})
                                    toolCutNow(end) = [];
                                    break
                                else
                                    toolRefVectorNow{end+1} = fun(toolRefVectorNow);
                                end
                            end
                        end
                        while true
                            if isempty(toolCutNow)
                                assert(isempty(toolRefVectorNow))
                                break
                            end
                            if isRightGear
                                temp = polyclip(rightProfile, toolCutNow{end}, 'int');
                            else
                                temp = polyclip(leftProfile, toolCutNow{end}, 'int');
                            end
                            intersection = [temp{1}{1}, temp{2}{1}];
                            if isRightGear
                                temp = polyclip(rightProfile, toolCutNow{end}, 'dif');
                                rightProfile = [temp{1}{1} temp{2}{1}];
                            else
                                temp = polyclip(leftProfile, toolCutNow{end}, 'dif');
                                leftProfile = [temp{1}{1} temp{2}{1}];
                            end
                            if machineRef
                                fun = @(x) rotPolygon(x, ...
                                    toolRefVectorNow{end}(2,2) - toolRefVectorNow{end}(1,2), ...
                                    toolRefVectorNow{end}(2,1) - toolRefVectorNow{end}(1,1));
                            else
                                fun = @(x) x;
                            end
                            replot(inh, fun(intersection));
                            if isRightGear
                                replot(rph, fun(rightProfile));
                            else
                                replot(lph, fun(leftProfile));
                            end
                            replot(zh, fun(zigZagNow));
                            replot(tch, fun(toolCutNow{end}));
                            replot(trh, fun(toolRefVectorNow{end}));
                            drawnow
                            toolCutNow(end) = [];
                            toolRefVectorNow(end) = [];
                        end
                    end
                end
            end
        end
    end


    function replot(h, x, y)
        if isRightGear
            x(:,1) = -x(:,1) + a;
        end
        if ~exist('y', 'var')
            set(h, 'XData', x(:,1), 'YData', x(:,2));
        else
            set(h, 'XData', x, 'YData', y);
        end
    end

    function res = rotPolygon ( polygon, cosrot, sinrot, move, pivot)
        if ~exist('pivot', 'var')
            pivot = [0 0];
        end
        if ~exist('move', 'var')
            move = [0 0];
        end
        res = (polygon - pivot) * [cosrot, sinrot; -sinrot, cosrot]/norm([sinrot cosrot]) + ...
            pivot + move;
    end
end
