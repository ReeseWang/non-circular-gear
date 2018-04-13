zSafe = 60.0;
roughToolNumber = '0606';
roughSpindleSpeed = 3000;
roughFeedRate = 100;
roughZOffsets = [9 6 3 0];
blankDia = 30;                                  % Blank material diameter
teethToolNumber = '0808';
teethToolDia = 3.175;
teethFeedRate = 600;
teethSpindleSpeed = 3000;
roughYCleariance = 5;
%roughXCoords = -4:0.8:5;
roughXCoords = -3;

teethXRange = [-2.55 -5.55];

gCodePrefix = 'G28 X0 Y0 Z0 CM0 M25\n';

if regen
    [ ...
        leftRoughingToolPath, ...
        leftRoughingToolPathExtra, ...
        leftTeethToolPath, ...
        rightRoughingToolPath, ...
        rightRoughingToolPathExtra, ...
        rightTeethToolPath ...
        ] = generation_v2('fp.txt', blankDia);
end

roughingCode('3001', leftRoughingToolPath, leftRoughingToolPathExtra, gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughZOffsets, roughXCoords, blankDia, roughYCleariance, roughFeedRate);
roughingCode('3003', rightRoughingToolPath, rightRoughingToolPathExtra, gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughZOffsets, roughXCoords, blankDia, roughYCleariance, roughFeedRate);
teethCode('3002', leftTeethToolPath, gCodePrefix, teethToolNumber, ...
    teethSpindleSpeed, zSafe, teethXRange, teethFeedRate, teethToolDia);
teethCode('3004', rightTeethToolPath, gCodePrefix, teethToolNumber, ...
    teethSpindleSpeed, zSafe, teethXRange, teethFeedRate, teethToolDia);

function roughingCode ( filename, toolPath, toolPathExtra, gCodePrefix, toolNumber, ...
        spindleSpeed, zSafe, zOffsets, xCoords, blankDia, yCleariance, feedRate)
    f = fopen(filename, 'w');
    gCodePrefix2 = '';
    fprintf(f, gCodePrefix);
    fprintf(f, 'T%s\n', toolNumber);
    fprintf(f, 'S%d M203\nM8\n', spindleSpeed);
    fprintf(f, gCodePrefix2);
    fprintf(f, 'G1 X%.3f F8000\n', zSafe);

    for i = 1:length(xCoords)
        for zOffset = zOffsets
            [YZ, A] = toMachineRef(toolPath{1});
            fprintf(f, 'G1 Z%.3f Y%.3f X%.3f %s F8000\n', ...
                xCoords(i), ...
                blankDia/2+yCleariance, ...
                zSafe + zOffset, ...
                ACoord(A, true));
            fprintf(f, 'G1 X%.3f F8000\n', YZ(2) + zOffset);
            fprintf(f, 'G93\n');                        % Inverse time feed rate mode
            fprintf(f, 'G1 Y%.3f F%.3f\n', YZ(1), feedRate/(YZ(1) + blankDia/2 + yCleariance));

            for j = 2:length(toolPath)
                [YZ, A] = toMachineRef(toolPath{j});
                fprintf(f, 'Y%.3f X%.3f %s F%.3f; %d\n', ...
                    YZ(1), YZ(2) + zOffset, ACoord(A), feedRate/equvDist(toolPath{j-1}, toolPath{j}, blankDia/2), j);
            end
            fprintf(f, 'G94\n');                        % Inverse time feed rate mode

            fprintf(f, 'G1 X%.3f F8000\n', zSafe + zOffset);

            for j = [1 3]
                [YZ, A] = toMachineRef(toolPathExtra(j:j+1,:));
                fprintf(f, 'G1 Y%.3f X%.3f %s F8000\n', ...
                    (blankDia/2 + yCleariance) * sign(YZ(1)), ...
                    zSafe + zOffset, ACoord(A, true));
                fprintf(f, 'X%.3f\n', YZ(2) + zOffset);
                fprintf(f, 'G1 Y%.3f F%.3f\n', ...
                    YZ(1), feedRate);
                fprintf(f, 'G1 X%.3f F8000\n', zSafe + zOffset);
            end
        end
    end
    fprintf(f, 'G94\nM5\nM9\nM2');
    fclose(f);
end

function teethCode ( filename, toolPath, gCodePrefix, toolNumber, ...
       spindleSpeed, zSafe, xRange, feedRate, toolDia )
    f = fopen(filename, 'w');
    fprintf(f, gCodePrefix);
    fprintf(f, 'T%s\n', toolNumber);
    fprintf(f, 'S%d M203\n', spindleSpeed);
    gCodePrefix2 = '';
    fprintf(f, gCodePrefix2);
    fprintf(f, 'G1 X%.3f F8000\n', zSafe);

    for idxPass = 1:length(toolPath)
        if idxPass == length(toolPath)
            xRange = xRange(end:-1:1);
        end

        for idxTooth = 1:length(toolPath{idxPass})
            [YZ, A] = toMachineRef(toolPath{idxPass}{idxTooth}{1});
            fprintf(f, 'G1 Z%.3f Y%.3f X%.3f %s F8000\n', ...
                xRange(1), ...
                YZ(1), ...
                zSafe, ...
                ACoord(A, true));
            fprintf(f, 'G1 X%.3f F8000\n', YZ(2) + 1);

            for idxCut = 1:length(toolPath{idxPass}{idxTooth})
                [YZ, A] = toMachineRef(toolPath{idxPass}{idxTooth}{idxCut});
                fprintf(f, 'G1 Z%.3f Y%.3f X%.3f %s F600\n', ...
                    xRange(1), YZ(1), YZ(2), ACoord(A));
                fprintf(f, 'Z%.3f F%.3f\n', xRange(2), feedRate);
                if idxCut ~= length(toolPath{idxPass}{idxTooth})
                    fprintf(f, 'G1 X%.3f F8000\nZ%.3f\n', YZ(2) + 0.2, xRange(1));
                else
                    fprintf(f, 'G1 X%.3f F8000\n', zSafe);
                end
            end
        end
    end
    fprintf(f, 'M5\nM9\nM2');
    fclose(f);
end


function [yz, adeg] = toMachineRef ( in )
    rotVector = in(2,:) - in(1,:);
    rotVector = rotVector/norm(rotVector);
    yz = in(1,:) * [rotVector(2), 2*rotVector(1); -rotVector(1), 2*rotVector(2)];
    adeg = (pi + atan2(-rotVector(2), -rotVector(1))) * 180/pi; % 0~360 degree
end

function [yz, adeg] = toWorkRef ( in )
    rotVector = in(2,:) - in(1,:);
    rotVector = rotVector/norm(rotVector);
    yz = in(1,:);
    adeg = (pi + atan2(-rotVector(2), -rotVector(1))) * 180/pi; % 0~360 degree
end

function dist = equvDist( toolRefVectorPrev, toolRefVectorNow, maxRadius)
    dist = norm(toolRefVectorNow(1,:) - toolRefVectorPrev(1,:));
    vectPrev = toolRefVectorPrev(2,:) - toolRefVectorPrev(1,:);
    vectPrev = vectPrev/norm(vectPrev);
    vectNow = toolRefVectorNow(2,:) - toolRefVectorNow(1,:);
    vectNow = vectNow/norm(vectNow);
    ang = atan2(dot(vectNow, vectPrev * [0 1; -1 0]), dot(vectNow, vectPrev));
    dist = dist + maxRadius * ang;
end

function str = ACoord ( ang, reset )
    persistent prevAng
    if ~exist('reset', 'var')
        reset = false;
    end
    if isempty(prevAng) || reset || round(ang, 3) ~= round(prevAng, 3)
        str = sprintf('CM%.3f', ang);
    else
        str = '';
    end
    prevAng = ang;
end

