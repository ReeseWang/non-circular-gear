zSafe = 50.0;
roughToolNumber = 1;
roughSpindleSpeed = 5000;
roughFeedRate = 600;
blankDia = 30;                                  % Blank material diameter
teethToolNumber = 3;
teethToolDia = 3.175;
teethFeedRate = 100;
teethSpindleSpeed = 5000;
roughYCleariance = 5;
roughXCoords = -4:0.8:5;

teethXRange = [-0.2 5.2];

gCodePrefix = 'G17 G40 G21 G90 G93 G54 G49 G97 G64 P0.01 Q0.001 G91.1\n';

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

roughingCode('leftRough.ngc', leftRoughingToolPath, leftRoughingToolPathExtra, gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughXCoords, blankDia, roughYCleariance, roughFeedRate);
roughingCode('rightRough.ngc', rightRoughingToolPath, rightRoughingToolPathExtra, gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughXCoords, blankDia, roughYCleariance, roughFeedRate);
teethCode('leftTeeth.ngc', leftTeethToolPath, gCodePrefix, teethToolNumber, ...
    teethSpindleSpeed, zSafe, teethXRange, teethFeedRate, teethToolDia);

function roughingCode ( filename, toolPath, toolPathExtra, gCodePrefix, toolNumber, ...
        spindleSpeed, zSafe, xCoords, blankDia, yCleariance, feedRate)
    f = fopen(filename, 'w');
    fprintf(f, gCodePrefix);
    fprintf(f, 'M6 T%d\nG43 H%d\n', toolNumber, toolNumber);
    fprintf(f, 'M3 S%d\n', spindleSpeed);
    fprintf(f, 'G0 Z%.3f\n', zSafe);
    fprintf(f, 'G93\n');                        % Inverse time feed rate mode

    for i = 1:length(xCoords)
        [YZ, A] = toMachineRef(toolPath{1});
        fprintf(f, 'G0 X%.3f Y%.3f Z%.3f %s\n', ...
            xCoords(i), ...
            -blankDia/2-yCleariance, ...
            zSafe, ...
            ACoord(A, true));
        fprintf(f, 'G0 Z%.3f\n', YZ(2));
        fprintf(f, 'G1 Y%.3f F%.3f\n', YZ(1), feedRate/(YZ(1) + blankDia/2 + yCleariance));

        for j = 2:length(toolPath)
            [YZ, A] = toMachineRef(toolPath{j});
            fprintf(f, 'Y%.3f Z%.3f %s F%.3f; %d\n', ...
                YZ(1), YZ(2), ACoord(A), feedRate/equvDist(toolPath{j-1}, toolPath{j}, blankDia/2), j);
        end

        fprintf(f, 'G0 Z%.3f\n', zSafe);

        for j = [1 3]
            [YZ, A] = toMachineRef(toolPathExtra(j:j+1,:));
            fprintf(f, 'G0 Y%.3f Z%.3f %s\n', ...
                (blankDia/2 + yCleariance) * sign(YZ(1)), ...
                zSafe, ACoord(A, true));
            fprintf(f, 'Z%.3f\n', YZ(2));
            fprintf(f, 'G1 Y%.3f F%.3f\n', ...
                YZ(1), feedRate/abs(YZ(1) - blankDia/2 - yCleariance));
            fprintf(f, 'G0 Z%.3f\n', zSafe);
        end
    end
    fprintf(f, 'G94\nM5\nM2');
    fclose(f);
end

function teethCode ( filename, toolPath, gCodePrefix, toolNumber, ...
       spindleSpeed, zSafe, xRange, feedRate, toolDia )
    f = fopen(filename, 'w');
    fprintf(f, gCodePrefix);
    fprintf(f, 'M6 T%d\nG43 H%d\n', toolNumber, toolNumber);
    fprintf(f, 'M3 S%d\n', spindleSpeed);
    fprintf(f, 'G0 Z%.3f\n', zSafe);

    for idxPass = 1:length(toolPath)
        if idxPass == length(toolPath)
            xRange = xRange(end:-1:1);
        end

        for idxTooth = 1:length(toolPath{idxPass})
            [YZ, A] = toMachineRef(toolPath{idxPass}{idxTooth}{1});
            fprintf(f, 'G0 X%.3f Y%.3f Z%.3f %s\n', ...
                xRange(1), ...
                YZ(1), ...
                zSafe, ...
                ACoord(A, true));
            fprintf(f, 'G0 Z%.3f\n', YZ(2) + 1);

            for idxCut = 1:length(toolPath{idxPass}{idxTooth})
                [YZ, A] = toMachineRef(toolPath{idxPass}{idxTooth}{idxCut});
                fprintf(f, 'G1 X%.3f Y%.3f Z%.3f %s F60\n', ...
                    xRange(1), YZ(1), YZ(2), ACoord(A));
                fprintf(f, 'X%.3f F%.3f\n', xRange(2), feedRate);
                if idxCut ~= length(toolPath{idxPass}{idxTooth})
                    fprintf(f, 'G0 Z%.3f\nX%.3f\n', YZ(2) + 0.2, xRange(1));
                else
                    fprintf(f, 'G0 Z%.3f\n', zSafe);
                end
            end
        end
    end
    fprintf(f, 'M5\nM2');
    fclose(f);
end


function [yz, adeg] = toMachineRef ( in )
    rotVector = in(2,:) - in(1,:);
    rotVector = rotVector/norm(rotVector);
    yz = in(1,:) * [-rotVector(2), rotVector(1); rotVector(1), rotVector(2)];
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
        str = sprintf('A%.3f', ang);
    else
        str = '';
    end
    prevAng = ang;
end

