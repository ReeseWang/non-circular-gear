zSafe = 50.0;
roughToolNumber = 1;
roughSpindleSpeed = 5000;
roughFeedRate = 3000;
blankDia = 30;                                  % Blank material diameter
teethToolNumber = 10;
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
        rightRoughToolPath, ...
        rightRoughToolPathExtra, ...
        rightTeethToolPath ...
        ] = generation_v2('fp.txt', blankDia);
end

roughingCode('leftRough.ngc', leftRoughingToolPath, leftRoughingToolPathExtra, gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughXCoords, blankDia, roughYCleariance, roughFeedRate);

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
            ACoord(A));
        fprintf(f, 'G0 Z%.3f\n', YZ(2));
        fprintf(f, 'G1 Y%.3f F%.3f\n', YZ(1), feedRate/(YZ(1) + blankDia/2 + yCleariance));

        for j = 2:length(toolPath)
            [YZ, A] = toMachineRef(toolPath{j});
            fprintf(f, 'Y%.3f Z%.3f %s F%.3f; %d\n', ...
                YZ(1), YZ(2), ACoord(A), feedRate/equvDist(toolPath{j-1}, toolPath{j}, blankDia/2), j);
        end

        fprintf(f, 'G0 Z%.3f\n', zSafe);
    end
    fprintf(f, 'G94\nM5\nM2');
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

function str = ACoord ( ang )
    persistent prevAng
    if isempty(prevAng) || round(ang, 3) ~= round(prevAng, 3)
        str = sprintf('A%.3f', ang);
    else
        str = '';
    end
    prevAng = ang;
end

