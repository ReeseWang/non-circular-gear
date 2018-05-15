zSafe = 60.0;
roughToolNumber = '0404';
roughSpindleSpeed = 5000;
roughFeedRate = 600;
roughZOffsets = [10 5 0];
blankDia = 25;                                  % Blank material diameter
teethToolNumber = '0606';
teethToolDia = 3.175; % Useless for now
teethFeedRate = 2000;
teethSpindleSpeed = 5000;
roughYCleariance = 5;
%roughXCoords = -4:0.8:5;
roughXCoords = -1;

teethXRange = roughXCoords + 1.1 * [1 -1];

gCodePrefix = 'M05\nM25\nCM 0.0\nG98\n';

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

roughingCode('O4591.nc', leftRoughingToolPath, leftRoughingToolPathExtra, "%%\nO4591\n" + gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughZOffsets, roughXCoords, blankDia, roughYCleariance, roughFeedRate);
roughingCode('O4592.nc', rightRoughingToolPath, rightRoughingToolPathExtra, "%%\nO4592\n" + gCodePrefix, ...
    roughToolNumber, roughSpindleSpeed, zSafe, roughZOffsets, roughXCoords, blankDia, roughYCleariance, roughFeedRate);
teethCode('O4593.nc', leftTeethToolPath, "%%\nO4593\n" + gCodePrefix, teethToolNumber, ...
    teethSpindleSpeed, zSafe, teethXRange, teethFeedRate, teethToolDia);
teethCode('O4594.nc', rightTeethToolPath, "%%\nO4594\n" + gCodePrefix, teethToolNumber, ...
    teethSpindleSpeed, zSafe, teethXRange, teethFeedRate, teethToolDia);

function roughingCode ( filename, toolPath, toolPathExtra, gCodePrefix, toolNumber, ...
        spindleSpeed, zSafe, zOffsets, xCoords, blankDia, yCleariance, feedRate)
    f = fopen(filename, 'w');
    gCodePrefix2 = '';
    fprintf(f, gCodePrefix);
    fprintf(f, 'T%s\n', toolNumber);
    fprintf(f, 'S%d M203\nM8\n', spindleSpeed);
    fprintf(f, gCodePrefix2);
    fprintf(f, 'G0 X%.3f\n', zSafe);

    for i = 1:length(xCoords)
        for zOffset = zOffsets
            [YZ, A] = toMachineRef(toolPath{1});
            fprintf(f, 'G0 Z%.3f Y%.3f X%.3f %s\n', ...
                xCoords(i), ...
                blankDia/2+yCleariance, ...
                zSafe + zOffset, ...
                ACoord(A, true));
            fprintf(f, 'G0 X%.3f\n', YZ(2) + zOffset);
            %fprintf(f, 'G93\n');                        % Inverse time feed rate mode
            fprintf(f, 'G1 Y%.3f F%.3f\n', YZ(1), feedRate);%/(YZ(1) + blankDia/2 + yCleariance));

            for j = 2:length(toolPath)
                [YZ, A] = toMachineRef(toolPath{j});
                fprintf(f, 'Y%.3f X%.3f %s\n', ...
                    YZ(1), YZ(2) + zOffset, ACoord(A));%/equvDist(toolPath{j-1}, toolPath{j}, blankDia/2));
            end
            %fprintf(f, 'G94\n');                        % Inverse time feed rate mode

            fprintf(f, 'G0 X%.3f\n', zSafe + zOffset);

            for j = [1 3]
                [YZ, A] = toMachineRef(toolPathExtra(j:j+1,:));
                fprintf(f, 'G0 Y%.3f X%.3f %s\n', ...
                    (blankDia/2 + yCleariance) * sign(YZ(1)), ...
                    zSafe + zOffset, ACoord(A, true));
                fprintf(f, 'X%.3f\n', YZ(2) + zOffset);
                fprintf(f, 'G1 Y%.3f\n', ...
                    YZ(1));
                fprintf(f, 'G0 X%.3f\n', zSafe + zOffset);
            end
        end
    end
    fprintf(f, 'G28 U0 V0 W0\nM26\nM205\nM9\nM99\n%%');
    fclose(f);
end

function teethCode ( filename, toolPath, gCodePrefix, toolNumber, ...
       spindleSpeed, zSafe, xRange, feedRate, toolDia )
    f = fopen(filename, 'w');
    fprintf(f, gCodePrefix);
    %fprintf(f, 'T%s\n', toolNumber);
    fprintf(f, 'S%d M203\n', spindleSpeed);
    gCodePrefix2 = '';
    fprintf(f, gCodePrefix2);
    fprintf(f, 'G0 X%.3f\n', zSafe);

    for idxPass = 1:length(toolPath)
        if idxPass == length(toolPath)
            xRange = xRange(end:-1:1);
        end

        for idxTooth = 1:length(toolPath{idxPass})
            [YZ, A] = toMachineRef(toolPath{idxPass}{idxTooth}{1});
            fprintf(f, 'G0 Z%.3f Y%.3f X%.3f %s \n', ...
                xRange(1), ...
                YZ(1), ...
                zSafe, ...
                ACoord(A, true));
            fprintf(f, 'G0 X%.3f\n', YZ(2) + 1);
            fprintf(f, 'F%.3f\n', feedRate);

            for idxCut = 1:length(toolPath{idxPass}{idxTooth})
                [YZ, A] = toMachineRef(toolPath{idxPass}{idxTooth}{idxCut});
                fprintf(f, 'G1 Z%.3f Y%.3f X%.3f %s\n', ...
                    xRange(1), YZ(1), YZ(2), ACoord(A));
                fprintf(f, 'Z%.3f\n', xRange(2));
                if idxCut ~= length(toolPath{idxPass}{idxTooth})
                    fprintf(f, 'G0 X%.3f\nZ%.3f\n', YZ(2) + 0.2, xRange(1));
                else
                    fprintf(f, 'G0 X%.3f\n', zSafe);
                end
            end
        end
    end
    fprintf(f, 'G28 U0exi V0 W0\nM26\nM205\nM9\nM99\n%%');
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

