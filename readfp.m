%% 导入文本文件中的数据。
% 用于从以下文本文件导入数据的脚本:
%
%    C:\Users\ruoxi\Documents\MATLAB\non-cirular-gear\fp.txt
%
% 要将代码扩展到其他选定数据或其他文本文件，请生成函数来代替脚本。

% 由 MATLAB 自动生成于 2017/11/14 11:08:06

%% 初始化变量。
filename = 'fp.txt';
delimiter = ',';

%% 每个文本行的格式:
%   列1: 双精度值 (%f)
%	列2: 双精度值 (%f)
% 有关详细信息，请参阅 TEXTSCAN 文档。
formatSpec = '%f%f%[^\n\r]';

%% 打开文本文件。
fileID = fopen(filename,'r');

%% 根据格式读取数据列。
% 该调用基于生成此代码所用的文件的结构。如果其他文件出现错误，请尝试通过导入工具重新生成代码。
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

%% 关闭文本文件。
fclose(fileID);

%% 对无法导入的数据进行的后处理。
% 在导入过程中未应用无法导入的数据的规则，因此不包括后处理代码。要生成适用于无法导入的数据的代码，请在文件中选择无法导入的元胞，然后重新生成脚本。

%% 创建输出变量
polData = [dataArray{1:end-1}];
%% 清除临时变量
clearvars filename delimiter formatSpec fileID dataArray ans;
