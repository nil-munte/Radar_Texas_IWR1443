delete(instrfind);

configFileName = '1443config.cfg';
textFile = 'dataMatlab.txt';

[DATA_sphandle,UART_sphandle, config] = serialConfig(configFileName);
ConfigParameters = parseConfigFile(config);

fileID=fopen(textFile,'w');

stopVal = 0;
numDetection=1;

figure();
set(gcf, 'Position', get(0, 'Screensize'));
H = uicontrol('Style', 'PushButton', ...
                    'String', 'Stop', ...
                    'Callback', 'stopVal = 1','Position',[100 600 100 30]);
h = scatter([],[],'filled');
axis([-0.5,0.5,0,1.5]);
xlabel('X (m)'); ylabel('Y (m)');
grid minor

while(1)
    [dataOk,numDetection] = update(ConfigParameters,DATA_sphandle,fileID,h,numDetection);

    if stopVal == 1
        stopVal = 0;
        break;
    end

end

fprintf(UART_sphandle, 'sensorStop');
fclose(fileID);
delete(instrfind);
close all

function [DATA_sphandle,UART_sphandle, config] = serialConfig(configFileName)

    %UART_sphandle = serial('/dev/tty.usbmodemR10310411','BaudRate',115200);
    %UART_sphandle = serial('COM7','BaudRate',115200);
    UART_sphandle = serial('COM14','BaudRate',115200);
    set(UART_sphandle,'Parity','none')
    set(UART_sphandle,'Terminator','LF')
    fopen(UART_sphandle);
    
    %DATA_sphandle = serial('/dev/tty.usbmodemR10310414','BaudRate',921600);
    %DATA_sphandle = serial('COM8','BaudRate',921600);
    DATA_sphandle = serial('COM13','BaudRate',921600);
    set(DATA_sphandle,'Terminator', '');
    set(DATA_sphandle,'InputBufferSize', 65536);
    set(DATA_sphandle,'Timeout',10);
    set(DATA_sphandle,'ErrorFcn',@dispError);
    set(DATA_sphandle,'BytesAvailableFcnMode','byte');
    set(DATA_sphandle,'BytesAvailableFcnCount', 2^16+1);%BYTES_AVAILABLE_FCN_CNT);
    set(DATA_sphandle,'BytesAvailableFcn',@readUartCallbackFcn);
    fopen(DATA_sphandle);
    
    config = cell(1,100);
    fid = fopen(configFileName, 'r');
    if fid == -1
        fprintf('File %s not found!\n', configFileName);
        return;
    else
        fprintf('Opening configuration file %s ...\n', configFileName);
    end
    tline = fgetl(fid);
    k=1;
    while ischar(tline)
        config{k} = tline;
        tline = fgetl(fid);
        k = k + 1;
    end
    config = config(1:k-1);
    fclose(fid);
    
    mmwDemoCliPrompt = char('mmwDemo:/>');
    
    %Send CLI configuration to IWR14xx
    fprintf('Sending configuration from %s file to IWR16xx ...\n', configFileName);
    for k=1:length(config)
        command = config{k};
        fprintf(UART_sphandle, command);
        fprintf('%s\n', command);
        echo = fgetl(UART_sphandle); % Get an echo of a command
        done = fgetl(UART_sphandle); % Get "Done"
        prompt = fread(UART_sphandle, size(mmwDemoCliPrompt,2)); % Get the prompt back
    end

end

function ConfigParameters = parseConfigFile(config)

    for i = 1:length(config)
        configLine = strsplit(config{i});
        
        % We are only interested in the channelCfg, profileCfg and frameCfg parameters:
        if strcmp(configLine{1},'channelCfg')
            channelCfg.txChannelEn = str2double(configLine{3});
            
            numTxAzimAnt = bitand(bitshift(channelCfg.txChannelEn,0),1) +...
                    bitand(bitshift(channelCfg.txChannelEn,-2),1);
            numTxElevAnt = bitand(bitshift(channelCfg.txChannelEn,-1),1);
                    
            channelCfg.rxChannelEn = str2double(configLine{2});
            numRxAnt = bitand(bitshift(channelCfg.rxChannelEn,0),1) +...
                bitand(bitshift(channelCfg.rxChannelEn,-1),1) +...
                bitand(bitshift(channelCfg.rxChannelEn,-2),1) +...
                bitand(bitshift(channelCfg.rxChannelEn,-3),1);
            numTxAnt = numTxElevAnt + numTxAzimAnt;
            
        elseif  strcmp(configLine{1},'profileCfg')
            profileCfg.startFreq = str2double(configLine{3});
            profileCfg.idleTime =  str2double(configLine{4});
            profileCfg.rampEndTime = str2double(configLine{6});
            profileCfg.freqSlopeConst = str2double(configLine{9});
            profileCfg.numAdcSamples = str2double(configLine{11});
            profileCfg.numAdcSamplesRoundTo2 = 1;
            while profileCfg.numAdcSamples > profileCfg.numAdcSamplesRoundTo2
                profileCfg.numAdcSamplesRoundTo2 = profileCfg.numAdcSamplesRoundTo2 * 2;
            end 
            profileCfg.digOutSampleRate = str2double(configLine{12}); %uints: ksps
            
        elseif strcmp(configLine{1},'frameCfg')
            frameCfg.chirpStartIdx = str2double(configLine{2});
            frameCfg.chirpEndIdx = str2double(configLine{3});
            frameCfg.numLoops = str2double(configLine{4});
            frameCfg.numFrames = str2double(configLine{5});
            frameCfg.framePeriodicity = str2double(configLine{6});
            
        end
    end  
    
    
    ConfigParameters.numChirpsPerFrame = (frameCfg.chirpEndIdx -...
        frameCfg.chirpStartIdx + 1) *...
        frameCfg.numLoops;
    ConfigParameters.numDopplerBins = ConfigParameters.numChirpsPerFrame / numTxAnt;
    ConfigParameters.numRangeBins = profileCfg.numAdcSamplesRoundTo2;
    ConfigParameters.rangeResolutionMeters = 3e8 * profileCfg.digOutSampleRate * 1e3 /...
        (2 * profileCfg.freqSlopeConst * 1e12 * profileCfg.numAdcSamples);
    ConfigParameters.rangeIdxToMeters = 3e8 * profileCfg.digOutSampleRate * 1e3 /...
        (2 * profileCfg.freqSlopeConst * 1e12 * ConfigParameters.numRangeBins);
    ConfigParameters.dopplerResolutionMps = 3e8 / (2*profileCfg.startFreq*1e9 *...
        (profileCfg.idleTime + profileCfg.rampEndTime) *...
        1e-6 * ConfigParameters.numDopplerBins * numTxAnt);
    ConfigParameters.maxRange = 300 * 0.9 * profileCfg.digOutSampleRate /(2 * profileCfg.freqSlopeConst * 1e3);
    ConfigParameters.maxVelocity = 3e8 / (4*profileCfg.startFreq*1e9 *(profileCfg.idleTime + profileCfg.rampEndTime) * 1e-6 * numTxAnt);

end

function [dataOk,detObj]=readAndParseData14xx(configParameters,DATA_sphandle)
    
    OBJ_STRUCT_SIZE_BYTES = 12;
    BYTE_VEC_ACC_MAX_SIZE = 2^15;
    MMWDEMO_UART_MSG_DETECTED_POINTS = 1;
    MMWDEMO_UART_MSG_RANGE_PROFILE = 2;
    maxBufferSize = 2^15;
    
    magicWord = [2, 1, 4, 3, 6, 5, 8, 7];
    
    magicOk = 0;
    dataOk = 0;
    detObj = [];
    frameNumber = 0;
      
    persistent byteBuffer
    if isempty(byteBuffer)
        byteBuffer = zeros(maxBufferSize,1);
    end
    
    persistent byteBufferLength
    if isempty(byteBufferLength)
        byteBufferLength = 0;
    end
    
    persistent magiNotOkCounter
    if isempty(magiNotOkCounter)
        magiNotOkCounter = 0;
    end
    
    bytesToRead = get(DATA_sphandle,'BytesAvailable');
    if (bytesToRead ~= 0)
        % Read the Data Serial Port
        [bytevec, byteCount] = fread(DATA_sphandle, bytesToRead, 'uint8');
        
         % Check if the buffer is not full, and then add the data to the buffer:
        if(byteBufferLength + byteCount < maxBufferSize)
            byteBuffer(byteBufferLength+1:byteBufferLength + byteCount) = bytevec(1:byteCount);
            byteBufferLength = byteBufferLength + byteCount;
        end
        
    end
 
    % Check that the buffer is not empty:
    if byteBufferLength > 16
        byteBufferStr = char(byteBuffer);
        
        % Search for the magic number inside the buffer and check that at least one magic number has been found:
        startIdx = strfind(byteBufferStr', char(magicWord));
        if ~isempty(startIdx)
            
            % Check the position of the first magic number and put it at
            % the beginning of the buffer
            if length(startIdx) >= 2
                if startIdx(end-1) > 1
                    byteBuffer(1:byteBufferLength-(startIdx(1)-1)) = byteBuffer(startIdx(1):byteBufferLength);
                    byteBufferLength = byteBufferLength - (startIdx(1)-1);
                end
            else
                if startIdx(1) > 1
                    byteBuffer(1:byteBufferLength-(startIdx(1)-1)) = byteBuffer(startIdx(1):byteBufferLength);
                    byteBufferLength = byteBufferLength - (startIdx(1)-1);
                end
            end
            if byteBufferLength < 0
                byteBufferLength = 0;
            end
            
            totalPacketLen = sum(byteBuffer(8+4+[1:4]) .* [1 256 65536 16777216]');
            if ((byteBufferLength >= totalPacketLen) && (byteBufferLength ~= 0)) 
                magicOk = 1;
            else
                magicOk = 0;
            end
        end
    end
    
    if (magicOk == 1)
        %%%%% HEADER
        word = [1 256 65536 16777216];
        idx = 0;

        magicNumber = byteBuffer(idx + [1:8]);
        idx = idx + 8;
        Header.version = dec2hex(matmul(byteBuffer(idx+(1:4)),word));
        idx = idx + 4;
        Header.totalPacketLen = matmul(byteBuffer(idx+(1:4)),word);
        idx = idx + 4;
        Header.platform = dec2hex(matmul(byteBuffer(idx+(1:4)),word));
        idx = idx + 4;
        Header.frameNumber = matmul(byteBuffer(idx+(1:4)),word);
        frameNumber = Header.frameNumber;
        idx = idx + 4;
        Header.timeCpuCycles = matmul(byteBuffer(idx+(1:4)),word);
        idx = idx + 4;
        Header.numDetectedObj = matmul(byteBuffer(idx+(1:4)),word);
        idx = idx + 4;
        Header.numTLVs = matmul(byteBuffer(idx+(1:4)),word);
        idx = idx + 4;
        %Header.subFrameNumber = matmul(byteBuffer(idx+(1:4)),word);
        %idx = idx + 4;
        
        
        %%%%% TLV
        
        % Analyze each of TLV messages:
        for tlvIdx = 1:Header.numTLVs
            word = [1 256 65536 16777216];
            % First, analyze the TLV header (TLV type and length):
            tlv.type = matmul(byteBuffer(idx+(1:4)),word);
            idx = idx + 4;
            tlv.length = matmul(byteBuffer(idx+(1:4)),word);
            idx = idx + 4;
            
            % Check that the TLV message is of the right type (Detected objects):
           if tlv.type == MMWDEMO_UART_MSG_DETECTED_POINTS
                detObj =[];

                word = [1 2^8];
                tlv_numObj = matmul(byteBuffer(idx+(1:2)),word);
                idx = idx + 2;
                tlv_xyzQFormat = 2.^(matmul(byteBuffer(idx+(1:2)),word));
                idx = idx + 2;

                rangeIdx= zeros(1,tlv_numObj,'int16');
                dopplerIdx= zeros(1,tlv_numObj,'int16');
                peakVal= zeros(1,tlv_numObj,'int16');
                x= zeros(1,tlv_numObj,'int16');
                y= zeros(1,tlv_numObj,'int16');
                z= zeros(1,tlv_numObj,'int16');
                
                for objectNum=1:tlv_numObj
                    rangeIdx(objectNum) = matmul(byteBuffer(idx+(1:2)),word);
                    idx = idx + 2;
                    dopplerIdx(objectNum) = matmul(byteBuffer(idx+(1:2)),word);
                    idx = idx + 2;
                    peakVal(objectNum) = matmul(byteBuffer(idx+(1:2)),word);
                    idx = idx + 2;
                    x(objectNum)= matmul(byteBuffer(idx+(1:2)),word);
                    idx = idx + 2;
                    y(objectNum)= matmul(byteBuffer(idx+(1:2)),word);
                    idx = idx + 2;
                    z(objectNum)= matmul(byteBuffer(idx+(1:2)),word);
                    idx = idx + 2;

                end
                
                rangeVal = rangeIdx * configParameters.rangeIdxToMeters;
                dopplerIdx(dopplerIdx > (configParameters.numDopplerBins/2 - 1)) = dopplerIdx(dopplerIdx > (configParameters.numDopplerBins/2 - 1)) - 65535;
                dopplerVal = dopplerIdx * configParameters.dopplerResolutionMps;
                
                x = double(x) / tlv_xyzQFormat;
                y = double(y) / tlv_xyzQFormat;
                z = double(z) / tlv_xyzQFormat;
                
                detObj.numObj = tlv_numObj;
                detObj.x = x;
                detObj.y = y;
                detObj.z = z;
                
                dataOk=1;
                
            end
            
        end
        %Remove processed data
        if idx > 0
            shiftSize = Header.totalPacketLen;
            byteBuffer(1: byteBufferLength-shiftSize) = byteBuffer(shiftSize+1:byteBufferLength);
            byteBufferLength = byteBufferLength - shiftSize;
            if byteBufferLength < 0
                byteBufferLength = 0;
            end
        end
        
    else
        magiNotOkCounter = magiNotOkCounter + 1;
    end

end

function v = matmul(byteBuffer,word)
    byteBuffer=byteBuffer';
    v = sum(byteBuffer.*word);
    if (v>32767)
        v=v-65536;
    end
end


function [dataOk,numDetection] = update(ConfigParameters,DATA_sphandle,fileID,h,numDetection)

    [dataOk, detObj] = readAndParseData14xx(ConfigParameters,DATA_sphandle);
    if and(dataOk==1,size(detObj) > 0)
        x = -detObj.x;
        y = detObj.y;
        z = detObj.z;

        % Plot the radar points
        h.XData = x;
        h.YData = y;
        drawnow limitrate;

        for i=1:numel(x)
            fprintf(fileID,'%d %d %d %d\n',numDetection,x(i), y(i), z(i));
        end
        
        numDetection=numDetection+1;
        pause(0.05);

    end
end