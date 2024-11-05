function [trackResults, channel]= tracking_MMT(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the PRN code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

% Remain code and carrier phase
trackResults.remCodePhase       = inf(1, settings.msToProcess);
trackResults.remCarrPhase       = inf(1, settings.msToProcess);

% Multipath Mitigation Technology variables
% Added by Sergio Vicenzo - 5 Nov 2024
trackResults.MMT_codeDelay_LOS  = nan(1, settings.msToProcess);
trackResults.MMT_codeDelay_NLOS = nan(1, settings.msToProcess);
trackResults.MMT_codeDelay_NLOS_relative  = nan(1, settings.msToProcess);
trackResults.MMT_time  = nan(1, settings.msToProcess);

%C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));



%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = settings.intTime;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
    settings.dllDampingRatio, ...
    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = settings.intTime;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);

% Start waitbar
hwb = waitbar(0,'Tracking...');

%Adjust the size of the waitbar to insert text
CNoPos=get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels

    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample)

        % Edited by Sergio Vicenzo to also handle "int16" files
        % 12 Feb 2024
        if strcmp(settings.dataType,'int16') 
            fseek(fid, ...
                dataAdaptCoeff*(settings.skipNumberOfBytes + ...
                (channel(channelNr).codePhase-1)*2), ...
                'bof');
        else
            fseek(fid, ...
                dataAdaptCoeff*(settings.skipNumberOfBytes + ...
                channel(channelNr).codePhase-1), ...
                'bof');
        end

        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1022) caCode(1023) caCode caCode(1) caCode(2)];

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;

        %C/No computation
        vsmCnt  = 0;CNo = 0;

        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods

            %To record MMT computation time
            % Added by Sergio Vicenzo - 5 Nov 2024
            tic 
            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 50) == 0)

                Ln=sprintf('\n');
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(settings.numberOfChannels),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt), ...
                    ' of ', int2str(codePeriods), ' msec',Ln...
                    'C/No: ',CNo,' (dB-Hz)'];

                try
                    waitbar(loopCnt/codePeriods, ...
                        hwb, ...
                        trackingStatus);
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            %% Read next block of data ------------------------------------------------
            % Record sample number (based on 8bit samples)
            % Edited by Sergio Vicenzo to also handle "int16" files
            % 12 Feb 2024
            if strcmp(settings.dataType,'int16') 
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff/2;
            else
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff;
            end

			% Find the size of a "block" or code period in whole samples
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;
            
            % Find the size of a "block" or code period in whole samples
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);

            % Read in the appropriate number of samples to process this
            % interation
            [rawSignal, samplesRead] = fread(fid, ...
                dataAdaptCoeff*blksize, settings.dataType);

            rawSignal = rawSignal';

            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;  % transpose vector
            end


            % If did not read in enough samples, then could be out of
            % data - better exit
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);
                return
            end

            %% Generate the carrier frequency to mix the signal to baseband -----------
            % Save remCarrPhase for current correlation
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;

            % Get the argument to sin/cos functions
            time    = (0:blksize) ./ settings.samplingFreq;
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));

            % Finally compute the signal to mix the collected data to
            % bandband
            carrsig = exp(1i .* trigarg(1:blksize));

            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            %% Compute multicorrelators -----------------------------------------------
            % Added by Sergio Vicenzo - Added on 5 Nov 2024
            chip_spacings = [-(flip(settings.chipspacing_dpe_precalc:...
                            settings.chipspacing_dpe_precalc:1)),0,...
                            settings.chipspacing_dpe_precalc:...
                            settings.chipspacing_dpe_precalc:1];
            count = 1;
            precalc_correlations    = ...
                zeros(length(chip_spacings),1);
            
            for spacing=(chip_spacings)
    
                % Define index into the code vector 
                delay_index = ...
                    remCodePhase-spacing : ...
                        codePhaseStep : ...
                        ((blksize-1)*codePhaseStep + ...
                       remCodePhase-spacing);

                caCode1 = caCode;
                tcodee = ceil(delay_index)+2;
                tcodee(tcodee==0)=tcodee(tcodee==0)+1;
    
                s = caCode1(tcodee);
    
                I = sum(s  .* iBasebandSignal);
                Q = sum(s  .* qBasebandSignal);
    
            % Store the correlations
                precalc_correlations(count,1) = ...
                    sqrt(I.^2 + Q.^2);
                precalc_correlations(count,2) = ...
                    I;
                precalc_correlations(count,3) = ...
                    Q;

                count=count+1;
            end % End for getting multicorr values
              
            % Get maximum corr value to initialize MMT
            [~,codedelay_index] = ...
                    max(precalc_correlations(:,1),[],'all','linear');

            LOS_codeDelay = (chip_spacings(codedelay_index));

            %% Multipath Mitigation Technology (MMT) ----------------------------------
            % Generate grid of LOS and reflected path delay
            codeDelay_NLOS = LOS_codeDelay-0.1:0.005:LOS_codeDelay+0.5;
            codeDelay_LOS  = LOS_codeDelay-0.1:0.005:LOS_codeDelay+0.1;
            MMT_cost_func3 = zeros(length(codeDelay_LOS),length(codeDelay_NLOS));

            NLOS_length=length(codeDelay_NLOS);
            LOS_length=length(codeDelay_LOS);
    
            % Use parallel computing to run MMT
            % Requires parallel computing toolbox
            parfor delay1_index=1:LOS_length
                delay1=codeDelay_LOS(delay1_index);

                for delay2_index=1:NLOS_length
                    delay2=codeDelay_NLOS(delay2_index);

                    if delay2 >= delay1 

                        delay_index_0 = ...
                        remCodePhase : ...
                            codePhaseStep : ...
                            ((blksize-1)*codePhaseStep + ...
                            remCodePhase);
    
               
                        delay_index_1 = ...
                        remCodePhase-delay1 : ...
                            codePhaseStep : ...
                            ((blksize-1)*codePhaseStep + ...
                            remCodePhase-delay1);
    
                        delay_index_2 = ...
                        remCodePhase-delay2 : ...
                            codePhaseStep : ...
                            ((blksize-1)*codePhaseStep + ...
                            remCodePhase-delay2);
    
                        delay_index_12 = ...
                        remCodePhase-(delay2-delay1)  : ...
                            codePhaseStep : ...
                            ((blksize-1)*codePhaseStep + ...
                            remCodePhase-(delay2-delay1) );

                        caCode1 = caCode;
                        tcodee_0 = ceil(delay_index_0)+2;
                        tcodee_0(tcodee_0==0)=tcodee_0(tcodee_0==0)+1;
                        tcodee_0(tcodee_0 > length(caCode1)) = length(caCode1);
                        
                        tcodee_1 = ceil(delay_index_1)+2;
                        tcodee_1(tcodee_1==0)=tcodee_1(tcodee_1==0)+1;
                        tcodee_1(tcodee_1 > length(caCode1)) = length(caCode1);
                

                        tcodee_2 = ceil(delay_index_2)+2;
                        tcodee_2(tcodee_2==0)=tcodee_2(tcodee_2==0)+1;
                        tcodee_2(tcodee_2 > length(caCode1)) = length(caCode1);
        
                        tcodee_12 = ceil(delay_index_12)+2;
                        tcodee_12(tcodee_12==0)=tcodee_12(tcodee_12==0)+1;
                        tcodee_12(tcodee_12 > length(caCode1)) = length(caCode1);
        
                        s_0 = caCode1(tcodee_0);
                        s_1 = caCode1(tcodee_1);
                        s_2 = caCode1(tcodee_2);
                        s_12 = caCode1(tcodee_12);
                        R_xm_1 = sum(iBasebandSignal.*s_1);
                        R_ym_1 = sum(qBasebandSignal.*s_1);

                        R_xm_2 = sum(iBasebandSignal.*s_2);
                        R_ym_2 = sum(qBasebandSignal.*s_2);
        
                        R_mm_12 = sum(s_0.*s_12);
                        R_mm_0  = sum(s_0.*s_0); 
        
                        R_1 = (sum(iBasebandSignal.*s_1) ...
                        + (sum(qBasebandSignal.*s_1)).*1i);
        
                        R_2 = (sum(iBasebandSignal.*s_2) ...
                        + (sum(qBasebandSignal.*s_2)).*1i);


                        % Apply amplitude constraint through Lagrange
                        % multiplier
                        X = (settings.MMT_const^4)*(norm(R_2,1)^2) ...
                            - (settings.MMT_const^2)*(norm(R_1,1)^2);
                            
                        Y = 2*(settings.MMT_const^2)*...
                         (R_mm_0*((norm(R_1,1)^2)+(norm(R_2,1)^2)) ...
                         - 2*R_mm_12*real(R_1*R_2));
        
                        Z = norm(R_mm_0*R_2 - R_mm_12*R_1,1)^2 - ...
                         (settings.MMT_const^2)*norm(R_mm_0*R_1 - ...
                         R_mm_12 - R_mm_12*R_2,1)^2;
    
             
                        gamma_1 = (-Y + sqrt(Y^2 - 4*X*Z))/(2*X);
        
                        a3 = (((R_mm_0 - gamma_1)* R_xm_1) - (R_xm_2*R_mm_12))/...
                            (R_mm_0^2 - gamma_1*R_mm_0*(1-settings.MMT_const^2)...
                            - (settings.MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                        b3 = (((R_mm_0 + ((settings.MMT_const^2) * gamma_1))...
                            * R_xm_2) - (R_xm_1*R_mm_12))/...
                            (R_mm_0^2 - gamma_1*R_mm_0*(1-settings.MMT_const^2)...
                            - (settings.MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                        c3 = (((R_mm_0 - gamma_1)* R_ym_1) - (R_ym_2*R_mm_12))/...
                            (R_mm_0^2 - gamma_1*R_mm_0*(1-settings.MMT_const^2)...
                            - (settings.MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                        d3 = (((R_mm_0 + ((settings.MMT_const^2) * gamma_1))...
                            * R_ym_2) - (R_ym_1*R_mm_12))/...
                            (R_mm_0^2 - gamma_1*R_mm_0*(1-settings.MMT_const^2)...
                            - (settings.MMT_const^2)*(gamma_1^2) - R_mm_12^2);


                        MMT_cost_func3(delay1_index,delay2_index) = ...
                            abs(2*real((a3-1i*c3)*...
                            (sum(iBasebandSignal.*s_1) + 1i*sum(qBasebandSignal.*s_1)))...
                            +2*real((b3-1i*d3)*...
                            (sum(iBasebandSignal.*s_2) + 1i*sum(qBasebandSignal.*s_2))) - ...
                            2*(a3*b3 + c3*d3)*R_mm_12 - (a3^2 + b3^2 + c3^2 +d3^2)*R_mm_0);

                    end 
   
                 end % End for delay2

            end % End for delay1

            % Maximise MMT cost function
            [MaxCorrValue,~] = max(MMT_cost_func3,[],'all','linear');
        
            [codeDelay_LOS_indx,codeDelay_NLOS_indx,~] = ...
                find(MMT_cost_func3==MaxCorrValue); 
                   
            % Obtain LOS code delay
            trackResults(channelNr).MMT_codeDelay_LOS(loopCnt)...
                 = mean(codeDelay_LOS(codeDelay_LOS_indx));
            LOS_delay = mean(codeDelay_LOS(codeDelay_LOS_indx));
    
            % Obtain reflected path code delay
            trackResults(channelNr).MMT_codeDelay_NLOS(loopCnt)...
                 = mean(codeDelay_NLOS(codeDelay_NLOS_indx));
    
            % Obtain relative delay of reflected path
            trackResults(channelNr).MMT_codeDelay_NLOS_relative(loopCnt)...
                =  abs(mean(codeDelay_NLOS(codeDelay_NLOS_indx)) - ...
                mean(codeDelay_LOS(codeDelay_LOS_indx)));
    
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
    
             % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 2;
            earlyCode   = caCode(tcode2);
    
            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 2;
            lateCode    = caCode(tcode2);
    
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 2;
            promptCode  = caCode(tcode2);
    
            remCodePhase = (tcode(blksize) + codePhaseStep) ...
                - settings.codeLength;

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);

            %% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);

            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            %% Find DLL error and update code NCO -------------------------------------
            % Disabled by Sergio Vicenzo - 5 Nov 2024
%             codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
%                 (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

            % DLL error now defined by LOS delay from MMT
            % Added by Sergio Vicenzo - 5 Nov 2024
            codeError = LOS_delay;

            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;

            %% Record various measures to show in postprocessing ----------------------
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            if (rem(loopCnt,settings.CNo.VSMinterval)==0)
                vsmCnt=vsmCnt+1;
                CNoValue=CNoVSM(trackResults(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    trackResults(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
                trackResults(channelNr).CNo.VSMValue(vsmCnt)=CNoValue;
                trackResults(channelNr).CNo.VSMIndex(vsmCnt)=loopCnt;
                CNo=int2str(CNoValue);
            end

            % Clear MMT variables after every loop
            % Added by Sergio Vicenzo - 5 Nov 2024
            clear corrs corrs2  index_forMaxValue index_LOS anIndex ...
                Early_corr_index  minVal LOS_delay codeDelay_NLOS_indx ...
                codeDelay_LOS_indx delay_indexxx x A b precalc_correlations...
                Late_corr_index 

            % Record MMT computation time for analysis
            % Added by Sergio Vicenzo - 5 Nov 2024
            trackResults(channelNr).MMT_time(loopCnt) = toc;

        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
