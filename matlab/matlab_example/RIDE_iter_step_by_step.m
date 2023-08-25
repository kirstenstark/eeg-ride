%function results = RIDE_iter(data,cfg)
%global temp555;


%% OUTCOMMENT IF RITER_iter- function is used
%        c = 62;
%        data =  squeeze(data(:,c,:));
%        cfg = cfg1;
%        j = 1;


%% Start of ride iter function

for sec = 1:1%default values and initial settings
    
        [d1,d2] = size(data);
        bd = cfg.bd;

        
        data0  = data;
        %the frequency band beyond 20Hz is not involved in the iterations
        %but will be preserved in the final iteration step
        
        
        data = filtering20(data,1,fix(10*20*cfg.re_samp*d1/1000)); % band-pass filters the data, but at 20x the frequency resolution as normal FFT
        % compare if this oversampling is really necessary or if we can use a standard filter function
        %
        for j = 1:cfg.comp_num
            max_latency(j) = max(cfg.comp.latency{j}(:));
            min_latency(j) = min(cfg.comp.latency{j}(:));
            length_c(j) = d1+max_latency(j)-min_latency(j);
            com_c(:,j) = zeros(d1,1);
            com_c1(:,:,j) = zeros(d1,d2);%stimulus-locked matrix
            amp_c(:,j) = zeros(d2,1);
        end
        
        % if isfield(cfg,'latency_a')
        %     ms = cfg.latency_a;
        %     max_latency_ms = max(ms(~isnan(ms)));
        %     min_latency_ms = min(ms(~isnan(ms)));
        %     length_ms = d1+max_latency_ms-min_latency_ms;
        %     com_ms = zeros(d1,1);
        %     com_ms1(:,:) = zeros(d1,d2);
        %     amp_ms = zeros(d2,1);
        %     d_ms = length(find(~isnan(ms)));
        %     ms_twd = cfg.ms_twd;
        % end

end

% Set stream_flow
trend_i = cfg.comp_num;stream_flow = 1:cfg.comp_num;
for c = 1:cfg.comp_num
   if cfg.comp.name{c}(1) == 'r' trend_i = c-1;end
end
if trend_i == cfg.comp_num stream_flow = [trend_i,1:trend_i-1];end % !!! might be needed if we have C component or for R, S
if trend_i == cfg.comp_num-1 stream_flow = [trend_i,1:trend_i-1,cfg.comp_num];end % !!! might be needed if we have C component or for R, S
if trend_i == 1 stream_flow = [2,1];end % !!! this is always the case for S, R

%% Second section: Code section (from here to end of function 
%for code = 1:1

    l1 = [];stop = 0;
    stop_c = zeros(cfg.comp_num,1);
            
 %%    
for iter = 1:cfg.inner_iter % !!! always 1:100 -> for loops ends around line 260
                fprintf('goes into loop')
                disp(iter);
                % %% ATTENTION: WHEN FOR-LOOP IS OUTCOMMENTED -> THIS SECTION NEEDS TO BE RUN TWICE, ONCE WITH ITER=1, ONCE WITH iter = 100
                % iter=1;
                % iter = 100; % only if loop is not run, first set to 1 to generate com_old, then set to 100

                % if iter == cfg.inner_iter stop = 1;end
                
                for track_conv = 1:1%track the convergence
                    if iter > 2 
                        for c = 1:cfg.comp_num
                            l1(iter-2,c) = sum(abs(com_c(:,c) - com_old(:,c)));
                            if iter >3 
                                if l1(iter-3,c)-l1(iter-2,c)<0.01*(l1(1,c)-l1(2,c)) stop_c(c) = 1;end%convergence is now defined as 0.01
                                if l1(iter-3,c)==l1(iter-2,c) stop_c(c) = 1;end
                            end
                        end
                        if cfg.comp_num == 1 stop_c(c) = 1;end
                    end
                    for c = 1:cfg.comp_num 
                        com_old(:,c) = com_c(:,c);
                         %decision of the termination of the interation of each component
                    end
                    
                end
                    
%                 if channel == 3 temp555(:,iter) = com_c(:,2);end
%                 temp555(iter) = com_c(500,2);
               
                % stream flow loop
              for c = stream_flow % - loop stops in line 178
              % c = 2; % only if loop is outcommented
                   if stop_c(c) == 0 % - loop stops in line 179
                        fprintf('goes into stream flow loop')
                        disp(iter)
                        temp = data;
                        for j = 1:cfg.comp_num
                            if j~=c temp = temp - com_c1(:,:,j);end
                            
                            % if isfield(cfg,'latency_a') temp = temp - com_ms1;end
                            
                        end

                            residue = temp;
                            temp = nan(length_c(c),d2);
                            for j = 1:d2
                                temp(1-cfg.comp.latency{c}(j)+max_latency(c):d1-cfg.comp.latency{c}(j)+max_latency(c),j) = residue(:,j);
                            end

                          %   if channel == 44 && iter==inner_iter&&c==2 figure;subplot(1,2,1);plot(temp);hold on;plot(mean_nan(temp,2),'linewidth',6);end

                           temp = RIDE_detrend(temp,[1 1+fix((cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1))*bd) ... % !!! onset and offset of window edges
                           fix((cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1))*(1-bd)) cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1)]+max_latency(c)+cfg.comp.twd{c}(1));
                           
%                                         % %% RIDE detrend step_by_step
%                                         % 
%                                         % %%%function f = RIDE_detrend(data,twd)
%                                         % 
%                                                                    data = temp
%                                                                    twd = [1 1+fix((cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1))*bd) ... % !!! onset and offset of window edges
%                                                                     fix((cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1))*(1-bd)) cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1)]+max_latency(c)+cfg.comp.twd{c}(1)
%                                         
%                                         % index = isnan(data);
%                                         % data(index) = 0;
%                                         % [d1,d2,d3] = size(data);
%                                         % for j = 1:d2
%                                         %     for k = 1:d3
%                                         %         a = mean(data(twd(1):twd(2),j,k));
%                                         %         b = mean(data(twd(3):twd(4),j,k));
%                                         %         temp = ([1:d1])*(b-a)/((twd(4)+twd(3))/2-(twd(2)+twd(1))/2);
%                                         %         temp = temp - mean(
%                                         %         data(:,j,k) = data(:,j,k) - temp';
%                                         %     end
%                                         % end
%                                         % f = data;
%                                         % f(index) = nan;
%                                         
%                                         
%                                         
%                                         
%                                          
%                                         index = isnan(data);
%                                         data(index) = 0;
%                                         [d1,d2,d3] = size(data);
%                                         d = ((twd(4)+twd(3))/2-(twd(2)+twd(1))/2);
%                                         
%                                         temp0 = repmat([1:d1]',[1,d2,d3]);
%                                         
%                                                 a = mean(data(twd(1):twd(2),:,:),1);
%                                                 b = mean(data(twd(3):twd(4),:,:),1);
%                                                 temp = (b-a)/d;
%                                                 data = data - temp0.*temp(ones(1,d1),:,:);
%                                         
%                                                 temp = mean(data(twd(1):twd(2),:,:),1);
%                                                 data = data - temp(ones(1,d1),:,:);
%                                                 
%                                         f = data;
%                                         f(index) = nan;
%                                         %

                        temp0 = median_2d(temp');
                        temp0([1:cfg.comp.twd{c}(1)+max_latency(c),cfg.comp.twd{c}(2)+max_latency(c):length_c(c)])=0;
                        
                        temp0((cfg.comp.twd{c}(1)+max_latency(c)):(cfg.comp.twd{c}(2)+max_latency(c)))=...
                            temp0((cfg.comp.twd{c}(1)+max_latency(c)):(cfg.comp.twd{c}(2)+max_latency(c))).*RIDE_tukey(cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1)+1,bd*2);
                        
                        temp1 = temp0(:,ones(1,d2));
                        temp1(isnan(temp)) = nan;
                        
                        %temp1(:, 45) = nan;
                        
                        temp1 = (reshape(temp1(~isnan(temp1)),d1,d2)); % We're trying to keep any NaNs in the data
                        com_c(:,c) = temp0(max_latency(c)+1:max_latency(c)+d1);
                        com_c1(:,:,c) = temp1; % Errors because for the r component, the 1st size of `temp1` (determined by `length_c`) is not 650
                    end % - end of loop (line 98)
             end % - end of stream flow loop (line 96)


             if mean(stop_c) == 1 stop = 1;end
             if stop == 1
                    for c = stream_flow
                        
                        temp = data;
                        for j = 1:cfg.comp_num
                            if j~=c temp = temp - com_c1(:,:,j);end
                            
                            % if isfield(cfg,'latency_a') temp = temp - com_ms1;end
                            
                        end
                        
                        residue = temp;
                        temp = nan(length_c(c),d2);
                        for j = 1:d2
                            temp(1-cfg.comp.latency{c}(j)+max_latency(c):d1-cfg.comp.latency{c}(j)+max_latency(c),j) = residue(:,j);
                        end
                        
                        if cfg.final == 1
                            temp(isnan(temp))=0;

                                amp_c(:,c) = mean(com_c(cfg.comp.twd{c}(1):cfg.comp.twd{c}(2),c*ones(1,d2)).*temp(max_latency(c)+[cfg.comp.twd{c}(1):cfg.comp.twd{c}(2)],:));

                        end
                        
                        
                    end
                    
                    % for e = 1:fix(isfield(cfg,'latency_a'))
                    %     temp = data;
                    %     for j = 1:cfg.comp_num
                    %         temp = temp - com_c1(:,:,j);
                    %     end
                        
                    %     temp = temp - com_ms1;
                        
                        
                        
                    %     residue = temp;
                    %     temp = nan(length_ms,d_ms);l=0;
                    %     for j = 1:d2
                    %         for jj = 1:length(find(~isnan(ms(:,j))))
                    %             l = l+1;
                    %             temp(1-ms(jj,j)+max_latency_ms:d1-ms(jj,j)+max_latency_ms,l) = residue(:,j);
                    %         end
                    %     end
                        
                        
                    %     temp = RIDE_detrend(temp,[1 fix((ms_twd(2)-ms_twd(1))*bd) ...
                    %         fix((ms_twd(2)-ms_twd(1))*(1-bd)) ms_twd(2)-ms_twd(1)]+max_latency_ms+ms_twd(1));



                        
                    %     temp0 = median_2d(temp');
                    %     temp0(max_latency_ms+1:max_latency_ms+d1) = temp0(max_latency_ms+1:max_latency_ms+d1) + com_ms;
                    %     temp0([1:ms_twd(1)+max_latency_ms,ms_twd(2)+max_latency_ms:length_ms])=0;
                    %     temp0((ms_twd(1)+max_latency_ms):(ms_twd(2)+max_latency_ms))=...
                    %         temp0((ms_twd(1)+max_latency_ms):(ms_twd(2)+max_latency_ms)).*RIDE_tukey(ms_twd(2)-ms_twd(1)+1,bd*2);
                        
                    %     temp1 = zeros(d1,d2);
                    %     for j = 1:d2
                    %         for jj = 1:length(find(~isnan(ms(:,j))))
                    %             temp1(:,j) = temp1(:,j)+temp0(1-ms(jj,j)+max_latency_ms:d1-ms(jj,j)+max_latency_ms);
                    %         end
                    %     end
                    %     com_ms = temp0(max_latency_ms+1:max_latency_ms+d1);
                    %     com_ms1(:,:) = temp1;
                    % end
%                    if cfg.channel == 32  figure;imagesc(com_ms1);figure;plot(mean(com_ms1,2));hold on;plot(com_ms,'r');end
                    
                    break;
                end 
                

    end % END of ITER LOOP (Line 67)
    
  %% 
            
            
            for last_iter = 1:1
                
                if cfg.final == 1
                    for ongoing = 1:1 %release time window function and detrending
                        %allocate trend to the last C component
                         for c = stream_flow
                                temp = data0;
                                for j = 1:cfg.comp_num
                                    if j~=c temp = temp - com_c1(:,:,j);end
                                    if exist('ms','var') temp = temp - com_ms1;end % !!! probably only for microsaccades
                                end


                                    residue = temp;
                                    temp = nan(length_c(c),d2);
                                    for j = 1:d2
                                        temp(1-cfg.comp.latency{c}(j)+max_latency(c):d1-cfg.comp.latency{c}(j)+max_latency(c),j) = residue(:,j);
                                    end
                                    
                                    

                                temp0 = mean_nan(temp,2);
                                
                                if c==stream_flow(1)
                                    temp0([1:cfg.comp.twd{c}(1)+max_latency(c)])=0;
                                    tem = RIDE_tukey(cfg.comp.twd{c}(2)-cfg.comp.twd{c}(1)+1,bd*2);tem = tem(1:fix(end/2));
                                    temp0((cfg.comp.twd{c}(1)+max_latency(c)):(cfg.comp.twd{c}(1)+max_latency(c))+length(tem)-1)=...
                                        temp0((cfg.comp.twd{c}(1)+max_latency(c)):(cfg.comp.twd{c}(1)+max_latency(c))+length(tem)-1).*tem;
                                end
                                
                                temp1 = temp0(:,ones(1,d2));
                                temp1(isnan(temp)) = nan;
                                temp1 = (reshape(temp1(~isnan(temp1)),d1,d2));
                                com_c(:,c) = temp0(max_latency(c)+1:max_latency(c)+d1);
                                com_c1(:,:,c) = temp1;

                         end
                         
                        %  for e = 1:fix(isfield(cfg,'latency_a'))
                        %      temp = data0;
                        %      for j = 1:cfg.comp_num
                        %          temp = temp - com_c1(:,:,j);
                        %      end
                             
                        %      temp = temp - com_ms1;
                             
                             
                             
                        %      residue = temp;
                        %      temp = nan(length_ms,d_ms);l=0;
                        %      for j = 1:d2
                        %          for jj = 1:length(find(~isnan(ms(:,j))))
                        %              l = l+1;
                        %              temp(1-ms(jj,j)+max_latency_ms:d1-ms(jj,j)+max_latency_ms,l) = residue(:,j);
                        %          end
                        %      end
                             
                             
                             
                             
                        %      temp0 = mean_nan(temp,2);
                        %      temp0(max_latency_ms+1:max_latency_ms+d1) = temp0(max_latency_ms+1:max_latency_ms+d1) + com_ms;
                             
                        %      temp1 = zeros(d1,d2);
                        %      for j = 1:d2
                        %          for jj = 1:length(find(~isnan(ms(:,j))))
                        %              temp1(:,j) = temp1(:,j)+temp0(1-ms(jj,j)+max_latency_ms:d1-ms(jj,j)+max_latency_ms);
                        %          end
                        %      end
                        %      com_ms = temp0(max_latency_ms+1:max_latency_ms+d1);
                        %      com_ms1(:,:) = temp1;
                        %  end
%                          if cfg.channel == 32  figure;imagesc(com_ms1);figure;plot(mean(com_ms1,2));hold on;plot(com_ms,'r');end
                         for c = stream_flow(1)
                                temp = data0;
                                for j = 1:cfg.comp_num
                                    if j~=c temp = temp - com_c1(:,:,j);end
                                    %  if exist('ms','var') temp = temp - com_ms1;end
                                end


                                    residue = temp;
                                    temp = nan(length_c(c),d2);
                                    for j = 1:d2
                                        temp(1-cfg.comp.latency{c}(j)+max_latency(c):d1-cfg.comp.latency{c}(j)+max_latency(c),j) = residue(:,j);
                                    end
                                    
                                    % !!! Here we don't apply the tukey window apparently

                                temp0 = mean_nan(temp,2);                                
                                temp1 = temp0(:,ones(1,d2));
                                temp1(isnan(temp)) = nan;
                                temp1 = (reshape(temp1(~isnan(temp1)),d1,d2));
                                com_c(:,c) = temp0(max_latency(c)+1:max_latency(c)+d1);
                                com_c1(:,:,c) = temp1;

                         end
                         
                         

                    end
                end
                

            end
            
            
  
            

            % if exist('ms','var') results.comp_ms = com_ms;results.comp_ms1 = mean(com_ms1,2);end
            results.amp = amp_c;
            results.comp = com_c;
            results.comp1 = permute(mean(com_c1,2),[1,3,2]);
            results.iter = iter;
            results.l1 = l1;
            results.trend_c = stream_flow(1);
            

%end


