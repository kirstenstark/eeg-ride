function results = RIDE_call(data,cfg)

% for copyright = 1:1
% % Code Author: Guang Ouyang, HKBU, 2010,
% % The RIDE method developers: Guang Ouyang, Werner Sommer, Changsong Zhou (alphabetical order)
% % Copyright (C) <2010>  Guang Ouyang, Werner Sommer, Changsong Zhou
% % 
% %     This program is free software: you can redistribute it and/or modify
% %     it under the terms of the GNU General Public License as published by
% %     the Free Software Foundation, either version 3 of the License, or
% %     (at your option) any later version.
% % 
% %     This program is distributed in the hope that it will be useful,
% %     but WITHOUT ANY WARRANTY; without even the implied warranty of
% %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% %     GNU General Public License for more details.
% % 
% %     You should have received a copy of the GNU General Public License
% %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% end


cfg0 = cfg; %save the initial configurations
if length(cfg.comp.name)==1 
    results.erp = mean(data,3);return;%if only separate one component, it is ERP
end

% %showing the progress
% if cfg.prg == 1
% figure;set(gcf,'menubar','none');axis([0,100,0.5,1.5]);axis off;temp = get(gcf,'position');set(gcf,'position',[temp(1),temp(2),temp(3),150]);pause(0.001);
% end


for section = 1:1  %preparation, down_sampling-----------------------------------------------------------------
    [d1,d2,d3] = size(data);
    epoch_length = d1;
    erp = mean(data,3);results.erp = erp;%original ERP
    results.latency0 = cfg.comp.latency;%save the original latency information, i.e., RT

    % %------------down sample data---------------------
    % rs = (cfg.re_samp/cfg.samp_interval);
    % data = data(round(linspace(1,d1,fix(d1/rs))),:,:);
    
    % %----------only if using specified template to measure C------------
    % if isfield(cfg,'template')
    %     if strcmp(cfg.template.method,'g_mean')
    %         load(cfg.temp_path);
    %         template = template(round(linspace(1,d1,fix(d1/rs))),:,:);
    %     end
    % end
    % %-------------------------------------------------------------

    [d1,d2,d3] = size(data);%new size after down samping
    
    
    % %----------only for microsaccade-------------------------------
    % if isfield(cfg,'latency_a') 
    %     cfg.latency_a = round(cfg.latency_a/cfg.re_samp);
    %     cfg.latency_a = round(cfg.latency_a-median(cfg.latency_a(~isnan(cfg.latency_a))));
    %     cfg.ms_twd = fix((cfg.ms_twd-cfg.epoch_twd(1))/cfg.re_samp)+[1,-1];
    % end
    % %------------------------------------------------------------------
    
    for j = 1:cfg.comp_num %unified the unit of time window and latency information
        if cfg.comp.latency{j} == 0 cfg.comp.latency{j} = zeros(d3,1);end
        if ~ischar(cfg.comp.latency{j})
            if strcmp(cfg.comp.name{j},'r')
                cfg.comp.twd{j} = cfg.comp.twd{j} + median(results.latency0{j});
                cfg.comp.twd{j}(cfg.comp.twd{j}<cfg.rwd) = cfg.rwd;%left boundary of twd for R not less then 200 ms % !!! warn the user about this
                cfg.comp.twd{j}(cfg.comp.twd{j}>cfg.epoch_twd(2)) = cfg.epoch_twd(2); % !!! warn the user about this
            end
            cfg.comp.latency{j} = ((cfg.comp.latency{j})/cfg.re_samp); % !!! cfg.re_samp is the sample interval
            cfg.comp.latency{j} = round(cfg.comp.latency{j} - median(cfg.comp.latency{j}));
        end 
            cfg.comp.twd{j} = fix((cfg.comp.twd{j}-cfg.epoch_twd(1))/cfg.re_samp)+[1,-1];%convert component time window to sampling unit % !!! divide by sample interval (e.g. 2 ms)
            % if cfg.comp.twd{j}(2)<cfg.comp.twd{j}(1) cfg.comp.twd{j} = [cfg.comp.twd{j}(2)-1 cfg.comp.twd{j}(2)]';end % !!! a bit weird, rather error if twd(2) < twd(1)
    end
    
    % %--------specify the searching duration-----------------------
    % if isfield(cfg,'dur')
    %     for j = 1:length(cfg.comp.name)
    %         cfg.dur{j} = fix(cfg.dur{j}/cfg.re_samp);
    %     end
    % else
    %     for j = 1:length(cfg.comp.name)
    %        cfg.dur{j} = round(((cfg.comp.twd{j}(2)-cfg.comp.twd{j}(1)))/2);
    %     end
    % end
    % %----------------------------------------------------------------

end
    
        

    
% for section = 1:1%initial estimation of the latency of C---------------------------------------
%     for initial_c = 1:1 
%             n_of_c = 0;c_i = 0;
%             for j = 1:cfg.comp_num
%                 if ischar(cfg.comp.latency{j})
%                         if cfg.prg == 1 disp(['woody_for_',cfg.comp.name{j}]);end
%                         n_of_c = n_of_c + 1;c_i(n_of_c) = j;
%                         temp = 1:d2;temp1 = cfg.comp.twd{j};
%                             if isfield(cfg,'template')
%                                 if strcmp(cfg.template.method,'g_mean')
%                                     cfg.temp = template(temp1(1):temp1(2),temp);
%                                 end
%                                 if isfield(cfg.template,'chan')
%                                     temp = cfg.template.chan;
%                                     if isfield(cfg.template,'hann_amp')
%                                         cfg.template.hann_amp = cfg.template.hann_amp(cfg.template.chan);
%                                     end
%                                 end
%                             end
%                             %-------using Woody's method by default
%                         cfg.comp.latency{j} = woody(data(temp1(1):temp1(2),temp,:),cfg,cfg.dur{j});%
                        
%                 end
%             end
%     end
% end



for section = 1:1 %RIDE iteration
    stop=0;
    
    % for j = 1:n_of_c% track latency evolution of C components
    %     latency_i{j} = cfg.comp.latency{c_i(j)};latency_i{j} = latency_i{j}(:);
    %     l_change(:,j) = ones(d3,1);%track for evolution of the latency in order to terminate the iteration
    %     c_change(:,j) = ones(d3,1);%track for evolution of the correlation in order to terminate the iteration
    % end

    % if cfg.prg == 1 disp('RIDE decomposition: ');end
    outer_iter = 4;if n_of_c == 0 outer_iter = 1;end%outer iteration is empericaly limited to 4, but if no c component, no need to do outer iteration
    % for iter = 1:outer_iter
        
        % %report the progress
        %  if n_of_c ~= 0 
        %      prog = fix(100*(1-sum(mean(l_change,2))/d3*(10-iter)/10));
        %      if cfg.prg == 1 
        %          barh(prog);axis([0,100,0.5,1.5]);axis off;text(10,1.5,strcat('iteration','--',num2str(prog),'%done'));pause(0.001);
        %          fprintf(strcat('iteration',num2str(iter),'--',num2str(prog),'%%done\n'));
        %      end
        %      if (fix(100*(1-sum(mean(l_change,2))/d3*(10-iter)/10)))>=99
        %          stop=1;%stop the iteration when more 99% of the single trial latency do not change
        %      end
        %  end
        % if cfg.prg == 1 fprintf('iteration step for each channel:\n');end
        % for j = 1:n_of_c        latency_i{j}(:,iter) = cfg.comp.latency{c_i(j)};         end
         if iter == outer_iter stop=1;end

         
         for sec_RIDE_inner_iter = 1:1%RIDE_RIDE_RIDE_RIDE_RIDE_RIDE_RIDE_RIDE
             cfg1 = cfg;
             cfg1.final = stop;
             cfg1.inner_iter = 100;%this is to safeguard the extreme
             %iterations (usually inner iteration stops before 20)
             
             c_l = zeros(d1,cfg.comp_num,d2);%c_l is the latency synced RIDE component
             c_sl = c_l;%c_l is the stimulus synced RIDE component
             for c = 1:d2 % loops over channels
                 
                 rst = RIDE_iter(squeeze(data(:,c,:)),cfg1);%RIDE decomposition
                 
                 c_l(:,:,c) = rst.comp;c_sl(:,:,c) = rst.comp1;
                 
                 if stop==1  amp(:,c,:) = permute(rst.amp,[1,3,2]);end
                 
                %  if isfield(cfg,'latency_a')
                %      comp_ms(:,c) = rst.comp_ms;comp_ms1(:,c) = rst.comp_ms1;
                %  end%only for microsaccade
                 
                 
                  if cfg.prg == 1 fprintf(strcat(num2str(rst.iter),'|'));end
             end
             if cfg.prg == 1 fprintf('\n');end
             comp = permute(c_l,[1,3,2]);comp1 = permute(c_sl,[1,3,2]);
         end
         %RIDE_RIDE_RIDE_RIDE_RIDE_RIDE_RIDE_RIDE
         
         
         if cfg.prg == 1 fprintf('\n');end

        

         for section1 = 1:1%%%%%%%%%%%latency estimation of C


                        for cc = 1:n_of_c
                            clear tem;
                            no_p{cc} = ones(d3,1);%check there are peaks in x-corr curves or not
                            %low pass filter the template only for
                            %calculation of cross-correlation curve
                            temp1 = filtering20(comp(:,:,c_i(cc)),1,round(10*cfg.high_cutoff*d1*cfg.re_samp/1000));
                            for l = 1:d3 %trial_index
                                if l_change(l,cc)==1
                                    temp2 = data(:,:,l);
                                    for j = 1:cfg.comp_num
                                        if j~=c_i(cc)
                                        temp2 = temp2 - move2(comp(:,:,j),cfg.comp.latency{j}(l),'1d');
                                        end
                                    end
                                    if isfield(cfg,'latency_a')%remove the ms component
                                        for a = 1:length(find(~isnan(cfg.latency_a(:,l))))
                                                temp2 = temp2 - move2(comp_ms,cfg.latency_a(a,l),'1d'); 
                                        end
                                    end
                                    
                                    %low pass filter the data only for calculation of cross-correlation curve
                                    temp2 = filtering20(temp2,1,round(10*cfg.high_cutoff*d1*cfg.re_samp/1000));

                                    %remove the linear trend to safeguard
                                    %that the latency estimation is not affected by drifting
                                    temp2 = RIDE_detrend(temp2,[cfg.comp.twd{c_i(cc)}(1),cfg.comp.twd{c_i(cc)}(1)+fix((cfg.comp.twd{c_i(cc)}(2)-cfg.comp.twd{c_i(cc)}(1))*cfg.bd),...
                                        cfg.comp.twd{c_i(cc)}(1) + fix((cfg.comp.twd{c_i(cc)}(2)-cfg.comp.twd{c_i(cc)}(1))*(1-cfg.bd)), cfg.comp.twd{c_i(cc)}(2)]);

                                    
                                    tem0 = 1:d2; 
                                    if isfield(cfg,'template') %only for specified tempalte matching
                                        if isfield(cfg.template,'chan')
                                            tem0 = cfg.template.chan;
                                        end
                                    end
                                    for c = 1:length(tem0)%cross correlation
                                        tem(:,c) = xcov(temp2(:,tem0(c)),temp1(:,tem0(c)),fix(size(temp2,1)/2),cfg.xc);
                                    end
%                                     tem0 = mean(tem(fix(size(tem,1)/2-cfg.dur{c_i(cc)})+1:fix(size(tem,1)/2+cfg.dur{c_i(cc)}),:));

                                        temp11{cc}(:,l) = mean(tem,2);

                                    %low pass filter the cross correlation
                                    %curve
                                    temp11{cc}(:,l) = filtering10(temp11{cc}(:,l),1,round(10*cfg.high_cutoff*size(temp2,1)*cfg.re_samp/1000));
                                end
                            end
                            %detrend the cross-correlation curve to make
                            %sure the latency will not be found on the boundaries
                            temp11{cc} = RIDE_detrend(temp11{cc}(:,:),[1,fix(size(temp11{cc},1)*cfg.bd),fix(size(temp11{cc},1)*(1-cfg.bd)),size(temp11{cc},1)]);
                            temp = [];
                            temp_11 = temp11{cc}(fix(size(temp11{cc},1)/2-cfg.dur{c_i(cc)})+1:fix(size(temp11{cc},1)/2+cfg.dur{c_i(cc)}),:);
                            
                            
                            if strcmpi(cfg.latency_search,'most_prob') %search the nearest peak from the most probable estimation
                                for j = 1:d3 temp(j) = nearest_latency(temp_11(:,j),cfg.comp.latency{c_i(cc)}(j)+fix(size(temp_11,1)/2));end
                            end
                            if strcmpi(cfg.latency_search,'all')%search the largest peak
                                for j = 1:d3 temp(j) = find_peak(temp_11(:,j));end
                            end

                            
                            % randomly assign the latency of the trials without x-corr peak
                            for j = 1:d3 
                                if no_peak(temp_11(:,j))==0 no_p{cc}(j) = 0;end
                            end
                                  
                            temp(no_p{cc}==0) = round(randn(length(find(no_p{cc}==0)),1)*std(temp(no_p{cc}==1))) + fix(size(temp_11,1)/2);
                            
                            temp(temp<1)=1;temp(temp>size(temp_11,1)) = size(temp_11,1);%make sure not exceed boundary
                           for j = 1:d3 corr_i{cc}(j,iter) = temp_11(temp(j),j);end %track the correlation values
                           

                           temp = round(temp-median(temp)); %covert the C latencies to relative values by subtracting the median


                           %tack the latency evolution and correlation values (if the evolution returns, then stop updating)   
                           if iter>1
                                for j = 1:d3
                                    if (temp(j)-latency_i{cc}(j,iter))*(latency_i{cc}(j,iter)-latency_i{cc}(j,iter-1))<=0&&l_change(j,cc)==1 l_change(j,cc) = 0;end
                                    c_change(j,cc)=1;
                                    for jj = 1:iter-1 
                                        if corr_i{cc}(j,jj)>=corr_i{cc}(j,iter)    c_change(j,cc)=0;        end
                                    end
                                end
                             end

                            index = (l_change(:,cc)==1&c_change(:,cc)==1);
                            cfg.comp.latency{c_i(cc)}(index) = temp(index);
                            cfg.comp.latency{c_i(cc)} = round(cfg.comp.latency{c_i(cc)}-median(cfg.comp.latency{c_i(cc)}));

                        end
                       

         end
            
            if stop==1;break;end
            
            
%             results.temp_c(:,iter) = comp(:,9,1);
    end
    if cfg.prg == 1 fprintf('100%%done\n');close(gcf);end
end




for section = 1:1%%%%%final data%%%%

    %if data has been down sampled, apply interpolation to restore the
    %original resolution
    
    %and re-baselining
    results.erp_new = 0;
    results.residue = erp;
    % if isfield(cfg,'latency_a') %only for microsaccades
    %     results.ms = baseline(interp2d(comp_ms,round(linspace(1,epoch_length,d1)),1:epoch_length,'spline'));
    %     results.ms_sl = baseline(interp2d(comp_ms1,round(linspace(1,epoch_length,d1)),1:epoch_length,'spline'));
    % end
    bl_wd = fix(-cfg.epoch_twd(1)/cfg.samp_interval)+1:fix(-cfg.epoch_twd(1)/cfg.samp_interval+cfg.bl/cfg.samp_interval);%baseline time window
    for j = 1:cfg.comp_num
        component(:,:,j) = interp2d(comp(:,:,j),round(linspace(1,epoch_length,d1)),1:epoch_length,'spline'); % !!! only relevant if data has been downsampled
        component(:,:,j) = baseline(component(:,:,j),bl_wd);
        component1(:,:,j) = interp2d(comp1(:,:,j),round(linspace(1,epoch_length,d1)),1:epoch_length,'spline');
        component1(:,:,j) = baseline(component1(:,:,j),bl_wd);
        results.residue = results.residue - component1(:,:,j);
        eval(['results.',cfg.comp.name{j},' = component(:,:,j);']);
        eval(['results.',cfg.comp.name{j},'_sl = component1(:,:,j);']);
        eval(['results.latency_',cfg.comp.name{j},' = cfg.comp.latency{j}*cfg.re_samp;']);
        eval(['results.amp_',cfg.comp.name{j},'=amp(:,:,j);']);
    end
    
    % if isfield(cfg,'latency_a') %only for microsaccades
    %     results.residue = results.residue - results.ms_sl;
    % end
    
    eval(['results.',cfg.comp.name{1},' = baseline(results.',...
        cfg.comp.name{1},' ,bl_wd) + repmat(mean(erp(bl_wd,:)),[epoch_length,1]);']); % !!! Un-does baseline correction by adding basline window from ERPs
    eval(['results.',cfg.comp.name{1},'_sl = baseline(results.',...
        cfg.comp.name{1},'_sl,bl_wd) + repmat(mean(erp(bl_wd,:)),[epoch_length,1]);']);
% %     

    eval(['results.',cfg.comp.name{rst.trend_c},' = baseline(results.',...
        cfg.comp.name{rst.trend_c},' + results.residue,bl_wd);']);
    eval(['results.',cfg.comp.name{rst.trend_c},'_sl = baseline(results.',...
        cfg.comp.name{rst.trend_c},'_sl + results.residue,bl_wd);']);
    
    if cfg.comp_num == 1 
        bl_wd = 1:-fix(cfg.epoch_twd(1)/cfg.samp_interval);
        eval(['results.',cfg.comp.name{1},' = baseline(results.',...
        cfg.comp.name{1},' + results.residue,bl_wd) + repmat(mean(erp(bl_wd,:)),[epoch_length,1]);']);
        eval(['results.',cfg.comp.name{1},'_sl = baseline(results.',...
            cfg.comp.name{1},'_sl + results.residue,bl_wd) + repmat(mean(erp(bl_wd,:)),[epoch_length,1]);']);
    end
    for j = 1:cfg.comp_num eval(['results.erp_new = results.erp_new + results.',cfg.comp.name{j},';']);end
    % if n_of_c~=0 results.latency_i = latency_i;results.no_p = no_p;end
    results.cfg = cfg0;
%     results.cfg1 = cfg;
%     results.l1 = l1;
    % if exist('corr_i','var') results.corr_i = corr_i;end



    
    

end


