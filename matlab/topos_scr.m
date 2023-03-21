function topos_scr(data,tpt,comp,chanlocs,varargin)



d1 = length(data);
d2 = length(tpt);


tpt1 = fix(linspace(1,size(data{1},1),d2));
for j = 1:d1
    for k = 1:d2
        subplot(d1,d2,k+(j-1)*d2);
        
                
        topoplot(data{j}(tpt1(k),:),chanlocs,varargin{:});
        
        
        axis on;set(gca,'DataAspectRatioMode','auto');
                set(gca,'PlotBoxAspectRatioMode','auto');
        axis([-0.7,0.7,-0.7,0.7]);axis off;axis tight;
        p = get(gca,'position');
        set(gca,'position',[p(1),p(2),0.06,0.19]);
                
                
        if j==1 text(-0.2,0.7,[num2str(fix(tpt(k))),'ms']);end
        if k==1 text(-1,0,comp{j});end
    end
end
