function [] = PlotRe(Results,t,xg,StructureName,mg,plots)

mvt = size(Results.Yta,2);
T = max(t);
switch plots
    case 'plot'
        %% Comparing accelerations
        figure(1)
        LB=min(min(Results.Ya/9.807))-0.2*max(max(Results.Ya/9.807));
        UB=max(max(Results.Ya/9.807))+0.2*max(max(Results.Ya/9.807));
        for i=1:mvt
            subplot(mvt,1,i),
            if i>1
                plot(t,Results.Ya(:,mvt-i+1)/9.807,'LineStyle','-','Color','b')
                hold on
                plot(t,Results.Yta(:,mvt-i+1)/9.807,'LineStyle','--','Color','r')
            else
                plot(t,Results.Yta(:,mvt-i+1)/9.807,'LineStyle','--','Color','r')
            end
            if i == 1
                title('Floors Acceleration')
            end
            if i==mvt
                xlabel('time [sec]','FontSize',8)
                legend('Without-TMD','With-TMD')
            end
            ylabel('acc[g]','FontSize',8)
            axis([0 T LB UB])
            grid on
        end
        %% Comparing displacements
        figure(2)
        LB=min(min(Results.Yd*100))-0.2*max(max(Results.Yd*100));
        UB=max(max(Results.Yd*100))+0.2*max(max(Results.Yd*100));
        for i=1:mvt
            subplot(mvt,1,i),
            if i>1
                plot(t,Results.Yd(:,mvt-i+1)*100,'LineStyle','-','Color','b')
                hold on
                plot(t,Results.Ytd(:,mvt-i+1)*100,'LineStyle','--','Color','r')
            else
                plot(t,Results.Ytd(:,mvt-i+1)*100,'LineStyle','--','Color','r')
            end
            if i == 1
                title('Floors Displacements')
            end
            if i==mvt
                xlabel('time [sec]','FontSize',8)
                legend('Without-TMD','With-TMD')
            end
            ylabel('dis[cm]','FontSize',8)
            axis([0 T LB UB])
            grid on
        end
        %% Plot Displacement and Acceleration on structure
        figure('Position',[500   100   560   800])
        SH = 3.3*(0:mvt)*mg;
        UB=max([1.2*max(max(abs(Results.Yd*100))) 1.2*max(max(abs(Results.Ytd*100)))]);
        LB = -UB;

        UBa=max([1.2*max(max(abs(Results.Ya./9.807))) 1.2*max(max(abs(Results.Yta./9.807)))]);
        LBa=-UBa;

        % ------------ Plot Displacement --------------------------
        subplot(1,2,1)

        plot([0 max(Results.Yd)]*100,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        plot([0 min(Results.Yd)]*100,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        plot([0 max(Results.Ytd)]*100,SH,'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        plot([0 min(Results.Ytd)]*100,SH,'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        axis([LB UB 0 max(SH)])
        yticks(SH)
        yticklabels({0:mvt-1,'TMD'})
        xlabel('Dis _{[Cm]}')
        ylabel('Storey')
        title('Displacement')
        grid on
        legend('Without-TMD','','With-TMD','')
        % ------------ Plot Acceleration --------------------------
        subplot(1,2,2)
        plot([0 max(Results.Ya)]./9.807,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        plot([0 min(Results.Ya)]./9.807,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');

        plot([0 max(Results.Yta)]./9.807,SH,'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        plot([0 min(Results.Yta)]./9.807,SH,'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        % legend('Without-TMD','','With-TMD','')
        axis([LBa UBa 0 max(SH)])
        yticks(SH)
        yticklabels({0:mvt-1,'TMD'})
        grid on
        xlabel('Acc _{[g]}')
        ylabel('Storey')
        title('Acceleration')

        %% Plort Drift and Force
        figure('Position',[500   100   560   800])
        SH = 3.3*(0:mvt)*mg;
        UB=max([1.2*max(max(abs(Results.Ydrift*100))) 1.2*max(max(abs(Results.Ydrift*100)))]);
        LB = -UB;

        UBa=max([1.2*max(max(abs(Results.Shear))) 1.2*max(max(abs(Results.Shear)))]);
        LBa=-UBa;

        % ------------ Plot Drift --------------------------
        subplot(1,2,1)

        plot([0 max(Results.Ydrift)]*100,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        plot([0 min(Results.Ydrift)]*100,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        plot([0 max(Results.Ytdrift)]*100,SH(1:mvt),'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        plot([0 min(Results.Ytdrift)]*100,SH(1:mvt),'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        axis([LB UB 0 max(SH(1:mvt))])
        yticks(SH(1:mvt))
        yticklabels({0:mvt-1})
        xlabel('Drift _{%}')
        ylabel('Storey')
        title('Drift')
        grid on
        % legend('Without-TMD','','With-TMD','')
        % ------------ Plot Shear --------------------------
        subplot(1,2,2)

        plot(max(Results.Shear),SH(2:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        hold on
        plot(min(Results.Shear),SH(2:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
        plot(max(Results.Sheart(:,1:mvt-1)),SH(2:mvt),'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
        plot(min(Results.Sheart(:,1:mvt-1)),SH(2:mvt),'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');

        yticks(SH(2:mvt))
        yticklabels({1:mvt-1})
        xlabel('Force _{[Kn]}')
        ylabel('Storey')
        title('Shear Force')
        axis([LBa UBa SH(2) max(SH(1:mvt))])
        grid on
        legend('Without-TMD','','With-TMD','')

        %% Plot Displacement and Acceleration on Anim

    case 'anim'
        close all
        figure('Position',[500   100   560   800])
        SH = 3.3*(0:mvt)*mg;
        UB=max([1.2*max(max(abs(Results.Yd*100))) 1.2*max(max(abs(Results.Ytd*100)))]);
        LB = -UB;

        UBa=max([1.2*max(max(abs(Results.Ya./9.807))) 1.2*max(max(abs(Results.Yta./9.807)))]);
        LBa=-UBa;

        for i = 1 : length(t)
            % ------------ Plot Displacement --------------------------
            subplot('Position',[0.13 0.39 0.334 0.534])
            plot(zeros(1,mvt+1),SH,"Color",[0.5 0.5 0.5],'LineStyle','-.','Marker','square','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);
            hold on
            plot([0 Results.Yd(i,:)]*100,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
            plot([0 Results.Ytd(i,:)]*100,SH,'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
            if i>1
                plot([0 max(Results.Yd(1:i,:))]*100,SH(1:mvt),'Color',[1 0.5 0.5],'LineWidth',0.6,'LineStyle','--','Marker','square','MarkerEdgeColor',[1 0.5 0.5],'MarkerFaceColor','none');
                plot([0 max(Results.Ytd(1:i,:))]*100,SH,'Color',[0.5 0.5 1],'LineWidth',0.7,'LineStyle','--','Marker','square','MarkerEdgeColor',[0.5 0.5 1],'MarkerFaceColor','none');
                plot([0 min(Results.Yd(1:i,:))]*100,SH(1:mvt),'Color',[1 0.5 0.5],'LineWidth',0.6,'LineStyle','--','Marker','square','MarkerEdgeColor',[1 0.5 0.5],'MarkerFaceColor','none');
                plot([0 min(Results.Ytd(1:i,:))]*100,SH,'Color',[0.5 0.5 1],'LineWidth',0.7,'LineStyle','--','Marker','square','MarkerEdgeColor',[0.5 0.5 1],'MarkerFaceColor','none');
            end
            axis([LB UB 0 max(SH)])
            yticks(SH)
            yticklabels({0:mvt-1,'TMD'})
            xlabel('Dis _{[Cm]}')
            ylabel('Storey')
            title('Displacement')

            grid on
            % ------------ Plot Acceleration --------------------------
            subplot('Position',[0.57 0.3857 0.3345 0.539])
            plot(zeros(1,mvt+1),SH,"Color",[0.5 0.5 0.5],'LineStyle','-.','Marker','square','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]);
            hold on
            plot([0 Results.Ya(i,:)]./9.807,SH(1:mvt),'Color','r','LineWidth',0.6,'LineStyle','-','Marker','square','MarkerEdgeColor','r','MarkerFaceColor','r');
            plot([0 Results.Yta(i,:)]./9.807,SH,'Color','b','LineWidth',0.7,'LineStyle','-','Marker','square','MarkerEdgeColor','b','MarkerFaceColor','b');
            if i>1
                plot([0 max(Results.Ya(1:i,:))]./9.807,SH(1:mvt),'Color',[1 0.5 0.5],'LineWidth',0.6,'LineStyle','--','Marker','square','MarkerEdgeColor',[1 0.5 0.5],'MarkerFaceColor','none');
                plot([0 max(Results.Yta(1:i,:))]./9.807,SH,'Color',[0.5 0.5 1],'LineWidth',0.7,'LineStyle','--','Marker','square','MarkerEdgeColor',[0.5 0.5 1],'MarkerFaceColor','none');
                plot([0 min(Results.Ya(1:i,:))]./9.807,SH(1:mvt),'Color',[1 0.5 0.5],'LineWidth',0.6,'LineStyle','--','Marker','square','MarkerEdgeColor',[1 0.5 0.5],'MarkerFaceColor','none');
                plot([0 min(Results.Yta(1:i,:))]./9.807,SH,'Color',[0.5 0.5 1],'LineWidth',0.7,'LineStyle','--','Marker','square','MarkerEdgeColor',[0.5 0.5 1],'MarkerFaceColor','none');
            end
            axis([LBa UBa 0 max(SH)])
            yticks(SH)
            yticklabels({0:mvt-1,'TMD'})
            grid on
            xlabel('Acc _{[g]}')
            ylabel('Storey')
            legend('','Without-TMD','With-TMD')
            title('Acceleration')

            subplot('Position',[0.13 0.11 0.78 0.16])
            plot(t,xg)
            hold on
            plot(t(i),xg(i),'o','MarkerFaceColor','red');
            axis([0 T -1.1*max(abs(xg)) 1.1*max(abs(xg))])
            xlabel('Time _{[sec]}')
            ylabel('Acceleration _{[g]}')
            pause(0.0001)
            clf

        end
end
end

