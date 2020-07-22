function vidStruc=FrameByFrame_Overlay(varargin)
if nargin==2, mVals=varargin{2}; else; mVals=[]; end
vidStruc=varargin{1};
imgAxH = axes('Position',[0 0 1 1],'YDir','normal','units','pixels','YColor','none','XColor','none'); %reverse
% valsAxH = axes('Position',[0 00 1 0.4],'YColor','none','XColor','none','Color', 'none');
% axis on
if ~isempty(mVals) && ~(min(size(mVals))==1)
    if max(abs(mVals(:,3)))>pi; mVals(:,3)=deg2rad(mVals(:,3)); end
end

for frameNum=1:numel(vidStruc)
    axes(imgAxH);
    %     set(imgAxH,'YDir','reverse','YColor','none','XColor','none');
    image(vidStruc(frameNum).cdata); hold on
    
    % plot frame number
    %     text(10,20,num2str(frameNum),'FontSize',30,'color','y')
    %     hold off
    
    %% plot associated values
    if ~isempty(mVals)
        if min(size(mVals))==1 % just a trace
            %         axes(valsAxH);
            %         set(valsAxH,'YColor','none','XColor','none','Color', 'none');
            frameVals=nan(1,300);
            earlyVals=mVals(max([1 frameNum-149]):frameNum-1);
            lateVals=mVals(frameNum:min([frameNum+149 numel(mVals)]));
            frameVals(151-numel(earlyVals):150)=earlyVals;
            frameVals(151:150+numel(lateVals))=lateVals;
            plot(1:300,-(frameVals*2)+560,'color',[121, 215, 111]/255,'linewidth',1.5); %set('xlim',[0 300]);
            %         hold on
            plot(151,-mVals(frameNum)*2+560,'Marker','o','MarkerFaceColor','r');
        elseif min(size(mVals))==2 % trace + spike times
            % assuming phase & spike times
            if exist('paH','var')
                polaraxes(paH); cla;
                set(paH,'visible','on')
            else
                polaraxes('Position',[0.05 0.1 0.3 0.3],'units','pixels','Color','none'); %reverse
            end
            if mVals(2,frameNum)
                
                polarplot(mVals(1,frameNum),1,'Marker','o','MarkerFaceColor','r',...
                    'MarkerEdgeColor','none','Color','none');
            end
            %             if mVals(2,frameNum) || ~exist('paH','var')
            paH = gca; hold on
            polarplot(ones(3,1)*mVals(1,frameNum),0:2,'Color',[121, 215, 111]/255,'linewidth',1.5);
            paH.Color='none';
            paH.ThetaZeroLocation='top';
            paH.ThetaTickLabel={'Protracted','','','\downarrow','','',...
                'Retracted','','','\uparrow','',''};
            paH.ThetaDir = 'counterclockwise';
            paH.FontSize = 15;
            paH.RGrid='off';
            paH.ThetaGrid='on';
            paH.RTickLabel='';
            %             GridLineStyle
            box off
            hold off
            %             end
        else % x/y + angle value
            lineLength=50;
            tipX=mVals(frameNum,1)+(lineLength*cos(mVals(frameNum,3)));
            tipY=mVals(frameNum,2)+(lineLength*sin(mVals(frameNum,3)));
            line([mVals(frameNum,1) tipX],[mVals(frameNum,2) tipY],'color',[121, 215, 111]/255,'linewidth',1.5);
        end
    end
    % save output
    vidStruc(frameNum)=getframe(gcf);
    if exist('paH','var')
        set(paH,'visible','off')
    end
    axes(imgAxH); hold off
end
