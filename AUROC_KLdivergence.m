% Analyze number of images needed to discriminate between lesion and
% control monkey
load('C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\KLdivergence-fixation-location_List.mat')

healthyKL = KLnorm(:,:,1);
lesionKL = KLnorm(:,:,2);
lesionKL(isnan(lesionKL(:,1)),:)=[];
ROC = cell(1,size(KLnorm,2));

thresh = 0:100;
setlength = 30;%not actual set size

titles = {
    'Fixations 1-5';
    'Fixations 6-10';
    'Fixations 11-15';
    'Fixations 16-20';
    'Fixations 21-25';
    'All Fixations 1-25';
    };


figure
for group = 1:6
    h = subplot(2,3,group);
    colorset = cool;
    set(gca,'ColorOrder',colorset(1:2:end,:))
    hold all
    for i = 1:19;
        subsetind = 1:30*i;
        TP = NaN(1,length(thresh)); %True positive
        FA = NaN(1,length(thresh)); %False alarm
        health = healthyKL(subsetind,group);
        health(isnan(health)) = [];
        lesion = lesionKL(subsetind,group);
        lesion(isnan(lesion)) = [];
        
        for ii = 1:length(thresh)
            TP(ii) = sum(health > thresh(ii))/length(health);
            FA(ii) = sum(lesion > thresh(ii))/length(lesion);
        end
        ROC{group}(i) = -trapz(FA,TP);
        plot(FA(1,:),TP(1,:))
    end
    plot([0 1],[0 1],'k--','linewidth',3)
    hold off
    xlabel('FA')
    ylabel('TP')
    xlim([0 1])
    ylim([0 1])
    title(titles{group})
end

figure
for group = 1:6
    subplot(2,3,group)
    plot(ROC{group})
    xlabel('Number of Sets of 30 images')
    ylabel('AUROC a.u.')
    title(titles{group})
end