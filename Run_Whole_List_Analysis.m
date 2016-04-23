%% Runs whole List task analysis
% written by Seth Konig 2013
% apparently this code got deleted by accident, realized on 9/11/13 but
% data for KL values still exists. Re wrote this section of code on June 21, 2014

% [1] Get the salience maps for multiple task sets
% [1.5] Get Average Saliency Map
% [2] Get Fixations and Saccades from Behavioral Files
% [3] Calculate Salience at Each Fixation
% [4] Calculate Average Salience at Each Fixation Across mutliple data sets
% [5] Calculate Similarity in Fixation Locations with KL divergence
% [6] Calculate Fixation Durations and Saccade Amplitudes across image sets
%%
%---[1] Get the salience maps for multiple task sets---%
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
for imset = 1:12;
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:length(fileList)
        bmps = strfind(fileList{i},'bmp');
        if ~isempty(bmps)
            if double(fileList{i}(bmps-2)) <= 57 %ascii for number
                imageindexes = [imageindexes i];
            end
        end
    end
    for i = 1:length(imageindexes)
        imagefile = fileList{imageindexes(i)};
        getSalienceMap(imagefile)
    end
end
%%
%---[1.5] Get Average Saliency Map---%%
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
imageX = 800;
imageY = 600;
avgmap = zeros(imageY,imageX);
gray = NaN(2,1);
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:90
        load([num2str(i) '-saliencemap'],'fullmap');
        if any(any(isnan(fullmap)))
            disp('nans :(')
        end
        avgmap = avgmap+fullmap;
    end
end
imagesc(avgmap)
%%
%---[2] Get Fixations and Saccades from Behavioral Files---%
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
imageX = 800;
imageY = 600;
CNDFile = [list_image_dir 'List.cnd'];
ITMFile = [list_image_dir 'List01.itm']; %same item file structure for all files in this task only directory for images changes
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    dirData = dir(dirName);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    eyeindexes = [];
    for i = 1:length(fileList)
        period = strfind(fileList{i},'.');
        if (length(fileList{i}) == period(end)+1) && ...
                (double(fileList{i}(period+1)) <=57)
            eyeindexes = [eyeindexes i];
        end
    end
    for i = 1:length(eyeindexes)
        cortexfile = fileList{eyeindexes(i)};
        getListeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY)
    end
end
%%
%---[3] Calculate Salience at Each Fixation---%
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
PLOTOPTIONS = 'none';
imageX = 800; imageY = 600;
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    for eyefile = eyedatafiles;
        fixationSalience_and_significanceList(matfiles.mat{eyefile},imageX,imageY,PLOTOPTIONS)
    end
end

%%
%%---[4] Calculate Average Salience at Each Fixation Across mutliple data sets---%
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
tags = {'MP','TT'};
minlen = 100;
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    matfiles = what;
    statfiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'FixationStatistics.mat');
        if ~isempty(str)
            statfiles = [statfiles i];
        end
    end
    for stat = statfiles;
        load(matfiles.mat{stat},'statistics')
        minlen = min([minlen size(statistics.numbervalues,3)]);
    end
end

allimages = NaN(180*2*12,1);
alldata = NaN(180*2*12,3,2,minlen);
allshuffled = NaN(180*2*12,3,2,minlen);
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    matfiles = what;
    statfiles = zeros(1,length(tags));
    for i = 1:length(matfiles.mat);
        if ~isempty(strfind(matfiles.mat{i},'FixationStatistics.mat'));
            for ii = 1:length(tags);
                if ~isempty(strfind(matfiles.mat{i},tags{ii}))
                    statfiles(ii) = i;
                end
            end
        end
    end
    
    for stat = 1:length(statfiles);
        if statfiles(stat) ~= 0
            i = length(tags)*(imset-1)+stat;
            load(matfiles.mat{statfiles(stat)})
            combineddata = shuffunshuffdata{2}{3};
            combinedshuffled = shuffunshuffdata{1}{3};
            virtind = 180*(i-1)+1:i*180;
            virtind = virtind(1:size(combineddata,1));
            allimages(virtind) = images;
            alldata(virtind,:,:,:) = combineddata(:,:,:,1:minlen);
            allshuffled(virtind,:,:,:) = combinedshuffled(:,:,:,1:minlen);
        end
    end
end

allmeanvals = NaN(3,2,minlen);
allstdvals = NaN(3,2,minlen);
allnumvals = NaN(3,2,minlen);
for i = 1:size(allmeanvals,1)
    for ii = 1:size(allmeanvals,2)
        for iii = 1:minlen
            allmeanvals(i,ii,iii) = nanmean(alldata(:,i,ii,iii));
            allstdvals(i,ii,iii) = nanstd(alldata(:,i,ii,iii));
            allnumvals(i,ii,iii) = sum(~isnan(alldata(:,i,ii,iii)));
        end
    end
end

% z-test of means agains random distributions assuming mean is larger
allzp = NaN(size(allmeanvals)); %p-values
allcI = NaN(size(allmeanvals)); %top confidence interval value, lowest is typticall 0/-Inf
for i = 1:size(allmeanvals,1)
    for ii = 1:size(allmeanvals,2)
        shuffledvals = allshuffled(:,i,ii,:);
        shuffledvals(isnan(shuffledvals)) = [];
        for iii = 1:size(allmeanvals,3)
            [~,p,ci] = ztest(shuffledvals,allmeanvals(i,ii,iii),std(shuffledvals),...
                0.05);
            allzp(i,ii,iii) = p;
            if i == 3;
                allcI(i,ii,iii) = ci(1);
            else
                allcI(i,ii,iii) = ci(2);
            end
        end
    end
end

allstatistics.meanvalues = allmeanvals;
allstatistics.stdvalues = allstdvals;
allstatistics.numbervalues = allnumvals;
allstatistics.pvalues = allzp;
allstatistics.confidenceintervals = allcI;

clrs = ['rb'];
legendlabels = {'Salience','Salience Contrast','Image Intensity'};
averaginglabels = {' at mean fixation location',' mean during fixation'};
for i = 1:size(allcI,1)
    figure
    hold on
    for ii = 1:size(allcI,2)
        p(ii)= plot(allcI(i,ii)*ones(1,size(allcI,3)),['--' clrs(ii)]);
        pp(ii) = errorbar(shiftdim(allmeanvals(i,ii,:)),...
            shiftdim(allstdvals(i,ii,:))./sqrt(shiftdim(allnumvals(i,ii,:))),clrs(ii));
    end
    hold off
    title([legendlabels{i} ' at a fixation by fixation number across all'...
        ' all monkeys and image sets'])
    legend([p pp],{['Chance ' legendlabels{i} averaginglabels{1}],...
        ['Chance ' legendlabels{i} averaginglabels{2}],...
        [legendlabels{i} averaginglabels{1}], ...
        [legendlabels{i} averaginglabels{2}]},'Location','NorthEastOutside')
    xlabel('Fixation Number')
    ylabel('Normalized Value')
end

allstatvariablenames = {
    'alldata: cell arrray containing combined values for salience, salience';
    'contrast, and image intensity at fixations across all monkeys and data files.';
    'Each of these cells are arranged by row,column,and z-column in the following';
    'manner:  Rows are arranged as Salience, salience contrast, and intensity;';
    'Columns indicate if data is these parmaters at the average fixation cooridante';
    '(col 1) or the average of the parameters during a fixatoin; Z-column is';
    'organized by fixation number';
    '';
    'allshuffled: same as alldata but shuffled data';
    '';
    'allstatistics: structure with results from a z-test to determine if fixations';
    'occur at salience, salience contrasts, and image intesntisy valeus at rates';...
    'higher than what would be expected by chance. Random distributions from each';
    'parameter is compared to the mean value at that parameter by fixation';
    'statistics contains means, std, number of fixations, p-values,and confidence intervals';
    'These variables are arranged by row,column,z-column. Rows are arranged as';
    'Salience, salience contrast, and intensity. Columns indicate if data is these';
    'paramaters at the average fixation cooridante (col 1) or the average of the';
    'parameters during a fixatoin. Z-column is organized by fixation number';
    };

save(['C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\'...
    'Combined-List-Salience'],'allimages','minlen',...
    'alldata','allshuffled','allstatistics','allstatvariablenames');

novsal = NaN(90*12,minlen);
repsal = NaN(90*12,minlen);
novsalTT = NaN(90*12,minlen);
repsalTT = NaN(90*12,minlen);
for imset = 1:12
    for STAT = 1:2
        i = length(tags)*(imset-1)+STAT;
        virtind = 180*(i-1)+1:i*180;
        tempsal = alldata(virtind,1,1,:);
        tempsal = reshape(tempsal,[size(tempsal,1),size(tempsal,4)]);
        images = allimages(virtind);
        nanind = find(isnan(images));
        images(nanind) = [];
        tempsal(nanind,:) = [];
        if ~isempty(images);
            if STAT == 1;
                for img = 1:90
                    ind = find(images == img);
                    if ~isempty(ind)
                        novsal(90*(imset-1)+img,:) = tempsal(ind(1),:);
                    end
                    if length(ind) == 2;
                        repsal(90*(imset-1)+img,:) = tempsal(ind(2),:);
                    end
                end
            else
                for img = 1:90
                    ind = find(images == img);
                    if ~isempty(ind)
                        novsalTT(90*(imset-1)+img,:) = tempsal(ind(1),:);
                    end
                    if length(ind) == 2;
                        repsalTT(90*(imset-1)+img,:) = tempsal(ind(2),:);
                    end
                end
                
            end
        end
    end
end
nans = find(isnan(novsal(:,1)));
novsal(nans,:) = [];
nans = find(isnan(repsal(:,1)));
repsal(nans,:) = [];
nans = find(isnan(novsalTT(:,1)));
novsalTT(nans,:) = [];
nans = find(isnan(repsalTT(:,1)));
repsalTT(nans,:) = [];

npnov = sqrt(size(novsal,1));
nprep = sqrt(size(repsal,1));
npnovTT = sqrt(size(novsalTT,1));
nprepTT = sqrt(size(repsalTT,1));
figure
hold on
errorbar(nanmean(novsal),nanstd(novsal)/npnov,'b')
errorbar(nanmean(repsal),nanstd(repsal)/nprep,'r')
errorbar(nanmean(novsalTT),nanstd(novsalTT)/npnovTT,'g')
errorbar(nanmean(repsalTT),nanstd(repsalTT)/nprepTT,'m')
hold off
legend('MP novel','MP repeated','TT novel','TT repeated')
ylabel('Salience')
xlabel('Fixation number')
box off
%%
%---[5] Calculate Similarity in Fixation Locations with KL divergence ---%
list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
tags = {'MP','TT'};
imageX = 800;
imageY = 600;

%row: by image
%column 1: fixations 1-5
%column 2: fixations 6-10
%column 3: fixations 11-15
%column 4: fixations 16-20
%column 5: fixations 21-25, most trials have at least 20 fixations so might have little data here
%column 6: all fixations from 1 up to 25
% zcolom is by monkey

%pre-allocate space for data. Must be at least as big as necessary (the number of images
% the monkey ran on) or get 0's that screw up the data.
KLshuff = NaN(1080,6,2); %shuffled data-so chance levels
KLnorm = NaN(1080,6,2); % observed data or normal data
count = [1 1]; %keep track of image number across sets

for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    
    for eye = 1:length(eyedatafiles)
        monkey_ind = find(~cellfun(@isempty,strfind(tags,matfiles.mat{eyedatafiles(eye)}(1:2)))); %which monkey
        load(matfiles.mat{eyedatafiles(eye)}) %load the fixation location data
        for img = 1:90;
            img_ind = find(images == img); %find image pairs (novel and repeat)
            if length(img_ind) == 2; %if both images were presented.
                %For this task there was an error in the block structure so it's possible that an image was displayed 1 or 3 times
                
                %these will contain PDFs (probability distribution functions) of fixation locations
                novelfixations = cell(1,6);
                repeatfixations = cell(1,6);
                shuffled_novelfixations = cell(1,6);
                shuffled_repeatfixations = cell(1,6);
                for i = 1:length(novelfixations); %fill with a matrix of zeros
                    novelfixations{i} = zeros(imageY,imageX);
                    repeatfixations{i}=zeros(imageY,imageX);
                    shuffled_novelfixations{i} = zeros(imageY,imageX);
                    shuffled_repeatfixations{i} = zeros(imageY,imageX);
                end
                
                nov_fixations = fixationstats{img_ind(1)}.fixations; %fixation locations from novel presentation
                rep_fixations = fixationstats{img_ind(2)}.fixations; %fixation locations from repeat presentation
                
                if ~isempty(nov_fixations) && ~isempty(rep_fixations)
                    %remove 1st fixation if this is a fixation on the cross-hair
                    if nov_fixations(1,1) > imageX/2-100 && nov_fixations(1,1) < imageX/2+100 &&...
                            nov_fixations(2,1) < imageY/2+100 && nov_fixations(2,1) > imageY/2-100
                        nov_fixations(:,1) = [];
                    end
                    nov_fixations = round(nov_fixations);
                    %remove 1st fixation if this is a fixation on the cross-hair
                    if rep_fixations(1,1) > imageX/2-100 && rep_fixations(1,1) < imageX/2+100 &&...
                            rep_fixations(2,1) < imageY/2+100 && rep_fixations(2,1) > imageY/2-100
                        rep_fixations(:,1) = [];
                    end
                    rep_fixations=round(rep_fixations);
                    
                    %we want to take the same number of fixations from the novel and repeat trials
                    maxfixations = min(size(nov_fixations,2),size(rep_fixations,2));
                    maxfixations(maxfixations > 25) = 25; %don't care if there are more than 25 fixations
                    
                    %since puting into groups of five remove the remainder of #/5
                    maxfixations = maxfixations-rem(maxfixations,5);
                    
                    if maxfixations >=5
                        for fixation = 1:maxfixations;
                            nov_fix_x = nov_fixations(1,fixation);%horizontal fixation position for novel presentation
                            nov_fix_y = nov_fixations(2,fixation);%vertical fixation position for novel presentation
                            rep_fix_x = rep_fixations(1,fixation);%horizonal fixation position for repeat presentation
                            rep_fix_y = rep_fixations(2,fixation);%vertical fixation position for repeat presentation
                            shuff_nov_fix_x = randi(imageX);
                            shuff_nov_fix_y = randi(imageY);
                            shuff_rep_fix_x = randi(imageX);
                            shuff_rep_fix_y = randi(imageY);
                            
                            %make sure fixations are within image borders
                            nov_fix_x(nov_fix_x < 1) = 1;
                            nov_fix_x(nov_fix_x > imageX) = imageX;
                            nov_fix_y(nov_fix_y < 1) = 1;
                            nov_fix_y(nov_fix_y > imageY) = imageY;
                            rep_fix_x(rep_fix_x < 1) = 1;
                            rep_fix_x(rep_fix_x > imageX) = imageX;
                            rep_fix_y(rep_fix_y < 1) = 1;
                            rep_fix_y(rep_fix_y > imageY) = imageY;
                            %shuffled may be 0 but not greater than imageX or imageY
                            shuff_nov_fix_x(shuff_nov_fix_x < 1) = 1;
                            shuff_nov_fix_y(shuff_nov_fix_y < 1) = 1;
                            shuff_rep_fix_x(shuff_rep_fix_x < 1) = 1;
                            shuff_rep_fix_y(shuff_rep_fix_y < 1) = 1;
                            
                            %put fixations in their appropriate PDF. Mark matrix with a
                            %1 where there was a fixation
                            if fixation <= 5
                                novelfixations{1}(nov_fix_y,nov_fix_x) = novelfixations{1}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{1}(rep_fix_y,nov_fix_x) = repeatfixations{1}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{1}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{1}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <= 10
                                novelfixations{2}(nov_fix_y,nov_fix_x) = novelfixations{2}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{2}(rep_fix_y,nov_fix_x) = repeatfixations{2}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{2}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{2}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <=15
                                novelfixations{3}(nov_fix_y,nov_fix_x) = novelfixations{3}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{3}(rep_fix_y,nov_fix_x) = repeatfixations{3}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{3}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{3}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <=20
                                novelfixations{4}(nov_fix_y,nov_fix_x) = novelfixations{4}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{4}(rep_fix_y,nov_fix_x) = repeatfixations{4}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{4}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{4}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            elseif fixation <=25
                                novelfixations{5}(nov_fix_y,nov_fix_x) = novelfixations{5}(nov_fix_y,nov_fix_x)+1;
                                repeatfixations{5}(rep_fix_y,nov_fix_x) = repeatfixations{5}(rep_fix_y,nov_fix_x)+1;
                                shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{5}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                                shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{5}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                            end
                            %for all fixations
                            novelfixations{6}(nov_fix_y,nov_fix_x) = novelfixations{6}(nov_fix_y,nov_fix_x)+1;
                            repeatfixations{6}(rep_fix_y,nov_fix_x) = repeatfixations{6}(rep_fix_y,nov_fix_x)+1;
                            shuffled_novelfixations{6}(shuff_nov_fix_y,shuff_nov_fix_x) = shuffled_novelfixations{6}(shuff_nov_fix_y,shuff_nov_fix_x)+1;
                            shuffled_repeatfixations{6}(shuff_rep_fix_y,shuff_nov_fix_x) = shuffled_repeatfixations{6}(shuff_rep_fix_y,shuff_nov_fix_x)+1;
                        end
                        
                        % Function below...
                        % 1. convoles fixation matrix with 24 pixel (~1dva) gaussian filter
                        % to account for variability in fixation locations and eye tracking error
                        % 2. Bins (the verb) convolved matrix
                        % 3. removes 0's and replaces with minimum defined value in matlab eps (2 ^-52)
                        % 4. creates fixation pdf by dividing matrix by the sum of the matrix
                        % 5. calculate KL-divergence (symmetric form)
                        for i = 1: maxfixations/5;
                            Distance = KL_Divergence(novelfixations{i},repeatfixations{i});
                            Shuffled_Distance = KL_Divergence(shuffled_novelfixations{i},shuffled_repeatfixations{i});
                            KLshuff(count(monkey_ind),i,monkey_ind) = Shuffled_Distance; %shuffled data-so chance levels
                            KLnorm(count(monkey_ind),i,monkey_ind) = Distance; % observed data or normal data
                        end
                        % for all fixations
                        Distance = KL_Divergence(novelfixations{6},repeatfixations{6});
                        Shuffled_Distance = KL_Divergence(shuffled_novelfixations{6},shuffled_repeatfixations{6});
                        KLshuff(count(monkey_ind),6,monkey_ind) = Shuffled_Distance; %shuffled data-so chance levels
                        KLnorm(count(monkey_ind),6,monkey_ind) = Distance; % observed data or normal data
                    end
                    count(monkey_ind) = count(monkey_ind)+1;
                end
            end
        end
    end
end
save('C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\KLdivergenceData','KLnorm','KLshuff')

chance = KLshuff(:,1:5,:);
chance = nanmean(chance(1:end));
figure
hold on
bar(nanmean(KLnorm(:,:,1)),'b')
bar(nanmean(KLnorm(:,:,2)),'r')
plot(0.5:5.5,chance*ones(1,6),'k--')
plot(5.5:6.5,nanmean(nanmean(KLshuff(:,6,:)))*ones(1,2),'k--')
hold off
title('Similiarity Between Novel and Repeat Fixation Locations')
ylabel('KL Divergence (Distance in Bits')
set(gca,'XTick',1:6);
set(gca,'XTickLabel',{'1-5','6-10','11-15','16-20','21-25','all 1-25'});
xlabel('Fixations')

% test if observed data is similar to chance
pvals_chance = NaN(2,6);
for i = 1:2
    for ii = 1:6
        chance = KLshuff(:,ii,i);
        chance = chance(1:end);
        data = KLnorm(:,ii,i);
        data = data(1:end);
        [~,p] = ttest2(data,chance);
        pvals_chance(i,ii) = p;
    end
end

% test if MP has different behavior than TT who has a hippocampal lesion
pvals_lesion = NaN(1,6);
for ii = 1:6
    mp = KLnorm(:,ii,1);
    mp = mp(1:end);
    tt = KLnorm(:,ii,2);
    tt = tt(1:end);
    [~,p] = ttest2(mp,tt);
    pvals_lesion(ii) = p;
end
%%
%--- [6] Calculate Fixation Durations and Saccade Amplitudes across image sets---%

%--- [6] Calculate Fixation Durations and Saccade Amplitudes across image sets---%

list_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\List Task Analysis\List Image Sets\';
tags = {'MP','TT'};
imageX = 800;
imageY = 600;

fixation_durations = cell(2,length(tags));% row novel vs repeat, col by monkey
sac_amps = cell(2,length(tags));% row novel vs repeat, col by monkey
for i = 1:size(fixation_durations,1);
    for ii = 1:size(fixation_durations,2)
        fixation_durations{i,ii} = NaN(1080,75);
        sac_amps{i,ii} = NaN(1080,75);
    end
end

count = ones(1,length(tags));
for imset = 1:12
    if imset < 10
        dirName = [list_image_dir 'List' '0' num2str(imset)];
    else
        dirName = [list_image_dir 'List' num2str(imset)];
    end
    cd(dirName)
    matfiles = what;
    eyedatafiles = [];
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            eyedatafiles = [eyedatafiles i];
        end
    end
    
    for eye = 1:length(eyedatafiles)
        monkey_ind = find(~cellfun(@isempty,strfind(tags,matfiles.mat{eyedatafiles(eye)}(1:2)))); %which monkey
        load(matfiles.mat{eyedatafiles(eye)}) %load the fixation location data
        for img = 1:90;
            img_ind = find(images == img); %find image pairs (novel and repeat)
            if length(img_ind) == 2; %if both images were presented.
                %For this task there was an error in the block structure so it's possible that an image was displayed 1 or 3 times
                
                %these will contain PDFs (probability distribution functions) of fixation locations
                
                nov_fixations = fixationstats{img_ind(1)}.fixations; %fixation locations from novel presentation
                rep_fixations = fixationstats{img_ind(2)}.fixations; %fixation locations from repeat presentation
                
                nov_fixationtimes = fixationstats{img_ind(1)}.fixationtimes; %fixation locations from novel presentation
                rep_fixationtimes = fixationstats{img_ind(2)}.fixationtimes; %fixation locations from repeat presentation
                nov_saccadetimes = fixationstats{img_ind(1)}.saccadetimes;
                rep_saccadetimes = fixationstats{img_ind(2)}.saccadetimes;
                nov_XY = fixationstats{img_ind(1)}.XY;
                rep_XY = fixationstats{img_ind(2)}.XY;
                
                if ~isempty(nov_fixations) && ~isempty(rep_fixations)
                    %remove 1st fixation if this is a fixation on the cross-hair
                    if nov_fixations(1,1) > imageX/2-100 && nov_fixations(1,1) < imageX/2+100 &&...
                            nov_fixations(2,1) < imageY/2+100 && nov_fixations(2,1) > imageY/2-100
                        nov_fixations(:,1) = [];
                        nov_fixationtimes(:,1) = [];
                    end
                    nov_fixations = round(nov_fixations);
                    %remove 1st fixation if this is a fixation on the cross-hair
                    if rep_fixations(1,1) > imageX/2-100 && rep_fixations(1,1) < imageX/2+100 &&...
                            rep_fixations(2,1) < imageY/2+100 && rep_fixations(2,1) > imageY/2-100
                        rep_fixations(:,1) = [];
                    end
                    rep_fixations=round(rep_fixations);
                    rep_fixationtimes(:,1) = [];
                end
                if ~isempty(nov_fixationtimes)
                    fixdur = 5*(nov_fixationtimes(2,:)-nov_fixationtimes(1,:))+5;
                    fixation_durations{1,monkey_ind}(count(monkey_ind),1:length(fixdur)) = fixdur;
                end
                
                if ~isempty(rep_fixationtimes)
                    fixdur = 5*(rep_fixationtimes(2,:)-rep_fixationtimes(1,:))+5;
                    fixation_durations{2,monkey_ind}(count(monkey_ind),1:length(fixdur)) = fixdur;
                end
                
                if ~isempty(nov_saccadetimes)
                    sac = [];
                    for s = 1:size(nov_saccadetimes,2);
                        sacx = nov_XY(1,nov_saccadetimes(2,s))-nov_XY(1,nov_saccadetimes(1,s));
                        sacy = nov_XY(2,nov_saccadetimes(2,s))-nov_XY(2,nov_saccadetimes(1,s));
                        sac(s) = sqrt(sacx^2+sacy^2);
                    end
                    sac_amps{1,monkey_ind}(count(monkey_ind),1:length(sac)) = sac/24;
                end
                
                if ~isempty(rep_saccadetimes)
                    sac = [];
                    for s = 1:size(rep_saccadetimes,2);
                        sacx = rep_XY(1,rep_saccadetimes(2,s))-rep_XY(1,rep_saccadetimes(1,s));
                        sacy = rep_XY(2,rep_saccadetimes(2,s))-rep_XY(2,rep_saccadetimes(1,s));
                        sac(s) = sqrt(sacx^2+sacy^2);
                    end
                    sac_amps{2,monkey_ind}(count(monkey_ind),1:length(sac)) = sac/24;
                end
                
                count(monkey_ind) = count(monkey_ind)+1;
            end
        end
    end
end

median_fix = NaN(size(fixation_durations));
median_sac = NaN(size(sac_amps));
for i = 1:size(fixation_durations,1);
    for ii = 1:size(fixation_durations,2)
        countf = sum(~isnan(fixation_durations{i,ii}'));
        countf(countf == 0) = [];
        median_fix(i,ii) = median(countf);
        
        counts = sum(~isnan(sac_amps{i,ii}'));
        counts(counts == 0) =[];
        median_sac(i,ii) = median(counts);
    end
end

clr = ['rg';'bm'];
figure
hold on
for ii = 1:size(fixation_durations,2)
    for i = 1:size(fixation_durations,1);
        errorbar(nanmean(fixation_durations{i,ii}(:,1:median_fix(i,ii))),...
            nanstd(fixation_durations{i,ii}(:,1:median_fix(i,ii)))...
            ./sqrt(sum(~isnan(fixation_durations{i,ii}(:,1:median_fix(i,ii))))),clr(i,ii));
    end
end
legend('MP Nov','MP Repeat','TT nov','TT Repeat')

%%
clr = ['rg';'bm'];
figure
hold on
for ii = 1:size(sac_amps,2)
    for i = 1:size(sac_amps,1);
        errorbar(nanmean(sac_amps{i,ii}(:,1:median_fix(i,ii))),...
            nanstd(sac_amps{i,ii}(:,1:median_fix(i,ii)))...
            ./sqrt(sum(~isnan(sac_amps{i,ii}(:,1:median_fix(i,ii))))),clr(i,ii));
    end
end
legend('MP Nov','MP Repeat','TT nov','TT Repeat')