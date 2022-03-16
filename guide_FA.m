% Code written by Alessandro Del Vecchio and Carina M. Germer
% alessandro.del.vecchio@fau.de 
% pre-print https://www.biorxiv.org/content/10.1101/2022.01.23.477379v1.full.pdf

clear all
close all
clc

% load  Pulses files for the muscles
load('Subj3_VLVM.mat')

% Parameters-------------------
fsamp = 2048; %set your fsamp;
Wind_s = 1;  % hanning window duration 
HanningW = 2/round(fsamp*Wind_s)*hann(round(fsamp*Wind_s)); %unitary area    
Factor = 2; %here the number of factors
DataLength = 60 * fsamp;    %number of samples for each recording
Muscles = {'VL','VM'};

CutStart = 20; % NOTE, SELECT ONLY STATE PART ! 
CutEnd = 4;


ColorS2 = brewermap(10,'Set1');

% choose trial: 1 or 2
Trial = 1; 
MU_Mall{1} = Subj3.VL_Pulses{Trial};
MU_Mall{2} = Subj3.VM_Pulses{Trial};


%discard silent MUs
ActiveMUstart=cellfun(@(y) cellfun(@(z) sum(z<(CutStart+1)*fsamp)>3,y),MU_Mall,'Uni',0);
ActiveMUend=cellfun(@(y) cellfun(@(z) sum(z>DataLength-(CutEnd+1)*fsamp)>3,y),MU_Mall,'Uni',0);
ActiveMU = cellfun(@(w,z) w & z, ActiveMUstart,ActiveMUend,'Uni',0);
MU_M = cellfun(@(w,z) w(z),MU_Mall,ActiveMU,'Uni',0);

%create the discharge trains and the smoothed discharge trains
for m=1:length(MU_M) %for each muscle
    aux = MU_M{m};

    for mu=1:length(aux)
        Firings = zeros(1,DataLength);
        Firings(aux{mu})=1;

        aux_smoothF = filtfilt(HanningW,1,Firings*fsamp);
        SmoothedFirings{m}(mu,:) = detrend(aux_smoothF(CutStart*fsamp+1:end-CutEnd*fsamp),0); %cut initial and final segments and remove mean
    end
end

SFirings_all = [SmoothedFirings{1}; SmoothedFirings{2}]; % here merges all the firings in one matrix. This matrix is the one used for factorization

[Lambda, Psi, rotMat, stats, F] = factoran(SFirings_all',2,'rotat','promax'); % this is the factor analysis, refer to the pre-print for more details

% Separate muscles
NMU_M1=size(SmoothedFirings{1},1);
Lambda_M{1} = Lambda(1:NMU_M1,:);
Lambda_M{2} = Lambda(1+NMU_M1:end,:);

Psi_M{1} = Psi(1:NMU_M1,:);
Psi_M{2} = Psi(1+NMU_M1:end,:);

[~, ModCorrFA] = cellfun(@(x) max(mean(x)),Lambda_M);

% Correlate the smoothed discharges with the two factors
for m=1:2
    firings = SmoothedFirings{m};
    for neuron= 1:size(firings)    
        for Fact = 1:2
            R = corrcoef(firings(neuron,:),detrend(F(:,Fact),0));
            Corrvalue{m}(neuron,Fact) = R(2);
        end
    end
end


figure(), hold on
m=1; h1=biplot(Corrvalue{m},'Color',ColorS2(2+m,:));
m=2; h2=biplot(Corrvalue{m},'Color',ColorS2(2+m,:)); 

grid on
set(gca,'FontSize',10,'FontWeight','bold')
set(gcf,'position',[2 2 150*2 152*2])


% now we need to cluster these values
AllCorrValue = [Corrvalue{1};Corrvalue{2}];
AbsAllCorrValue = abs(AllCorrValue);
[~,idx] = pdist2([.65 0.1; 0.4 0.4; 0.1 .65],AbsAllCorrValue,'euclidean','Smallest',1); %classify with fixed centroid
idxM{1}=idx(1:NMU_M1);
idxM{2}=idx(NMU_M1+1:end);


figure('position',[100 100 300 300]), hold on
m=1; h1=biplot(Corrvalue{m},'Color',ColorS2(2+m,:),'LineWidth',1);
m=2; h2=biplot(Corrvalue{m},'Color',ColorS2(2+m,:),'LineWidth',1); 
scatter(AllCorrValue(idx==1,1),AllCorrValue(idx==1,2),20,ColorS2(1,:),'filled')
scatter(AllCorrValue(idx==2,1),AllCorrValue(idx==2,2),20,ColorS2(9,:),'filled')
scatter(AllCorrValue(idx==3,1),AllCorrValue(idx==3,2),20,ColorS2(2,:),'filled')
grid on
lgndHandles = [h1(1) h2(1)];
lgndLabels = {Muscles{1}, Muscles{2}};
legend(lgndHandles,lgndLabels,'location','southwest')
legend boxoff
xlabel([Muscles{ModCorrFA==1} ' Module'])
ylabel([Muscles{ModCorrFA==2} ' Module'])
set(gca,'FontSize',10,'FontWeight','bold')

%% Total variance

X = SFirings_all';
covariance=(X'*X)/(length(X));

[V,D]=eig(covariance);
[e,i]=sort(diag(D),'descend');
sorted=V(:,i);

e_perc=(e./sum(e))*100;

figure('position',[2 2 150*3 152*2]), hold on
plot(cumsum(e_perc),'ok:','MarkerFaceColor','k')
ylabel('Variance Explained (%)')
xlabel('Factors')
set(gca,'FontSize',12,'FontWeight','bold')
