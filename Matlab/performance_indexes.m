clear all; clc; close all;

% For the induction phase:
% • TT: observed time-to-target (in seconds) required for reaching
% the first time the target interval of 45–55 BIS values;
% • BIS-NADIR: the lowest observed BIS value;
% • ST10: settling time, defined as the time interval for the BIS to
% reach and steady within the BIS range between 45 and 55 (that
% is, the target value of 50 ± 5);
% • ST20: the same of ST10 but it considers a BIS range of 40 and 60;
% • US: undershoot, defined as the difference between the lower
% threshold of 45 and the minimum value of BIS below this threshold.

% Ffor maintencace phase:
% • TT: observed time-to-target (in seconds) required for reaching
% the first time the target interval of 45–55 BIS values;
% • BIS-NADIR: the lowest observed BIS value;

TT=zeros(13,1);
BIS_NADIR=zeros(13,1);
ST10=zeros(13,1);
ST20=zeros(13,1);
US=zeros(13,1);
TTp=zeros(13,1);
BIS_NADIRp=zeros(13,1);
TTn=zeros(13,1);
BIS_NADIRn=zeros(13,1);

% simulation data following the control scenario from Orignal Paper

Tsim=20*60;
SP=50;
tstep0=0;
astep=10;
tstep1=11*60;
tstep2=16*60;

% load(strcat('data_gpc_dist_2'));
load(strcat('data_gpc_dist_profile_2'));
istep0=1;

for Id=1:13
    istep1=find(t{Id}==tstep1);
    istep2=find(t{Id}==tstep2);
    BIS_NADIR(Id)=min(y{Id}(1:istep1));
    BIS_NADIRp(Id)=min(y{Id}(istep1:istep2-1));
    BIS_NADIRn(Id)=max(y{Id}(istep2:end));
    if BIS_NADIR(Id)<45
        US(Id)=45-BIS_NADIR(Id);
    else
        US(Id)=0;
    end
    flag=0;
    for i=1:1:(istep1-1)
        if y{Id}(i)<=55 & flag==0
            TT(Id)=(t{Id}(i)-t{Id}(istep0))/60;
            flag=1;
        end
        if abs(y{Id}(i)-50)>5
            ST10(Id)=(t{Id}(i+1)-t{Id}(istep0))/60;
        end
        if abs(y{Id}(i)-50)>=10
            ST20(Id)=(t{Id}(i+1)-t{Id}(istep0))/60;
        end
    end
    flag=0;
    for i=istep1:1:(istep2-1)
        if y{Id}(i)<55 & flag==0
            TTp(Id)=(t{Id}(i+1)-t{Id}(istep1))/60;
            flag=1;
        end
    end
    flag=0;
    for i=istep2:1:(length(y{Id}))
        if y{Id}(i)>45 & flag==0
            TTn(Id)=(t{Id}(i+1)-t{Id}(istep2))/60;
            flag=1;
        end
    end

    IAE(Id)=trapz(t{Id},abs(y{Id}-SP));

    figure(1);
    ax1=subplot(2,1,1);
    plot(t{Id}/60,y{Id});
    hold on;
    ax2=subplot(2,1,2);
    plot(t{Id}/60,u{Id}(:,1))
    hold on;
    grid on;

end


figure(1);
subplot(2,1,1);
ylabel('BIS');
xlabel('Time [min]');
grid on;
tlimite=tstep0+5*60;
plot([0 Tsim]/60,[50 50],'k:');
ylim([30 100]);
xlim([0 Tsim]/60);
subplot(2,1,2);
ylabel('Infusion');
xlabel('Time [min]');
grid on;
linkaxes([ax1,ax2],'x')

INDEXES=[TT BIS_NADIR ST20 ST10 US TTp BIS_NADIRp TTn BIS_NADIRn]
IAEmax=max(IAE)

