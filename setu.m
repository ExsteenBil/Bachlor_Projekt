function [TA,hTid,Localtid,Yf,dag]=setu(F1,F2,Tstart,Tslut,n);
%F1: Filnavn af forsøg fra datalogger 1
%F2: Filnavn af forsøg fra datalogger 2
%n: Hvorvidt der skal benyttes et specifikt tidsintervalg
%Tstart: Start tidspunktet for det specifikt tidsintervalg
%Tslut: slut tidspunktet for det specifikt tidsintervalg

%TA: Datalogger 1 og 2 samlet i et dataset
%htid: tiden i timer
%Localtid: Datetime format at hele forsøgset
%Yf: året
%dag: antal dage siden 1 Jan
%Written by: Ejnar Exsteen & Peter Garlet Rank - S193842 & s193851  2020
%date: 01/06/2022
T1 = readtable(F1,'PreserveVariableNames',true);    
T2 = readtable(F2,'PreserveVariableNames',true);
%fjerner nan values
T1= rmmissing(T1);
T2= rmmissing(T2);
%Navn for mit table
NameList1 = ["Time","C1","C2","C3","C4","C5","C6","C7","C8","C9","Solar Gain","Event"];
NameList2 = ["Time","FA","AAT","AAB","FB","ABT","ABB","P1","U1","Event"];
%Giver dataloggerns table nyt navn
allVars1 = 1:length(NameList1);
newNames1 = append(NameList1);
T1 = renamevars(T1,allVars1,newNames1);
allVars2 = 1:length(NameList2);
newNames2 = append(NameList2);
T2 = renamevars(T2,allVars2,newNames2);

%Matcher størrelsen på de 2 dataset, således at den mindste bestemmer
%størelsen
if size(T2,1) < size(T1,1)
    T1 = T1(1:size(T2,1),:);
    
else if size(T1,1) < size(T2,1)
    T2 = T2(1:size(T1,1),:);
else
   disp("Match of Dataset")
end
end
% fjerner "event" fra dataset
T1 = removevars(T1,size(T1,2));
T2 = removevars(T2,size(T2,2));
%fjerner tid fra Datalogger 1
T1 = removevars(T1,1);
%Samle dataloggerne i et dataset
TA = [T1,T2];
%Glidende gennemsnit
%TA.AAT = movmean(TA.AAT,12*5);
%TA.AAB = movmean(TA.AAT,12*5);
%TA.ABT = movmean(TA.ABT,12*5);
%TA.ABB = movmean(TA.ABB,12*5);
%TA.FB = movmean(TA.FB,12*5);
%TA.FA = movmean(TA.FA,12*5);
%TA.C1 = movmean(TA.C1,12*5);
%TA.C2 = movmean(TA.C2,12*5);
%TA.C3 = movmean(TA.C3,12*5);
%TA.C4 = movmean(TA.C4,12*5);
%TA.C5 = movmean(TA.C5,12*5);
%TA.C6 = movmean(TA.C6,12*5);
%TA.C7 = movmean(TA.C7,12*5);
%TA.C8 = movmean(TA.C8,12*5);
%TA.C9 = movmean(TA.C9,12*5);
%TA.P1 = movmean(TA.P1,12*5);

%Format af forsøgsdags dato/tid (bemærk det er sluttidspunkt)
FileInfo = dir(F1);
[Yf, Mf, dag, Hf, MNf, Sf] = datevec(FileInfo.datenum);
Date = strjoin({sprintf('%.0f',Yf),sprintf('%.0f',Mf),sprintf('%.0f',dag)},'-');
Tiden = strjoin({sprintf('%.0f',Hf),sprintf('%.0f',MNf),sprintf('%.0f',Sf)},':');
Datetid =strjoin({Date,Tiden},'\t');
tid1 = datetime(Datetid,'InputFormat','yyyy-MM-dd HH:mm:ss');
%Korrigeres for med længde af forsøget for at finde starten af forsøget
Localtid = flip(tid1 - seconds(TA.Time));
[h,m,s] =hms(Localtid);
%tid i timer
hTid = h + m/60 + s/3600;

%Valg af tidsintervalg
if n==1
  %intersect mellem valgte tidsintervalg og dataset
[hs,ms,ss]  =hms(datetime(Tstart));
[ht,mt,st]  =hms(datetime(Tslut));
[val1,pos1] = intersect(find(hs==h),find(ms==m));
[val2,pos2] = intersect(find(ht==h),find(mt==m));
Localtid = Localtid(val1(1):val2(1),:);
[h,m,s] =hms(Localtid);
%tid i timer
hTid = h + m/60 + s/3600;
TA =TA(val1(1):val2(1),:);
%Antal dag det er i året
dag = day(tid1,"dayofyear");
end
