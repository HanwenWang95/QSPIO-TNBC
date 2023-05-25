% Reference: Allen, R J et al. “Efficient Generation and Selection of Virtual Populations
%                               in Quantitative Systems Pharmacology Models.”
%             CPT: pharmacometrics & systems pharmacology vol. 5,3 (2016): 140-6.

function [hist_score,selected,pOut]=OptimizeVPgeneration(p,ProbabilityofInclusion,VPs,runs,observed_variables)

sf=10.^p(1);
threshold=0.1;

dim=numel(VPs(1,:));

pVal=zeros(dim,1);
ksstat=zeros(dim,1);
hist_score=zeros(runs,1);

NumOfVPs=size(VPs,1);
selectedsum=zeros(NumOfVPs,1);

for j=1:runs
    p=rand(NumOfVPs,1);

    selectedtemp=p<ProbabilityofInclusion*sf;

    %selectedsum=selectedsum+double(selectedtemp)/runs;

    selected=selectedtemp>threshold;

    FinalNum=sum(selected);

    if FinalNum>1

        for i=1:dim
            %test if dist. is standard normal
            [~,pVal(i),ksstat(i)]=kstest2(VPs(selected,i),observed_variables(:,i));
        end

        hist_score(j)=sum(ksstat);%max([ksstat1,ksstat3,ksstat4,ksstat5]);%
        pOut=pVal;
    else
        hist_score(j)=1*dim;
    end

end
hist_score=sum(hist_score)/runs;

end
