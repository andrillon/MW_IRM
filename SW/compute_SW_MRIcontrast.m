slow_Waves(:,5)=slow_Waves(:,5)-(25*500);
slow_Waves(:,8)=slow_Waves(:,8)-(25*500);

%%
box_waves=[];
for nP=1:40
    for nEl=unique(slow_Waves(:,3))'
        temp_waves=slow_Waves(slow_Waves(:,2)==nP & slow_Waves(:,3)==nEl,:);
        temp_vec=zeros(1,10*500);
        for nW=1:size(temp_waves,1)
            temp_vec((temp_waves(nW,5):temp_waves(nW,8))+10*500)=abs(temp_waves(nW,9));
        end
        box_waves(nP,nEl,:)=temp_vec;
    end
end