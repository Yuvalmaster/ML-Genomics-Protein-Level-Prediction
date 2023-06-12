function CpG_feature = CpG(variants,dataset)
sequences = variants.VariantSequence;
precent=[];
for i=1:length(sequences)
    seq=char(sequences(i));
    seq=seq(27:50);
    counter=0;
    for k=1:23
        if seq(k)=='C'
            if seq(k+1)=='G'
                counter=counter+1;
            end
        end
    end
    precent(i)=counter/23;

end
CpG_feature=[];
for i=1:length(table2array(dataset(:,1)))
    indexes=table2array(dataset(i,2:27));
    indexes = indexes(~isnan(indexes));
    CpG_feature(i)=mean(precent(indexes));
end
CpG_feature=transpose(CpG_feature);
