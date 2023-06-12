function [hydrofobic_feature,polar_feature] = hydrofobic_and_polar(variants,dataset)
sequences = variants.VariantSequence;
precent_hydrofobic=[]; precent_polar=[];
for i=1:length(sequences)
    seq=char(sequences(i));
    num_hydrofobic=0; num_polar=0;
    counter=codoncount(seq(27:50));
    num_hydrofobic=num_hydrofobic+counter.GCG+counter.GCC+counter.GCA+counter.GCT+counter.GTT+counter.GTG+counter.GTC...
        +counter.GTA+counter.AAA+counter.AAG+counter.ATT+counter.ATC+counter.TTT+counter.TTC+counter.TGG...
        +counter.CCG+counter.CCT+counter.CCA+counter.CCC+counter.ATA+counter.GGC+counter.GGT+counter.GGA+counter.GGG;
    num_polar=num_polar+counter.TCT+counter.TCC+counter.TCG+counter.TCA+counter.ACT+counter.ACA+counter.ACC...
        +counter.ACG+counter.TGT+counter.TGC+counter.AAT+counter.AAC+counter.CAA+counter.CAG+counter.TAT+counter.TAC;
   precent_hydrofobic(i)=num_hydrofobic/8;
   precent_polar(i)=num_polar/8;

end

hydrofobic_feature=[];
polar_feature=[];
for i=1:length(table2array(dataset(:,1)))
    indexes=table2array(dataset(i,2:27));
    indexes = indexes(~isnan(indexes));
    hydrofobic_feature(i)=mean(precent_hydrofobic(indexes));
    polar_feature(i)=mean(precent_polar(indexes));
end
hydrofobic_feature=transpose(hydrofobic_feature);
polar_feature=transpose(polar_feature);
