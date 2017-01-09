
# cd to directory where test vcf files for annotation are stored
cd /media/MokaNAS/projects/161208_vlookup/vcf

# make output directory for annotated vcfs
mkdir output
 
filelist=$(ls *gz)

for file in ${filelist[@]}
do
cd /home/kevin/GATK
java -jar GenomeAnalysisTK.jar -R /media/MokaNAS/projects/161012_genome.fa/genome.fa -T VariantAnnotator -o /media/MokaNAS/projects/161208_vlookup/vcf/output/${file%.vcf.gz}_plus_inhouse.vcf.gz --resource:Inhouse /home/kevin/Documents/HGVS/HGVS_VLookup_GenomicCoordinates_RemoveNonFloats_IngenuityCompatible.vcf.gz --resourceAlleleConcordance -V /media/MokaNAS/projects/161208_vlookup/vcf/$file -E Inhouse.PreviousClassification

cd /media/MokaNAS/projects/161208_vlookup/vcf

done