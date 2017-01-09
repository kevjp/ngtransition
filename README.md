# ngtransition


You start with
1 - a list of classified variants from the scientists
2 - a list of gene symbols and the relevant transcript, one transcript per gene symbol (with panel name also specified but not used)

These two inputs are combined to generate `HGVS_VLookup.csv`, which is then put through Mutalyzer (via website)

The Mutalyzer output `Recognised_by_mutalyser.csv` is then augmented with classifications from `HGVS_VLookup.csv`

This is then converted into a VCF file

sed is used to remove any classifications that aren't Class1, Class2, Class3, Class4 or Class5
it also gets rid of the `Class` so that the string becomes an integer as required by IVA for filtering

Finally GATK VariantAnnotator is used to annotate VCFs with this classification information.
