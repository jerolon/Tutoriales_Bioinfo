# How to design in situ probes:
## First step, decide on a sequence or set of sequences and extract them from the transcriptome file
Say you have decided on the Transcript that is identified in our assembly by TRINITY_DN28_c4_g2_i16, or maybe you want all the isoforms of that gene.
Navigate to folder 29-Stranded_Aseembly/9cdhitest, or wherever the transcriptome fasta file is located. 
Load the module seqtk/1.2-r94
```bash
module load seqtk/1.2-r94
#Grep the StrandedTrinityAss_95.fasta transcriptome file for the lines containing the IDs:
#This will get youu all isoforms
grep TRINITY_DN28_c4_g2 StrandedTrinityAss_95.fasta
#This would get you a specific one
grep TRINITY_DN28_c4_g2_i16 StrandedTrinityAss_95.fasta
```
The output will look something like this. We have to get rid of the `>´ thingies, though
```bash
>TRINITY_DN28_c4_g2_i16 len=4059 path=[0:0-286 1:287-288 2:289-556 3:557-712 4:713-2292 5:2293-2359 7:2360-3254 8:3255-3303 10:3304-4058]
>TRINITY_DN28_c4_g2_i3 len=4032 path=[0:0-286 1:287-288 2:289-556 3:557-712 4:713-2292 6:2293-2332 7:2333-3227 9:3228-3276 10:3277-4031]
```

For this, we will pipe the output into sed, then redirect into a file containing the IDs of the desired genes
```bash
grep TRINITY_DN28_c4_g2_i16 StrandedTrinityAss_95.fasta  |
sed 's/>//' > fileofIDs.txt
```

Make sure you loaded the seqtk module. With this program, we can extract the fasta sequences whose ID in the fileofIDs.txt:
```bash
seqtk subseq StrandedTrinityAss_95.fasta fileofIDs.txt > mytranscripts.fasta
```
Use this sequences to get a list of candidate primers in [Primer Blast](https://www.ncbi.nlm.nih.gov/tools/primer-blast/) that give products of the desired length, etc.. And come back when you get those

## How to check for primer specificity with in silico PCR
Once you get a list of candidate primers, we will check that they will only give us hits to the intended targets.
Load the module with the UCSC executables and Add isPcr to your path. **Note: make sure to write blat, not blast**
```bash
module load UCSC-executables/12may2022
echo $PATH
export PATH=/cm/shared/apps/UCSC-executables/12may2022/blat:$PATH
```

You have now created a file, say candidate_primers.tsv that has three columns, separated by TAB: the first is the primer pair id, the . The sequences are invented, and the primer Ids are arbitrary:
```csv
PrimerPair1  ATGCGGAGAGGATG GTGTGCTGGCGG
ParimerPair2  GTAGCCCGTGCAC  CTGTGAGTTCTCT
```
Run your primers against the transcriptome. A fasta file with predicted PCR products will appear for all your primers
```bash
isPcr StrandedTrinityAss_95.fasta candidate_primers.tsv pcr_predicted.fasta -out=fa
```

Check the pcr_predicted.fasta file. It should have the PCR products that you expect. If it has something else, discard that pair of primers.


## Check the PROBE sequence for specificity against the transcriptome
Now we will use the expected PCR sequence to check that our probes will hit only the specific mRNA sequence(s) that we are interested in.
We will keep the absolute path to the output file in the $pcrfa variable like so:
```
export pcrfa=$(pwd)/pcr_predicted.fasta
echo $pcrfa
```

Now we will move to the folder 27-Transcriptome_Blast/
Once you have navigated there, load the module blast+
```bash
module load blast+/2.7.1
blastn -query $pcrfa -db Dlaevis_BlastDB/Dlaevis_RNA -out resultados.out -evalue 1e06
less resultados.out
```
Ignore the warnings from FASTA-Reader, it is just warning you that the primer sequences are mentioned in the fasta identifiers.
Check your results, you should only see alignments to your intended targets or closely related sequences (isoforms) in the Plus/Plus orientation.
If not, then your probe will probably have unspecific binding, therefore, discard the specific pair of primers that made that hit

Pick at least two pairs of primers and order them. This will save you headaches if a particular pair does not amplify.
Good luck with your experiments.

To come: [how to check for sequence and Antisense orientation](Tutoriales/checarOri.md) when you get back the Sanger sequences from your primers

