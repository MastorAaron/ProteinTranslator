# ProteinTranslator
Rosalind.info coding project

**Environment Information:**
  - Python version: 3.11.5
  - Biopython version: 1.85

**Development Notes:**
This program is made to fulfill the problem set from  "https://rosalind.info/problems/ptra/" 

This is most likely meant to be an application of a naive algorithm and brute forcing the usage of several Codon Tables implemented in BioPython's toolset.

As of now, 6 Codon tables are deprecated: Tables 7,8, and 17-20 consistently raised exceptions when brute forcing single codon recipes.  A try-except block was implemented to ignore these. 

These are likely either broken or much too niche for this project. However, after implementing the Brute Force "probeTables()" function to test each Codon generated recursively in the "generateCombos(n)" function, it is most likely that these tables are nonfunctional. This could also be due to my version of Python 3.11.5 and BioPython 1.85, respectively, however, this is also unlikely since I have recently updated both.

Initially, I got the sample set working by setting the BioPython, SeqIO translate function's stop_symbol to '', this blank coupled with using a HammingDistance score between the given translation and the generated translation and yielded the most ideal score 0 for two of the 6 tables' scores.

However, this was before I became aware that there were more than 6 Codon tables available, considering I used range(1,10), and 7 and 8 were broken. Only after seeing the maxTableID did I determine that there are 27 tables actually available, (of the total 33 with the 6 broken).

An initially strict HammingDistance scorer was implemented and then a more loose version was implemented after determined that in some instances StopCodons ought to be Q (glutamates) instead of stops. I then implemented the ProbeTables() to use all 33 tables for all combos that generate Q as an amino acid.

  To rank scores of each codon table, each translation was compared to the given expected amino acid string. The project initially used a strict HammingDistance score, which was made looser to account for ambiguous **STOP codon** that are sometimes translated as Q (glutamine) or W (Tryptophan) depending on the domain of life.

**Biological Considerations:**
The next day, I determined I ought to probe all possible codons with each table available in order to determine if there were any other "Ambiguous" Codon recipes. This is due to the nature of well _Nature_. Different domains of life have different available and useful recipes for codons, and sometimes they've shifted what each recipe yields. 

In this case, StopCodons are a known wildcard in other domains of life. There is also an extended Codon table that makes some StopCodon recipes into O and U, etc, however, these two were among the list that are not used to generate Amino acid strings for the problem (i.e., amino acids except for B, J, O, U, X, and Z).

After using the Probe to determine if there were any other Ambiguous Codon Recipes, 13 Codons stood out as being atypical across all 33 Codon Table Domains. 

These 13 codons have been marked and are going to be used to implement a function that flags their presence and hopefully finds the appropriate Codon Table faster.

Before this new strategy, the only thing I could think of was changing the stop_symbol and the strictness of the hammingDistance to match StopCodons to known alternates (e.g., Q, W).

**Next Steps:**
Integrate ambiguous codon flags into the core strategy for early table elimination or prioritization.

Use the expected protein string to validate translations more efficiently before resorting to full brute-force.
