### <font color="#4bacc6">Sort-seq  vs Growth based assay scoring method</font>

Problem:
- Many tools exist to score variants in high-throughput growth based assays. None handle the Sort-seq functional scores produced by a <u>weighted average corresponding to median fluorescence and normalized sums of counts</u> of a variant in each bin.

Solution:
- Produce a <u>ratio score</u> that is compatible with the assumptions made in these programs. 
	- Test these tools with a ratio of the fold change in counts between the highest and lowest bins for each variant.
	- [ ] Is it a good assumption that cells with benign variants are most likely to sorted into the top bin? This will make it difficult if not impossible to <u>distinguish gain of function </u>variants.
- Compare performance of DESeq2, Enrich2, DimSum STEAM, and Rosace against results of using our functional scoring method incorporating all binned data
	- [ ] metric?

### <font color="#4bacc6">Normalization</font>

Problem:
- Between technical replicates
- Averaging scores and errors across biological replicates
- Batch effects

Solution:
- Weighting based on read count in technical replicates. This assumes higher counts produce more accurate functional score estimates.
- Double centering on mean of distribution of scores for synonymous and nonsense scores using local and global scores.
	- [ ] median?

### <font color="#4bacc6">Interpretable Quantitative Scores</font>

Problem:
- After normalization, it is difficult to orient yourself to determine how a variant performed in the assay.
- We want to assess if there are possible gain of function variants. 
Solution:
- Set magnitude of difference between observed median scores for nonsense and synonymous variants to be the unit of measure.
	- syn avg = 1
	- nonsense avg = 0
### <font color="#4bacc6">Characterizing Multimodal Score Distributions</font>

Problem: 
- Functional scores in Sort-seq assays exhibit multiple distinct modes/peaks by design. <u>Summary statistics</u> assuming a normal distribution will not accurately represent data of this shape. 
- Clinically, we want to estimate the distribution of subgroups and assign thresholds to <u>classify variants </u>(pathogenic/benign/gain of function).

Solution:
- We will fit mixed distributions to the functional score data after normalizing. 
- As a quality check, run statistical tests to confirm the distribution is multimodal.
	- Hartigan's dip test (R diptest package)
- We could also consider non-parametric methods.
- Classification
	- Barcoded data - KS test


### <font color="#4bacc6">Testing package on other Sort-seq data</font>

- GLI2
- OTX2/CRX comparison
- See list on trello for MAVE Sort-seq datasets
- See paper Dustin sent for applicability in other adjacent fields

---

As time allows... 


### <font color="#4bacc6">Error Analysis</font>

Problem: 
- Error is not estimated (and so inherently underestimated) by s.d. based method currently.

Solution:
- We want to estimate errors that can arise from various points in the experiment.
- Use nonparametric methods when possible.
- Compare error between replicates, returning average.

### <font color="#4bacc6">Producing Counts: dropout and quality </font>

Problem:
- We lose a large proportion of reads through the processing of raw fastq into variant sequence counts for each sample.
- For quality filtering, assess if there is a need for minimum read count cutoff and determine selection methodology. Likely not warranted for strictly aligned paired-end reads, but it's difficult to intuit the expected paired error rate at the scale of tens of millions of reads per sample.

Solution:
- Test DimSum WRAP
	- Allow for mismatches in flanking regions
- Trim flanking regions instead of adapters so there is a shorter sequence to merge
- Allow for some mismatches to be resolved when merging paired-end reads instead of automatically tossing?
- Could also deliver score even if variant was not represented in all samples (such as in DimSum)
### <font color="#4bacc6">Plotting</font>

Need to see if other tools fill this gap in a way that I can use.
- tacking on python existing code
### <font color="#4bacc6">Workflows</font>

Nextflow module.