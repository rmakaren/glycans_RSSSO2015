Samples are technical replicates from one biological sample.

To test for possible source of variation 5 different types of samples
were prepared:

- A - samples pooled after deglycoslylation (after 1st day/step); 24 samples
- B - samples are pooled after dyeing (after 2nd step); 24 samples
- C - samples are pooled after clean-up (3rd step) - 24 samples
- D - samples are not pooled at all - 24 samples
- E - samples are pooled after clean-up (3rd step ) - 96 samples

*Notice*: Sample marked with C and E are pooled at the same time (before running
in the UPLC machine) but will be on different Plates in preparation.

## Problems

In preparation of samples:

- A_24 - after samples were pooled this one contained smaller amount of "liquid"
since it was the last one taken from the pool

In UPLC run:

- pl1_pt1: machine wen wild (big oscillations in pressure) - A_01, A_02, A_03
  were repeated later because of that
- pl1_pt1_R1_pt2:  A_01, A_02, A_03 were repeated together with 2. part normally
- pl1_pt3: everything went ok
- pl1_R2_pl2_pt1: A_04 was repeated and the first part of plate 2 with it
- pl2_pt2: everything went ok
- pl2_pt2_R1_pt3: part 3 from plate 2 run.  Machine stopped working on sample E_87
- pl2_pt2_pt3_R2: after a week when machine started to work normally 
  the remaining samples were run (from E_87) + some strange samples (E_35_R2, E_38_R2, E_39_R2, E_40_R2)
  Strange was defined by lab people looking at chromatograms.

Header information for `./data/areas-20141202.csv` :

- SampleName - represents unique sample names
- glycan - represents glycan structures
- area - area as obtained after integration from chromatograms
- narea - normalized area to percent values

Header information for `./data/chromatograms-20141202.csv` :

- SampleName - represents unique sample names
- x - time
- y - UPLC intensity


Header information for `./data/chromatograms-20141202.csv` :

- SampleName - represents unique sample names
- Sample.Set.Name - internal sample set name (contains part name)
- Sample.Set.Start.Date - when were these samples run on a UPLC machine
- type - type of replicates (A,B,C,D or E)
- id - id :)
- run - if this samples was repeated it should have greater run number
- strange - mark given by lab people after exploring chromatograms

