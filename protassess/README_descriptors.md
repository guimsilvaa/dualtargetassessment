# Descriptor Definitions

## Surface-to-Volume Ratio (SVR)
As importantly cited in Volkamer, A., Kuhn, D., Rippmann, F., & Rarey, M. (2012). DoGSiteScorer: a web server for automatic binding site prediction, analysis and druggability assessment. Bioinformatics, 28(15), 2074-2075. 
Higher values generally indicate:

* more exposed cavities
* shallower grooves
* irregular surfaces

Lower values generally indicate:

* more enclosed pockets
* compact orthosteric sites
```text
SVR = Surface Area / Volume
```
---

## Compactness

Measures cavity compactness/buriedness. As referenced in: Hajduk, P. J., Huth, J. R., & Fesik, S. W. (2005). Druggability indices for protein targets derived from NMR-based screening data. Journal of medicinal chemistry, 48(7), 2518-2525.

Higher values:

* more compact pockets
* enclosed orthosteric cavities

Lower values:

* elongated channels
* solvent-exposed grooves
```text
Compactness = Volume / Surface Area
```
---

## Sphericity

Approximates cavity globularity, as importantly referenced in: Coleman, R. G., & Sharp, K. A. (2009). Shape and evolution of thermostable protein structure. Proteins: Structure, Function, and Bioinformatics, 78(2), 420–433.

Values closer to 1:

* more spherical cavities

Lower values:

* elongated tunnels
* irregular binding grooves
```text
Sphericity = π^(1/3) (6V)^(2/3) / A
```
---

# Residue Composition Metrics

## Hydrophobic Fraction

Fraction of cavity-lining residues classified as hydrophobic:

* ALA
* VAL
* LEU
* ILE
* MET
* PRO
* PHE
* TRP

---

## Aromatic Fraction

Fraction of aromatic cavity residues:

* PHE
* TYR
* TRP
* HIS

This descriptor is particularly informative for:

* aromatic stacking environments
* ACHE-like aromatic gorges
* CNS-target orthosteric sites

---

## Charged Fraction

Fraction of charged residues:

* ASP
* GLU
* LYS
* ARG
* HIS

Useful for:

* electrostatic characterization
* ATP-binding pockets
* polar recognition environments
