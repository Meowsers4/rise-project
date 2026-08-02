# Citations — SOD1 FEP Variant Pipeline

Curated, freely-viewable literature relevant to this project. Every entry was
verified against PubMed/PMC on 2026-08-01.

**Access legend**
- **Free full text** — verified free on PubMed Central (PMC) or open access at the publisher.
- **Abstract-only** — verified record exists, but I could NOT confirm a free full text. Check
  Google Scholar / Europe PMC / your SCC library proxy before relying on it.

> ⚠️ **Unverified source:** `data/variants.csv` cites `Kumar2017_TableS1` as a control source
> for many apo-monomer ΔΔG values. I could not locate a "Kumar 2017" SOD1-stability paper on
> PubMed. **Confirm this citation before publishing anything that rests on it.** The Lindberg
> 2005 and Nordlund/Oliveberg 2006 entries below are confirmed and carry the same numbers.

---

## §1 SOD1 target and stability (the biology)

1. **Wells NGM, Tillinghast GA, O'Neil AL, Smith CA (2021).** *Free energy calculations of
   ALS-causing SOD1 mutants reveal common perturbations to stability and dynamics along the
   maturation pathway.* Protein Sci 30(9):1804–1817.
   - PMID [34076319](https://pubmed.ncbi.nlm.nih.gov/34076319/); PMCID [PMC8376412](https://pmc.ncbi.nlm.nih.gov/articles/PMC8376412/); DOI 10.1002/pro.4132. **Free full text.**
   - **Relevance:** the closest published analog to this project. Uses **pmx** (the same engine
     as our Stage 3), computes ΔΔG for **disulfide-reduced apo** SOD1 (our §2.1 species), and
     averages over OPLSAA/AMBER force fields (see insight below). Their ΔΔG_dimerization values
     agree with ITC/SEC experiment. Read this first.

2. **Lindberg MJ, Byström R, Boknäs N, Andersen PM, Oliveberg M (2005).** *Systematically
   perturbed folding patterns of amyotrophic lateral sclerosis (ALS)-associated SOD1 mutants.*
   Proc Natl Acad Sci USA 102(28):9754–9759.
   - PMID [15987780](https://pubmed.ncbi.nlm.nih.gov/15987780/); PMCID [PMC1174986](https://pmc.ncbi.nlm.nih.gov/articles/PMC1174986/); DOI 10.1073/pnas.0501957102. **Free full text.**
   - **Relevance:** experimental gate source (`Lindberg2005_PNAS_102_9754` in `variants.csv`).
     Defines class-1 (monomer destabilization) / class-2 (dimer interface) perturbations, and
     shows **net-charge-changing mutations behave differently** (ΔΔG vs survival correlation,
     R = 0.91, only holds for charge-conserving). This is the experimental basis for our
     charge-neutral gate subset (config/pipeline.yaml:161-166).

3. **Nordlund A, Oliveberg M (2006).** *Folding of Cu/Zn superoxide dismutase suggests
   structural hotspots for gain of neurotoxic function in ALS: parallels to precursors in
   amyloid disease.* Proc Natl Acad Sci USA 103(27):10218–10223.
   - PMID [16798882](https://pubmed.ncbi.nlm.nih.gov/16798882/); DOI 10.1073/pnas.0601696103. **Free PMC.**
   - **Relevance:** experimental gate source (`NordlundOliveberg2006_PNAS_103_10218`).
   Apo-monomer folding analysis — supports apo-first + monomer default (§2.1, README §4 Stage 1).

4. **Nordlund A, Oliveberg M (2008).** *SOD1-associated ALS: a promising system for
   elucidating the origin of protein-misfolding disease.* HFSP J 2(6):354–364.
   - PMID [19436494](https://pubmed.ncbi.nlm.nih.gov/19436494/); DOI 10.2976/1.2995726. **Free PMC.**
   - **Relevance:** accessible review of SOD1 misfolding/aggregation — good background for the
     introduction and for Stage 4 mechanism interpretation.

5. **Nordlund A, Leinartaite L, Saraboji K, Aisenbrey C, Gröbner G, Zetterström P, Danielsson J,
   Logan DT, Oliveberg M (2009).** *Functional features cause misfolding of the ALS-provoking
   enzyme SOD1.* Proc Natl Acad Sci USA 106(24):9667–9672.
   - PMID [19497878](https://pubmed.ncbi.nlm.nih.gov/19497878/); DOI 10.1073/pnas.0812046106. **Free PMC.**
   - **Relevance:** apo/cofactor-loss promotes misfolding — context for why the apo, reduced form
     is the disease-relevant species to simulate.

6. **Svensson AKE, Bilsel O, Kayatekin C, Adefusika JA, Zitzewitz JA, Matthews CR (2010).**
   *Metal-free ALS variants of dimeric human Cu,Zn-superoxide dismutase have enhanced
   populations of monomeric species.* PLoS ONE 5(4):e10064.
   - PMID [20404910](https://pubmed.ncbi.nlm.nih.gov/20404910/); DOI 10.1371/journal.pone.0010064. **Free (open access).**
   - **Relevance:** experimental support for apo monomer being the destabilized, aggregation-prone
     species — justifies our apo/monomer design choice.

7. **Hsueh SCC, Nijland M, Peng X, Hilton B, Plotkin SS (2022).** *First Principles Calculation
   of Protein-Protein Dimer Affinities of ALS-Associated SOD1 Mutants.* Front Mol Biosci 9:845013.
   - PMID [35402516](https://pubmed.ncbi.nlm.nih.gov/35402516/); DOI 10.3389/fmolb.2022.845013. **Free (open access).**
   - **Relevance:** first-principles dimer-binding ΔΔG for SOD1 variants — relevant if/when we
     extend to dimer legs (interface-adjacent variants, per README §4 Stage 1).

---

## §2 Method — GROMACS + pmx alchemical FEP (our Stage 3 engine)

8. **Gapsys V, Michielssens S, Seeliger D, de Groot BL (2015).** *pmx: Automated protein
   structure and topology generation for alchemical perturbations.* J Comput Chem 36(5):348–354.
   - PMID [25487359](https://pubmed.ncbi.nlm.nih.gov/25487359/); DOI 10.1002/jcc.23804. **Free PMC.**
   - **Relevance:** the pmx method paper — the engine our Stage 3 (`src/fep/pmx_engine.py`) is
     built on (hybrid topology, mutate/gentop workflow, mutff45 force fields).

9. **Seeliger D, de Groot BL (2010).** *Protein thermostability calculations using alchemical
   free energy simulations.* Biophys J 98(10):2309–2316.
   - PMID [20483340](https://pubmed.ncbi.nlm.nih.gov/20483340/); DOI 10.1016/j.bpj.2010.01.051. **Free PMC.**
   - **Relevance:** the alchemical ΔΔG-of-folding protocol our thermodynamic cycle (folded leg vs
     unfolded tripeptide leg) follows.

10. **Gapsys V, Michielssens S, Seeliger D, de Groot BL (2016).** *Accurate and Rigorous
    Prediction of the Changes in Protein Free Energies in a Large-Scale Mutation Scan.*
    Angew Chem Int Ed 55(26):7364–7368.
    - PMID [27122231](https://pubmed.ncbi.nlm.nih.gov/27122231/); PMCID [PMC5074281](https://pmc.ncbi.nlm.nih.gov/articles/PMC5074281/); DOI 10.1002/anie.201510054. **Free PMC.**
    - **Relevance:** 762-mutation validation of pmx FEP. AUE ≈ 0.8–1.0 kcal/mol with consensus
      force fields; error splits ~equally between sampling, force-field, and experimental
      uncertainty. **Includes an explicit charge-conserving vs charge-changing error breakdown** —
      direct support for our charge-neutral gate decision. Their "average over ≥2 force fields"
      conclusion is the basis for idea-02 and a recommended M3 enhancement.

11. **Gapsys V, Seeliger D, de Groot BL (2012).** *New Soft-Core Potential Function for
    Molecular Dynamics Based Alchemical Free Energy Calculations.* J Chem Theory Comput
    8(7):2373–2382.
    - PMID [26588970](https://pubmed.ncbi.nlm.nih.gov/26588970/); DOI 10.1021/ct300220p. **Abstract-only** (ACS).
    - **Relevance:** the gapsys soft-core alternative to the beutler form currently configured
      (config/pipeline.yaml:116-127). Relevant for the post-benchmark soft-core switch (idea-15).

12. **Gapsys V, de Groot BL (2017).** *pmx Webserver: A User Friendly Interface for
    Alchemistry.* J Chem Inf Model 57(2):109–114.
    - PMID [28181802](https://pubmed.ncbi.nlm.nih.gov/28181802/); DOI 10.1021/acs.jcim.6b00498. **Abstract-only** (ACS).
    - **Relevance:** alchemical setup database and protocol examples — cross-check for λ-window
      count, soft-core settings, and equilibration choices.

---

## §3 Prescreen (Stage 2) and validation context

13. **Tokuriki N, Stricher F, Schymkowitz J, Serrano L, Tawfik DS (2007).** *The stability
    effects of protein mutations appear to be universally distributed.* J Mol Biol
    369(5):1318–1332.
    - PMID [17482644](https://pubmed.ncbi.nlm.nih.gov/17482644/); DOI 10.1016/j.jmb.2007.03.069. **Abstract-only** (Elsevier).
    - **Relevance:** FoldX-based ΔΔG study over many proteins — useful for interpreting expected
      accuracy/ranking behavior of our Stage 2 prescreen.

14. **Schymkowitz J, Borg J, Stricher F, Nys R, Rousseau F, Serrano L (2005).** *The FoldX web
    server: an online force field.* Nucleic Acids Res 33:W382–W388.
    - PMID [15980494](https://pubmed.ncbi.nlm.nih.gov/15980494/); DOI 10.1093/nar/gki387. **Abstract-only** (NAR may be open — check).
    - **Relevance:** primary reference for FoldX BuildModel (Stage 2 prescreen). Record the
      binary path in config when installed; FoldX is a licensed/registered tool, not conda.

15. **Tiberti M, Terkelsen T, Degn K, Beltrame L, Cremers TC, da Piedade I, Di Marco M,
    Maiani E, Papaleo E (2022).** *MutateX: an automated pipeline for in silico saturation
    mutagenesis of protein structures and structural ensembles.* Brief Bioinform 23(3):bbac074.
    - PMID [35323860](https://pubmed.ncbi.nlm.nih.gov/35323860/); DOI 10.1093/bib/bbac074. **Free article.**
    - **Relevance:** automates FoldX saturation mutagenesis across an ensemble — directly usable
      for the prescreen and for project-idea-01 (full stability landscape scan).

---

## §4 Experimental controls — mapping to this repo

The gate (`validation.gate_on: positive_control`, `ddg_gate_column: exp_ddg`) compares our FEP
ΔΔG against the **apo-monomer** ΔΔG column in `data/variants.csv`. Two of the three confirmed
sources are above:

| `exp_source` tag in variants.csv | Citation | Status |
|---|---|---|
| `Kumar2017_TableS1` | **UNVERIFIED — not found on PubMed** | ⚠️ confirm |
| `Lindberg2005_PNAS_102_9754` | Lindberg et al. 2005, PNAS (entry 2) | confirmed |
| `NordlundOliveberg2006_PNAS_103_10218` | Nordlund & Oliveberg 2006, PNAS (entry 3) | confirmed |

---

## Key take-aways for the project right now

1. **Consensus force-field averaging** (Gapsys 2016; Wells 2021) is the single biggest accuracy
   lever and requires no engine change — the pmx `mutff45` dir ships `amber99sbmut`,
   `amber99sb-star-ildn-mut`, `charmm22starmut`, `oplsaamut` (config/pipeline.yaml:107-108).
   Run the gate subset under 2–3 FFs and average.
2. **`max_rmse_kcal: 1.5`** (config:174) is defensible against the ~0.8–1.0 kcal/mol AUE Gapsys
   reports; possibly tightenable to ~1.2 after the first benchmark.
3. **Charge-neutral gate is literature-supported** (Gapsys 2016 error breakdown; Lindberg 2005
   charge effect) — cite these in `validation_gate.json` / config comments.
