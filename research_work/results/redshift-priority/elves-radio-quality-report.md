# ELVES radio quality and source links

9 September 2026. The predeclared metadata-only join finds 14 positional matches among all 26 pending targets. Nine matches additionally have an explicit AGC/UGC alias on the previously supported SIMBAD counterpart. Five remain positional associations only. Twelve targets have no entry inside the fixed 30-arcsecond cone; absence does not establish that no radio measurement exists.

## Results retained for every matching target

| ELVES target | AGC | Separation, arcsec | Source notes | Explicit supported alias |
|---|---:|---:|---|---|
| dw0243p1643 | 121135 | 1.78 | 1 | No |
| DDO047 | 3974 | 4.66 | 1 | Yes |
| UGC04115 | 4115 | 6.85 | 1 | Yes |
| LV J0913+1937 | 739005 | 0.68 | 1 | Yes |
| LSBC D634-03 | 191791 | 1.17 | 2 | No |
| KDG056 | 191706 | 3.78 | 1 | Yes |
| SexB | 5373 | 2.36 | 1 | Yes |
| LV J1017+2922 | 200232 | 0.90 | 1 | No |
| AGC208399 | 208399 | 4.60 | 1 | Yes |
| LV J1030+0607 | 203709 | 2.97 | 1 | Yes |
| LV J1000+3032 | 205590 | 1.07 | 1 | Yes |
| dw1138p0036 | 213171 | 0.71 | 1 | No |
| AGC740112 | 740112 | 1.91 | 1 | Yes |
| KKH86 | 231980 | 0.72 | 1 | No |

The saved JSON retains all 26 targets, including the twelve unmatched cases. All nine explicitly supported matches have code 1. No matched entry carries c, p or r. This establishes absence of those catalog flags, not absence of every possible contaminating source or measurement bias. None of the targets is excluded or certified fresh by this pass.

## Why these flags matter

The [Yu source documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJS/261/21?format=html&tex=true) marks overlapping hydrogen emission with c, a projected neighbor without declared confusion with p, and low spectral weight at an emission edge with r. Code 1/2 refers back to the ALFALFA detection classification. In the [ALFALFA documentation](https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/ApJ/861/49?format=html&tex=true), class 2 is a lower-significance radio detection associated with an already known optical redshift.

Thus the possible LSBC D634-03 match has a selection dependency worth retaining: its inclusion is partly supported by prior spectroscopy. That does not invalidate its measurement or prove the ELVES row uses this radio value. We must resolve the identity and adopted source first. The other thirteen matches carry class 1, which also does not imply a random or complete sample of galaxies.

The [Yu paper](https://arxiv.org/abs/2203.13404), page 2, describes the original survey band using heliocentric velocities and says it uses the final ALFALFA spectra. A targeted search found no explicit optical/radio convention statement for its reprocessed central-velocity column. The earlier Haynes optical-convention definition remains useful source evidence, but this pass does not certify every downstream conversion. Local PDF SHA256: f88b967c812ab6fb21629e8e7be2dd626dee362cc24d62c4ba33a4b7065918cd. Only paragraphs matching frame/convention terms were displayed; no target table values were inspected.

## Reproduction and scope

The compressed CDS source is `https://cdsarc.cds.unistra.fr/ftp/J/ApJS/261/21/table2.dat.gz`. Decompress it to a scratch file, then run:

```sh
python research_work/results/redshift-priority/elves_radio_quality.py /path/to/yu2022-table2.dat
```

The script pins the decompressed table hash to 816431e3254754d86fc2fccba6fd76b0347d0357f959b57bbfa34ae9d4283225, verifies 29,958 unique AGC entries, and reads only identifiers, optical counterpart positions and notes. For every source row, replacing all bytes outside those fields with X leaves the parsed result unchanged. This verifies that the extraction does not depend on numeric distance, velocity, uncertainty, flux or mass fields. Full table bytes are downloaded into the scratch cache, but numeric outcomes are neither parsed nor printed. This is an audit of access and usage, not a claim that the outcome bytes are unavailable to the machine.

The position calculation uses the established spherical great-circle formula, not a proposed physical law. It uses the source's optical counterpart position, not the broad radio-beam center. Matching has no adjustable radius per galaxy and no residual-based tie-breaking. Input hashes include the frozen protocol, publisher features and current identity evidence. The source table and filtered output support a reproducible join; they do not certify the original ELVES measurement or historical nonexposure.

## Consequence for the next step

Nine targets now have explicit links to a radio source with a documented detection class. Continue the provenance audit using those links, and retain the five positional cases as unresolved. In particular, the LV J1017+2922 match to AGC200232 still lacks an explicit PGC4231240 bridge, and AGC740112 still lacks an explicit PGC5808772 bridge. Quality flags do not resolve those historical group overlaps. A source with no confusion flag is not a measured void path, and no environmental rate term is justified by these data.

Current count remains 26 pending and zero certified fresh. No prediction, residual, interval or fitted rate changed. The next numerical comparison must still wait for a frozen common model and uncertainty account; these flags alone do not supply a peculiar-motion distribution.
