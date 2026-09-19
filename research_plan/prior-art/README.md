# Prior-art audit: gravitational memory and collective orbital dynamics

**Audit date:** 19 September 2026. **Scientific baseline:** `3a80fec8ae7d45cc6906c20cb9c537c3ac2c9698`.

**Conclusion:** several central ingredients are established, and the PM interpolation function is an exact rediscovery. This review has not established that the complete two-stage RUT model or its verified numerical results were previously published. It also has not established that they are original. Novelty must attach to a precisely stated contribution, not to the words *gravity*, *memory*, or *collective stabilization*.

This package consolidates the prior conversation's findings and expands the comparison into active lattices, chemotaxis, inertial aggregation, transported memory, system-bath models, stellar dynamics, and numerical response methods. It contains **34 inspected source records plus three separately labeled, unretrieved leads**. Inspection depth varies: exact equation matches use primary full text; abstract-level comparisons are labeled as such. It is not a systematic review with guaranteed coverage or a patent opinion.

## Read in this order

1. [Detailed equation and mechanism audit](audit.md): the exact matches, the partial mappings, and what does **not** follow from them.
2. [Claim register](claims.json): permitted wording, prohibited overclaims, supporting sources and remaining work.
3. [Sources and inspection levels](sources.json), with [BibTeX](references.bib).
4. [Comparison tasks for the local research agent](next-comparisons.md).
5. [Search record and limitations](search-log.md).
6. [Offline checks](verify.py) and [recorded verification](verification.json).

## Most consequential findings

| Finding | Consequence |
|---|---|
| Famaey and Binney (2005), equation 5, prints the Bekenstein toy interpolation in a form identical to `mu_PM`. | Attribute the function and associated spherical acceleration relation. Do not call either a new law. Preserve the original source's regime qualifications. |
| The AQUAL field equation and its action/constitutive derivative relation are established in Bekenstein and Milgrom (1984). | A new implementation or choice of response function is not the invention of this field-equation architecture. |
| Thomson, Durey and Rosales (2020) already distinguish geometric from memory-dependent collective instability in a particle-written field. | That distinction alone cannot be the novelty claim of RUT stage 5. Its particular kernel, dynamics and quantitative predictions still require comparison. |
| Inertial chemotaxis models already combine moving sources, produced/decaying fields, gradient response and kinetic populations. | Search and benchmark beyond gravitational terminology. |
| Chavanis and Sire (2008), and ordinary stellar-dynamics work, already study suppression of attractive collective collapse by pressure or velocity dispersion. | The particular rotating memory-mediated threshold may be a result; the general warming-stabilizes principle is not new. |
| Li and Valani's 2026 preprint reports persistent rotating/translating clusters from individually unstable droplets. | Do not claim first collective stabilization through a shared wavefield. Their bounce/coalescence instability and driven apparatus differ from our orbital instability. |
| The RUT cascade is exactly a second-order linear response with two poles; its stationary mass-sourced field reduces to attractive Gaussian pair interactions. | Locate possible novelty in the coupled dynamics and predictive consequences, not in auxiliary variables or stationary attraction alone. |

## What stays open

The specific RUT model combines an attractive footprint, inertial orbital populations, two response times and reciprocal source/force coupling. No complete equivalence to another physical theory was established in this search. Candidate contributions include a general support-torque relation, a controlled stability/phase diagram and a predictive history signature. Each remains **originality unverified** until an equation-level benchmark against the closest models is completed.

The fictional-universe premise remains unchanged. Identifying prior mathematics neither imports another theory's cosmology nor proves this model describes measured gravity. Novelty, numerical correctness and observational success are three independent questions.

## Preservation and reproduction

This is a literature/provenance update, not a rerun or reinterpretation of scientific gates. No file under `research_work/results/`, numerical archive, failed gate or frozen experiment protocol is edited. Stage 8's verification limitations remain exactly what its own records state. The audit is linked from `../formula-provenance.md` without rewriting historical declarations.

From the repository root, with Python, SymPy, NumPy and SciPy available:

```sh
python -B research_plan/prior-art/verify.py
```

The checks verify algebraic identities, citation-register integrity and this package's manifest. They do not certify originality, a complete physical equivalence, or the numerical correctness of the gravity simulations. The historical full suite was not rerun for this documentation update.

**Standing publication wording:** “We investigate a specified gravitational-memory hypothesis assembled from established mathematical structures. Several components have exact antecedents. The originality of the complete construction and its particular results remains to be demonstrated against the identified prior models.”
