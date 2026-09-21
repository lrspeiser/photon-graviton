# Area/stress construction checkpoint

Research parent: `0a44c76afb4210e5b28f581a18c862dd4f76bd0a`.
Protocol: `0efc843585608d2503f3fa9ef674d30f019723dd`.
Frozen source/results: `be1a6608726d19b933cd6f1039f1023e3a137456`.

Independent construction quick/full workflow `35563757918` passed both modes.
The independent regeneration/freeze workflow `35563775851` passed and committed
all six result files plus provenance. This checkpoint itself triggers ordinary
construction, validation and general regression checks; their final status must
be read from GitHub rather than inferred from this document.

## Actual result

**Constructed:** a coherent area-reference move `B(W-I)/(2i)` using the inherited
completed-move coefficient x. It has an exact linear-curvature term, keeps its
quadratic completion, and passes finite Hermiticity/covariance/locked-evolution
checks. The equal-weight reference path and quadrature phase are explicit new
architecture, not a uniqueness or symmetry-protection proof.

**Derived:** the classical quadratic frame kernel from that operator, retaining
both quadratures of all frame and connection components. An explicit local
placement of the inherited kinetic trace preserves the fresh linear constraint
algebra. There are two positive tensor modes per signed nonzero momentum without
inserting a TT projector. Both polarizations and three propagation directions
extrapolate toward `c_g^2/(U*x)=1/8`. This is a conditional classical tensor
result, not a full quantum-gravity pole or a photon/gravity speed equality.

**Included:** volume normalization, full spatial Wilson metric, spin links,
electromagnetic links, wall/slab terms, and complete candidate source derivatives.
Applying the same completion AFTER second quantization generates quartic
interactions and quartic stress. They are present in the model definition and
verified in a complete eight-mode/256-state Fock control. They are not silently
included in the larger-volume Gaussian one-loop response.

**Measured:** the same six-species Gaussian matter response to photon and frame
sources on 3D L=12,16,24, with several momenta, all slab bands, mirrors/doublers,
fixed filling controls and no separate coefficient adjustment. The completed
candidate has nonzero raw uniform frame response; no contact was discarded.
These raw values are not physical graviton masses.

**Failed:** the new Gaussian coordinate-stress Ward gate. The analytic small-q
limit is `x^2 ||{O,beta}+3I||`, approximately 0.1881249332 for the declared test.
The measured defect divided by q approaches 0.1880930648 at L=256. An independent
finite graph derivative verifies the source formula to about 2e-11. This is a
leading source defect, not merely a higher-order cutoff error. Internal spin and
charge covariance do not fix it. The failure belongs to this declared Gaussian
candidate; full interacting and alternative transport/projector completions have
not been ruled out.

## Physical gate status

- finite area-reference construction: pass at declared scope;
- placement-matched classical quadratic frame constraints: pass;
- internal U(1)/Spin(3) source consistency and small Fock controls: pass;
- leading Gaussian coordinate-stress consistency: fail;
- full nonlinear gravitational scalar/lapse constraint construction: not established;
- full interacting 3D photon and frame pole matching: not performed / prerequisites unmet;
- sector-specific fitted compensation: none;
- historical full multipair result: byte-identical, checked by SHA-256.

The next admissible revision must derive how the regulator's local transport
projectors and the induced interactions contribute to stress. It must fix the
actual coordinate identity while retaining the successful area and linear mode
results, rather than selecting a photon or gravity coefficient to restore an
apparent numerical match. See README.md and the exact executable outputs.
