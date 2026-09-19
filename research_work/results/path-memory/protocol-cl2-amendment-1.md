# CL-2 stage 1, amendment 1: the held-out lens is withdrawn, before any run

Declared 19 September 2026, after the first attempt to run `cl2.py` stopped while building its inputs and before any result was produced. Nothing in `protocol-cl2.md` is edited.

**What was wrong.** The protocol said J1538+5817 had "light profile, masses and KCWI profile all prepared" and would be scored under every fitted law as a held-out lens. Its KCWI profile and population masses are prepared; its light profile is not. The [light-profile audit](../slacs-light-profile-audit/report.md) records it as `missing_published_components` and says a fallback radius must not be substituted and a second component must not be invented. I declared a held-out system without checking the audit that governs it.

**What changes.** J1538+5817 is withdrawn from stage 1. No replacement is declared: the other KCWI lens outside CL-1's six, J0330-0020, has no light profile either, and no validation or test lens has prepared inputs. The reserved information of this stage is therefore only what section 0(6) also listed — the SPARC validation and test splits, fitted on the training split alone, and the cross-observable predictions (lensing-only fits predicting stellar motions, motion-only fits predicting Einstein radii). Every reading that the protocol wrote for the held-out lens is reported as "not available", never as a pass.

**What does not change.** The family, the widths, the blocks, the decision rule, the gates and their tolerances, the sensitivities and the transfer tests are as declared. The driver is changed only to tolerate the absence of the held-out lens and to record it.
