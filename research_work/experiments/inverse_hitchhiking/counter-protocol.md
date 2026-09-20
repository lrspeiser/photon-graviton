# IH-1 signed-steering follow-up

Declared after the positive-only results, before this fit. Preserve evidence-v1.
The primary inverse table has 190 rotation rows with required additional
acceleration below zero under the frozen ordinary-matter model. Test an
outward steering component as the owner's permitted sign variation.

Keep primary T's five-parameter positive response a_plus=a_sat*expit(F theta).
Add a bounded opposite coherent response:
S=expit(k*(ln(g_b/1e-10 m s^-2)-x_flip)),
a_H=a_plus-epsilon*g_b*S,
g_total=g_b*(1-epsilon*S)+a_plus.
epsilon in [0,0.95], x_flip in [-8,8], k fixed {0.5,1,2}.
The total ordinary-source acceleration remains positive; the additional
response can change sign. There is no negative gravitating mass inserted.
Interpretation is an effective opposite turning contribution, not established
binding microphysics or a repair of the ordinary-matter observations.

For each of the four original capacities and three k values, run galaxy-only
and joint training fits: 24 fits, three starts each. Start from the matching
primary T solution, with (epsilon,x_flip)=(0.1,0),(0.5,-2),(0.8,2).
Other parameter bounds, max500 evaluations, residuals, data and split stay as
declared. Select separately by the original validation criteria before scoring
these candidates on the test targets. The test set has already been exposed by
the primary campaign; this follow-up is not a fresh blind validation.

The epsilon=0 limit must reproduce the primary T predictions. Verify the
positive-total-acceleration bound. Source manifests pin any compatible utility
changes, and preserve the old source hashes against their old Git commits.
Project the signed winner to the same Coma primary and reach/eta sensitivities;
do not select on Coma. Check pressure600/1200 and projection128/256 as before.
Publish improvement or degradation, both signed and positive-only formulas,
and retain every failed or boundary-limited optimizer.
