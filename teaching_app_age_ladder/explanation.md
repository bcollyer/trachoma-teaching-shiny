## What is this model?

This app combines the two extensions built in the companion Shiny
prototypes into the **full model**:

- the **infection-history ladder** (partial immunity via rung-dependent
  recovery and infectiousness) from `teaching_app_ladder`, and
- **age structure** (18 age groups, age-assortative mixing, ageing,
  age-specific mortality, births) from the full model.


### Two nested structures

**Within each age group**, the ladder works exactly as in
`teaching_app_ladder`: people move `S_1 -> I_1 -> S_2 -> I_2 -> ... -> S_n
-> I_n`, with recovery rate $\nu_j$ increasing and infectiousness $l_j$
decreasing at higher rungs, and irreversible scarring (**TS**) then
trichiasis (**TT**) accumulating from repeat infection:

<svg width="100%" height="300" viewBox="-40 0 1080 300" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Within an age group: S1, I1, S2, I2, further rungs, Sn, In, TS, TT boxes with arrows">
  <defs>
    <marker id="a1" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#333333"/>
    </marker>
    <marker id="a2" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#777777"/>
    </marker>
  </defs>
  <rect x="15" y="60" width="100" height="55" rx="8" fill="#eaf2fb" stroke="#2b6cb0" stroke-width="2"/>
  <text x="65" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#2b6cb0">S_1</text>
  <rect x="165" y="60" width="100" height="55" rx="8" fill="#fdecea" stroke="#c0392b" stroke-width="2"/>
  <text x="215" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#c0392b">I_1</text>
  <rect x="315" y="60" width="100" height="55" rx="8" fill="#f1eafb" stroke="#6b46c1" stroke-width="2"/>
  <text x="365" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#6b46c1">S_2</text>
  <rect x="465" y="60" width="100" height="55" rx="8" fill="#fdf0e0" stroke="#b7791f" stroke-width="2"/>
  <text x="515" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#b7791f">I_2</text>
  <text x="632" y="98" text-anchor="middle" font-size="26" font-family="sans-serif" fill="#777777">&#8230;</text>
  <rect x="700" y="60" width="100" height="55" rx="8" fill="#f1eafb" stroke="#6b46c1" stroke-width="2"/>
  <text x="750" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#6b46c1">S_n</text>
  <rect x="850" y="60" width="100" height="55" rx="8" fill="#fdf0e0" stroke="#b7791f" stroke-width="2"/>
  <text x="900" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#b7791f">I_n</text>
  <rect x="850" y="190" width="100" height="55" rx="8" fill="#f0f0f0" stroke="#555555" stroke-width="2"/>
  <text x="900" y="223" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#555555">TS</text>
  <rect x="850" y="280" width="100" height="55" rx="8" fill="#f8e0e0" stroke="#7a1f1f" stroke-width="2"/>
  <text x="900" y="313" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#7a1f1f">TT</text>
  <path d="M115,87 L163,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <path d="M265,87 L313,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <text x="290" y="72" text-anchor="middle" font-size="11" font-family="sans-serif">&#957;1 (up ladder)</text>
  <path d="M415,87 L463,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <path d="M565,100 C 620,125 660,125 698,100" fill="none" stroke="#777777" stroke-width="2" stroke-dasharray="3,4" marker-end="url(#a2)"/>
  <path d="M800,87 L848,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <path d="M850,110 C 820,150 780,150 750,117" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <text x="800" y="160" text-anchor="middle" font-size="11" font-family="sans-serif">&#957;n (capped: SIS cycle)</text>
  <path d="M900,117 L900,188" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <text x="960" y="155" text-anchor="middle" font-size="11" font-family="sans-serif">k_TS</text>
  <path d="M900,247 L900,278" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#a1)"/>
  <text x="960" y="270" text-anchor="middle" font-size="11" font-family="sans-serif">k_TT</text>
</svg>

**Across age groups**, everyone (S, I, TS and TT alike) ages forward
through 18 bands at rate $1/w_a$ (where $w_a$ is the width in years of age
band $a$), with births entering the
youngest naive-susceptible rung and age-specific mortality removing people
from every compartment:

<svg width="100%" height="180" viewBox="0 0 900 180" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Age chain: age group 1 through 18, with ageing arrows and births">
  <defs>
    <marker id="b1" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#333333"/>
    </marker>
  </defs>
  <rect x="10" y="60" width="140" height="55" rx="8" fill="#eef7ee" stroke="#2e7d32" stroke-width="2"/>
  <text x="80" y="83" text-anchor="middle" font-size="13" font-family="sans-serif" fill="#2e7d32">Age 1 (0-1y)</text>
  <text x="80" y="100" text-anchor="middle" font-size="11" font-family="sans-serif" fill="#2e7d32">ladder inside</text>
  <rect x="220" y="60" width="140" height="55" rx="8" fill="#eef7ee" stroke="#2e7d32" stroke-width="2"/>
  <text x="290" y="83" text-anchor="middle" font-size="13" font-family="sans-serif" fill="#2e7d32">Age 2 (1-2y)</text>
  <text x="290" y="100" text-anchor="middle" font-size="11" font-family="sans-serif" fill="#2e7d32">ladder inside</text>
  <text x="440" y="98" text-anchor="middle" font-size="26" font-family="sans-serif" fill="#777777">&#8230;</text>
  <rect x="500" y="60" width="140" height="55" rx="8" fill="#eef7ee" stroke="#2e7d32" stroke-width="2"/>
  <text x="570" y="83" text-anchor="middle" font-size="13" font-family="sans-serif" fill="#2e7d32">Age 18 (70-90y)</text>
  <text x="570" y="100" text-anchor="middle" font-size="11" font-family="sans-serif" fill="#2e7d32">ladder inside</text>
  <path d="M-40,87 L8,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#b1)"/>
  <text x="-38" y="75" text-anchor="start" font-size="11" font-family="sans-serif">births: &#945; N_total</text>
  <path d="M150,87 L218,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#b1)"/>
  <text x="185" y="72" text-anchor="middle" font-size="10" font-family="sans-serif">ageing: 1/w_1</text>
  <path d="M360,87 L498,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#b1)"/>
  <path d="M700,140 L700,160" fill="none" stroke="#7a1f1f" stroke-width="2" marker-end="url(#b1)"/>
  <text x="700" y="175" text-anchor="middle" font-size="10" font-family="sans-serif" fill="#7a1f1f">exits (death, all ages)</text>
</svg>

$\varepsilon$ (age-assortative mixing) controls how much contact - and
therefore infection - happens within your own age band versus randomly
across the whole population: $\varepsilon=0$ is fully random mixing,
$\varepsilon=1$ is entirely within-age contact only.

### Equations

Within age group $a$, rung $j$ (mirroring `teaching_app_ladder`, but now
per age group and including ageing and mortality):

$$
\lambda_a = \beta \sum_{b=1}^{18} w_{ab} \frac{\sum_j l_j I_{b,j}}{N_b}
\qquad
w_{ab} = (1-\varepsilon)\frac{N_b}{N_{\text{tot}}} + \varepsilon\,\mathbb{1}[a=b]
$$

$$
\frac{dS_{a,j}}{dt} = (\text{ageing in}) - \frac{S_{a,j}}{w_a} - \lambda_a S_{a,j} + (\text{recovery from rung } j-1) - \mu_a S_{a,j}
$$

$$
\frac{dI_{a,j}}{dt} = (\text{ageing in}) - \frac{I_{a,j}}{w_a} + \lambda_a S_{a,j} - \nu_j I_{a,j} - \mu_a I_{a,j}
$$

where $w_a$ is the width (years) of age band $a$ and $\mu_a$ its
age-specific mortality rate. Births ($\alpha N_{\text{tot}}$) enter
$S_{1,1}$ (youngest age group, rung 1). $\nu_j$ and $l_j$ are exactly as in
`teaching_app_ladder`:

$$
\nu_j = (\nu_0 - \nu_1) e^{-\text{nuExponent}(j-1)} + \nu_1
\qquad\qquad
l_j = e^{-\text{lExponent}(j-1)}
$$

TS/TT accumulate from repeat-infection burden exactly as before, but now
per age group, with ageing and mortality added the same way as S/I.

| Symbol | Meaning |
|---|---|
| $\beta$ (transRate) | transmission intensity |
| $\varepsilon$ | age-assortative mixing (0 = random, 1 = within-age only) |
| $\nu_0,\nu_1,$ nu-exponent | recovery-rate-by-rung parameters |
| l-exponent | infectiousness-decay-by-rung parameter |
| $\alpha$ | birth rate (scales total population) |
| k_TS, k_TT | scarring / trichiasis accumulation rates |
| $w_a$, $\mu_a$ | age band width and mortality (fixed, not user-adjustable here) |

### Reading the plots

- **Infection prevalence over time** aggregates across *all* ages and
  rungs into one overall trace - useful for seeing how long the system
  takes to approach its equilibrium, and the effect of any MDA rounds
  (dashed lines).
- **Equilibrium infection prevalence by age** is a bar chart of the
  final (tail-window-averaged) age-specific prevalence - this is the
  classic trachoma signature of a peak in young children (ages 1-9) that
  falls with age as partial immunity accumulates, driven by the age
  structure and the ladder together.
- **Equilibrium disease burden by age** shows scarring (TS) and trichiasis
  (TT) prevalence by age band - since these accumulate irreversibly, they
  rise monotonically with age even as infection prevalence itself falls.

This app can be slow to update for large numbers of infection-history
rungs or long simulation horizons, since the underlying system has
$18 \times (2n+2)$ state variables where $n$ is the number of rungs.

See `trachoma_first_principles_300_history.R` in this repository for the
underlying model code (used here directly, not reimplemented), and
`teaching_app_ladder` for the same ladder mechanics without age structure.
