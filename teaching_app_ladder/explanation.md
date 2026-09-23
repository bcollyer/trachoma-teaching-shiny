## What is this model?

This app extends the basic SIS trachoma model (see the companion
`teaching_app`) with a **partial-immunity "ladder"** of infection-history
rungs, using the same continuous formulation as the full model in
`trachoma_first_principles_300_history.R` (which uses up to 100+ rungs per
age group) — but with a **configurable, small number of rungs** (default
5) and **no age structure** (a single well-mixed population) so the ladder
mechanics are easier to see on their own.

### Why a ladder?

Four features of trachoma biology motivate the extra compartments:

1. **No lasting immunity (SIS, not SIR)** — people are reinfected
   throughout life; there is no permanent recovered state.
2. **Partial acquired immunity builds up gradually** — rather than a single
   jump from "naive" to "immune", each completed infection episode nudges a
   person one rung further along a ladder. Immunity here shows up as
   **faster recovery** at higher rungs (governed by $\nu_0,\nu_1,$
   nu-exponent) and **lower infectiousness** at higher rungs (governed by
   the l-exponent) — exactly the mechanism used in the full model, just
   without an age dimension.
3. **Cumulative, irreversible disease** — scarring (**TS**) and trichiasis
   (**TT**) accumulate with repeat infection and do not resolve on their
   own (only **TT** can be resolved here, e.g. by surgery).
4. **No age structure** — instead of age groups, this app uses a single
   background birth/death rate $\mu = 1/\text{life expectancy}$ so that the
   naive-susceptible rung (**S_1**) keeps being replenished (playing the
   role that "newborns" would in the full age-structured model).

### Compartments

- **S_1 ... S_n** — susceptible, by number of past infections (S_1 = never
  infected; S_n = has recovered from many episodes)
- **I_1 ... I_n** — currently infected, by number of past infections
- **TS** — trachomatous scarring (irreversible)
- **TT** — trachomatous trichiasis (irreversible, but can be resolved by
  surgery/treatment)

$\sum_j (S_j + I_j) = N$ is held constant (births balance deaths). **TS**
and **TT** are a separate disease-burden ladder driven by repeat infection
(rungs $j \ge 2$) and do not feed back into transmission.

### How a person moves through the ladder

- Getting infected moves you from $S_j$ to $I_j$ (same rung).
- Recovering from $I_j$ moves you **up** to $S_{j+1}$ — except recovery
  from the **top** rung $I_n$, which loops back to $S_n$ (there is nowhere
  higher to go, so the ladder is "capped").
- **MDA** treatment instead cures $I_j \to S_j$ at the **same** rung, for
  every rung — matching `make_mda_events()` in the full model.

### Diagram

<svg width="100%" height="420" viewBox="-40 0 1080 420" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Ladder diagram: S_1, I_1, S_2, I_2, further rungs, S_n, I_n, TS, TT boxes with arrows for infection, immunity progression, MDA, and irreversible disease progression driven by every repeat-infection rung I_2 through I_n">
  <defs>
    <marker id="lArrowBlack" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#333333"/>
    </marker>
    <marker id="lArrowGrey" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#777777"/>
    </marker>
    <marker id="lArrowGreen" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#1a7f37"/>
    </marker>
  </defs>

  <!-- S1 -->
  <rect x="15" y="60" width="100" height="55" rx="8" fill="#eaf2fb" stroke="#2b6cb0" stroke-width="2"/>
  <text x="65" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#2b6cb0">S_1</text>
  <!-- I1 -->
  <rect x="165" y="60" width="100" height="55" rx="8" fill="#fdecea" stroke="#c0392b" stroke-width="2"/>
  <text x="215" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#c0392b">I_1</text>
  <!-- S2 -->
  <rect x="315" y="60" width="100" height="55" rx="8" fill="#f1eafb" stroke="#6b46c1" stroke-width="2"/>
  <text x="365" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#6b46c1">S_2</text>
  <!-- I2 -->
  <rect x="465" y="60" width="100" height="55" rx="8" fill="#fdf0e0" stroke="#b7791f" stroke-width="2"/>
  <text x="515" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#b7791f">I_2</text>
  <!-- ellipsis -->
  <text x="632" y="98" text-anchor="middle" font-size="26" font-family="sans-serif" fill="#777777">&#8230;</text>
  <!-- Sn -->
  <rect x="700" y="60" width="100" height="55" rx="8" fill="#f1eafb" stroke="#6b46c1" stroke-width="2"/>
  <text x="750" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#6b46c1">S_n</text>
  <!-- In -->
  <rect x="850" y="60" width="100" height="55" rx="8" fill="#fdf0e0" stroke="#b7791f" stroke-width="2"/>
  <text x="900" y="93" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#b7791f">I_n</text>
  <!-- TS -->
  <rect x="850" y="220" width="100" height="55" rx="8" fill="#f0f0f0" stroke="#555555" stroke-width="2"/>
  <text x="900" y="253" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#555555">TS</text>
  <!-- TT -->
  <rect x="850" y="330" width="100" height="55" rx="8" fill="#f8e0e0" stroke="#7a1f1f" stroke-width="2"/>
  <text x="900" y="363" text-anchor="middle" font-size="16" font-family="sans-serif" fill="#7a1f1f">TT</text>

  <!-- births into S1 -->
  <path d="M-30,87 L13,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="-38" y="80" text-anchor="start" font-size="12" font-family="sans-serif">births: &#956;N</text>

  <!-- S1 -> I1 -->
  <path d="M115,87 L163,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="140" y="72" text-anchor="middle" font-size="12" font-family="sans-serif">&#955;</text>

  <!-- I1 -> S2 -->
  <path d="M265,87 L313,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="290" y="72" text-anchor="middle" font-size="11" font-family="sans-serif">&#957;1 (up ladder)</text>

  <!-- S2 -> I2 -->
  <path d="M415,87 L463,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="440" y="72" text-anchor="middle" font-size="12" font-family="sans-serif">&#955;</text>

  <!-- I2 -> ... -> Sn -->
  <path d="M565,100 C 620,125 660,125 698,100" fill="none" stroke="#777777" stroke-width="2" stroke-dasharray="3,4" marker-end="url(#lArrowGrey)"/>
  <text x="630" y="140" text-anchor="middle" font-size="11" font-family="sans-serif" fill="#777777">(rungs 3 .. n-1)</text>

  <!-- Sn -> In -->
  <path d="M800,87 L848,87" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="825" y="72" text-anchor="middle" font-size="12" font-family="sans-serif">&#955;</text>

  <!-- In -> Sn (capped, SIS cycle at top rung) -->
  <path d="M850,110 C 820,150 780,150 750,117" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="800" y="160" text-anchor="middle" font-size="11" font-family="sans-serif">&#957;n (capped: SIS cycle)</text>

  <!-- all repeat-infection rungs (I_2 .. I_n, including the omitted ones) feed the scarring accumulator -->
  <path d="M515,115 L515,195" fill="none" stroke="#777777" stroke-width="2" stroke-dasharray="3,4"/>
  <path d="M632,133 L632,195" fill="none" stroke="#777777" stroke-width="2" stroke-dasharray="3,4"/>
  <path d="M900,117 L900,195" fill="none" stroke="#777777" stroke-width="2" stroke-dasharray="3,4"/>
  <path d="M515,195 L900,195" fill="none" stroke="#777777" stroke-width="2"/>
  <text x="700" y="212" text-anchor="middle" font-size="11" font-family="sans-serif" fill="#777777">every repeat-infection rung I_2 .. I_n feeds scarring</text>

  <!-- collected repeat-infection burden -> TS -->
  <path d="M900,195 L900,218" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="960" y="210" text-anchor="middle" font-size="12" font-family="sans-serif">r_TS</text>

  <!-- TS -> TT -->
  <path d="M900,277 L900,328" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="960" y="305" text-anchor="middle" font-size="12" font-family="sans-serif">r_TT</text>

  <!-- TT resolution -->
  <path d="M957,357 L1000,357" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#lArrowBlack)"/>
  <text x="965" y="345" text-anchor="start" font-size="12" font-family="sans-serif">&#945; (surgery)</text>

  <!-- MDA, applies to every rung -->
  <path d="M215,60 C 215,20 65,20 65,60" fill="none" stroke="#1a7f37" stroke-width="2" stroke-dasharray="6,4" marker-end="url(#lArrowGreen)"/>
  <text x="140" y="15" text-anchor="middle" font-size="11" font-family="sans-serif" fill="#1a7f37">MDA: cures I_j&#8594;S_j (every rung j)</text>
</svg>

Background births/deaths at rate $\mu = 1/\text{life expectancy}$ apply to
every rung (not drawn, to keep the diagram readable) — this is what keeps
a supply of naive susceptibles (**S_1**) available, playing the role that
ageing into the population does in the full age-structured model.

The dashed grey lines converging below the top row show that **every**
repeat-infection rung — $I_2$, the omitted rungs 3..n-1, and $I_n$ — feeds
the scarring accumulator, not just the last rung. First-ever infection
($I_1$) is the only rung that does **not** contribute to TS.

### Equations

For rung $j = 1, \dots, n$, define the rung-specific recovery rate and
relative infectiousness exactly as in the full model:

$$
\nu_j = (\nu_0 - \nu_1)\, e^{-\text{nuExponent} \cdot (j-1)} + \nu_1
\qquad\qquad
l_j = e^{-\text{lExponent} \cdot (j-1)}
$$

so recovery gets **faster** and infectiousness gets **weaker** at higher
rungs. The shared force of infection (same for every susceptible rung — the
ladder does not reduce *susceptibility*, only infectiousness and recovery
speed) is:

$$
\lambda = \beta \frac{\sum_{j=1}^{n} l_j I_j}{N}
$$

$$
\frac{dS_1}{dt} = \mu N - \lambda S_1 - \mu S_1
\qquad
\frac{dI_1}{dt} = \lambda S_1 - (\nu_1 + \mu) I_1
$$

$$
\frac{dS_j}{dt} = \nu_{j-1} I_{j-1} - \lambda S_j - \mu S_j
\qquad (1 < j < n)
$$

$$
\frac{dS_n}{dt} = \nu_{n-1} I_{n-1} + \nu_n I_n - \lambda S_n - \mu S_n
\qquad \text{(top rung is capped)}
$$

$$
\frac{dI_j}{dt} = \lambda S_j - (\nu_j + \mu) I_j
\qquad (1 < j \le n)
$$

$$
\frac{dT_S}{dt} = r_{TS} \frac{\sum_{j \ge 2} I_j}{N} (N - T_S) - \mu T_S
\qquad\qquad
\frac{dT_T}{dt} = r_{TT} \frac{T_S}{N} (N - T_T) - (\alpha + \mu) T_T
$$

| Symbol | Meaning | Units |
|---|---|---|
| $\beta$ | transmission intensity | per year |
| $\nu_0$ | recovery rate at rung 1 (first infection) | per year |
| $\nu_1$ | asymptotic recovery rate at high rungs | per year |
| nu-exponent | speed of the $\nu_0 \to \nu_1$ transition across rungs | dimensionless |
| l-exponent | decline in relative infectiousness across rungs | dimensionless |
| $n$ | number of infection-history rungs | count |
| $r_{TS}$ | rate of scarring accumulation, driven by repeat infection ($j\ge 2$) | per year |
| $r_{TT}$ | rate of progression from scarring to trichiasis | per year |
| $\alpha$ | rate of trichiasis resolution (e.g. surgery) | per year |
| $\mu$ | background birth/death rate $= 1/\text{life expectancy}$ | per year |

### Reading the plots

- **Infection prevalence** (all infected rungs combined, $\sum_j I_j / N$)
  shows the overall proportion currently infected.
- **Total susceptible vs total infected** shows the aggregate split.
- **Composition of infected rungs** is a stacked area of $I_1, \dots, I_n$
  as a proportion of $N$ over time — watch the mix shift from mostly
  first-time infections ($I_1$) early on towards a mix dominated by
  higher, less-infectious rungs once the ladder has had time to fill up.
- **TS and TT** show the slowly-accumulating, irreversible disease burden.
  Because these only depend on cumulative exposure via repeat infections,
  they keep rising long after infection prevalence has reached its steady
  state, and only fall when MDA (via reduced repeat infection) or
  trichiasis resolution ($\alpha$) outweigh new accumulation.

### Basic reproduction number

The displayed $R_0$ is only the **invasion** number — the expected number
of secondary infections from a single case introduced into a fully naive
(all-S_1) population, $R_0 = \beta / (\nu_0 + \mu)$ (using $l_1 = 1$).
Once the ladder fills up, the effective reproduction number depends on the
mix of rungs and no longer has this simple closed form — use the simulator
to see the resulting dynamics directly rather than relying on a single
number.

This is still a **simplified teaching model**: no age structure, a single
averaged life expectancy standing in for real demography, and typically far
fewer rungs than the full model. 
