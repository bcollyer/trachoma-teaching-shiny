## What is this model?

This app is a **teaching tool** for exploring the transmission dynamics of
ocular *Chlamydia trachomatis* (trachoma) using a simple **SIS
(Susceptible-Infected-Susceptible)** compartmental model.

Trachoma is modelled as SIS rather than SIR because infection does **not**
confer lasting protective immunity. Longitudinal field studies (e.g. Bailey
et al. 1999; Grassly et al. 2008) show that individuals are repeatedly
re-infected throughout life.

This is a **simplified teaching model** — it has **no age structure** and no **partial immunity**.

### Compartments

- **S** — susceptible individuals (currently uninfected, can become infected)
- **I** — infected individuals (currently infected, infectious to others)

The population size $N = S + I$ is held constant (no births/deaths in this
simplified teaching version).

### Diagram

<svg width="100%" height="230" viewBox="0 0 560 230" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Compartment diagram: S and I boxes with arrows for infection, recovery and MDA">
  <defs>
    <marker id="arrowBlack" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#333333"/>
    </marker>
    <marker id="arrowGreen" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#1a7f37"/>
    </marker>
  </defs>

  <!-- S box -->
  <rect x="50" y="95" width="120" height="60" rx="8" fill="#eaf2fb" stroke="#2b6cb0" stroke-width="2"/>
  <text x="110" y="132" text-anchor="middle" font-size="22" font-family="sans-serif" fill="#2b6cb0">S</text>

  <!-- I box -->
  <rect x="390" y="95" width="120" height="60" rx="8" fill="#fdecea" stroke="#c0392b" stroke-width="2"/>
  <text x="450" y="132" text-anchor="middle" font-size="22" font-family="sans-serif" fill="#c0392b">I</text>

  <!-- S -> I: force of infection -->
  <path d="M170,105 C 250,70 350,70 390,105" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#arrowBlack)"/>
  <text x="280" y="65" text-anchor="middle" font-size="14" font-family="sans-serif">force of infection: &#946; I / N</text>

  <!-- I -> S: spontaneous recovery -->
  <path d="M390,145 C 350,185 250,185 170,145" fill="none" stroke="#333333" stroke-width="2" marker-end="url(#arrowBlack)"/>
  <text x="280" y="205" text-anchor="middle" font-size="14" font-family="sans-serif">natural clearance: &#947;</text>

  <!-- I -> S: MDA treatment (dashed, higher arc) -->
  <path d="M450,95 C 450,20 110,20 110,95" fill="none" stroke="#1a7f37" stroke-width="2" stroke-dasharray="6,4" marker-end="url(#arrowGreen)"/>
  <text x="280" y="15" text-anchor="middle" font-size="14" font-family="sans-serif" fill="#1a7f37">MDA treatment (coverage per round)</text>
</svg>

**S** flows to **I** via the force of infection ($\beta I/N$); **I** flows
back to **S** either through natural clearance ($\gamma$, always on) or
through periodic MDA rounds (dashed green arrow, only at the round times you
set with the sliders).

### Equations

$$
\frac{dS}{dt} = -\beta \frac{I}{N} S + \gamma I
$$

$$
\frac{dI}{dt} = \beta \frac{I}{N} S - \gamma I
$$

| Symbol | Meaning | Units |
|---|---|---|
| $\beta$ | effective transmission rate | per year |
| $\gamma$ | recovery / spontaneous clearance rate | per year |
| $N$ | population size | people |
| $R_0 = \beta/\gamma$ | basic reproduction number | dimensionless |

### Equilibrium prevalence

Ignoring MDA, this SIS system has a stable **endemic equilibrium** whenever
$R_0 > 1$, at prevalence:

$$
I^*/N = 1 - \frac{1}{R_0} = 1 - \frac{\gamma}{\beta}
$$

If $R_0 \le 1$ the infection cannot sustain itself and dies out
($I \to 0$).

### Mass Drug Administration (MDA)

Trachoma control programmes rely on periodic **Mass Drug Administration**
(azithromycin) rather than a permanent cure for susceptibility. In this app,
each MDA round instantaneously treats a fixed **coverage** fraction of
currently infected people, moving them from **I** back to **S** — but
because there is no immunity, treated individuals can be reinfected again
before the next round. This is why repeated rounds (not a single round) are
usually needed to interrupt transmission, and why coverage and round
frequency both matter for elimination.

### How to use the simulator tab

Adjust the sliders and watch:

- how **prevalence** ($I/N$) rises to (or falls from) its equilibrium value;
- how increasing $\beta$ (more transmission) or decreasing $\gamma$ (slower
  clearance) raises $R_0$ and the equilibrium prevalence;
- how adding MDA rounds (coverage and interval) can push prevalence down,
  and what happens to prevalence once MDA stops (it rebounds towards the
  endemic equilibrium if $R_0$ is still above 1).


