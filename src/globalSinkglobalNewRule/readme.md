### Global New Rules, and a Global Sink

The model is as follows:

We simulate a 2D grid of size `LxL`, where each site contains molecules of various chemicals (`0` to `N_CHEMICALS`, with `0` being empty space). Each site has a total of `MOLECULES_PER_SITE` molecules, which is conserved (note: not the mass, but the number of molecules).

Chemicals are represented by integers from `0` (empty) to `N_CHEMICALS`. Each chemical has a "mass" equal to its number and a "volatility" randomly assigned at initialization. Chemicals interact through "rules": chemical reactions of the form `A + B + C -> D + E + F`. Rules must conserve or decrease mass (decreasing mass acts like a sink). Chemical `A` (with a mass of 1) is supplied form the top row of the lattice.

The crucial aspect is **rule addition**. The rules (reactions) are not set in stone - instead, they are randomly chosen based on the concentrations currently in the lattice. If an excess of molecule `D` is built up, a rule is chosen that converts `D` into something else. This "evolution" of rules may lead to more complex behaviour than a simple CRN.

### Algorithm

At each timestep, the following steps are performed in sequence:

1. **Supply** (parallel): Convert up to `SUPPLY_RATE` empty molecules to chemical 1 at each site in the top row of the lattice.

2. **Reactions** (parallel): For each site, identify all eligible rules where sufficient reactants are present. Compute a weight for each eligible rule as the product of `(concentration[c] * volatility[c])` for each reactant `c` (accounting for multiplicity). Select one rule proportionally to its weight and apply it once, consuming reactants and producing products.

3. **Diffusion** (parallel): Perform molecule swaps between neighboring sites in a 6-phase cyclic pattern to ensure no race conditions. In other words, the lattice is labelled into sites like `a,b,c,d,e,f`, in the following pattern:
    ```
    a1 b1 c1 d1 e1 f1 a2 b2 c2 d2 e2 f2
    c3 d3 e3 f3 a3 b3 c4 d4 e4 f4 a4 b4
    e5 f5 a5 b5 c5 d5 e6 f6 a6 b6 c6 d6
    ```
    Thus the neighbours of `a`s is:
    ```
    a1: b1, c3, 
    a2: f1, b2, c4, 
    a3: e1, f3, b3, c5
    a4: e2, f4, b4, c6
    a5: e3, f5, b5
    a6: e4, f6, b6
    ```
    Which ensures no overlap. We first parallely swap for all `A`s, then all `b`s, etc
    Each site attempts `D` swaps per phase, exchanging molecules if both sites have the required chemicals.

4. **Rule Addition** (sequential): Count global chemical abundances across the lattice. Sample three reactants based on global counts weighted by volatilities. If no matching rule exists, generate random products that conserve or decrease mass, and add the new rule if valid.

We perform the algorithm for `N_STEPS` and periodically record both the state of the lattice as well as the rules.
