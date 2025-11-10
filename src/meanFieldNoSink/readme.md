# Mean-field

This code implements a mean-field version of the model. It's quite simple: we have two parameters:

- The number of molecules in the system (including empty space)
- The 'supply rate'

Each timestep, we do the following:
- Supply: Try to make `SUPPLY_RATE` `0`s (empty sites) into `1`s (adding mass to the system). If we have fewer zeros, all of them are converted instead.
- React: Choose 3 random molecules from the system (which could include `0`s). If a rule exists for them, implement it, and convert as much as possible according to the rule. If no rule exists, create one, and react accordingly.