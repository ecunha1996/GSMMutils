import cobra


class CobrapyDFBA:
    def __init__(self, model):
        self.model = model
        self.pbar = None
        self.epsilon = 1E-4
        self.direction = 1
        self.terminal = True

    def add_dynamic_bounds(self, y):
        """Use external concentrations to bound the uptake flux of glucose."""
        biomass, phosphate = y  # expand the boundary species
        phospahte_max_import = -0.04*phosphate/(0.00185+phosphate)
        self.model.reactions.EX_C00009__dra.lower_bound = phospahte_max_import

    def dynamic_system(self, t, y):
        """Calculate the time derivative of external species."""

        biomass, phosphate = y  # expand the boundary species

        # Calculate the specific exchanges fluxes at the given external concentrations.
        with self.model:
            self.add_dynamic_bounds(y)

            cobra.util.add_lp_feasibility(self.model)
            feasibility = cobra.util.fix_objective_as_constraint(self.model)
            lex_constraints = cobra.util.add_lexicographic_constraints(self.model, ['e_ActiveBiomass__cytop', 'EX_C00009__dra'], ['max', 'max'])

        # Since the calculated fluxes are specific rates, we multiply them by the
        # biomass concentration to get the bulk exchange rates.
        fluxes = lex_constraints.values
        fluxes *= biomass

        # This implementation is **not** efficient, so I display the current
        # simulation time using a progress bar.
        if self.pbar is not None:
            self.pbar.update(1)
            self.pbar.set_description('t = {:.3f}'.format(t))
        return fluxes

    def infeasible_event(self, t, y):
        """
        Determine solution feasibility.

        Avoiding infeasible solutions is handled by solve_ivp's built-in event detection.
        This function re-solves the LP to determine whether or not the solution is feasible
        (and if not, how far it is from feasibility). When the sign of this function changes
        from -epsilon to positive, we know the solution is no longer feasible.

        """

        with self.model:
            self.add_dynamic_bounds(y)
            cobra.util.add_lp_feasibility(self.model)
            feasibility = cobra.util.fix_objective_as_constraint(self.model)

        return feasibility - self.epsilon