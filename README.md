This Shiny app uses the `simID` function from the `SemiCompRisks` package to simulate data from an illness-death model. The hazard models are parametric (Weibull) with a semi-Markov dependence between the non-terminal event and the hazard of the terminal event experienced after the non-terminal event. This first tab gives a description of each parameter's role in the data generating process.

The pdf of the terminal event and the mixed pdf/pmf of the non-terminal event is shown in the second tab.

The time-varying composition of the survival-determined principal strata of "Always-alive," "Control-killed," "Treatment-killed," and "Doubly-dead" is shown in the third tab. 

The fourth tab shows one potential causal effect of interest, a difference in cumulative incidence of the non-terminal event among those whose principal state at that point in time was "always-alive."
