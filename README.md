# Stage-structured-Jolly-Seber-model

Objective: To estimate survival and population growth for a population of marked individuals that differ in infection status

Code: Learn to simulate and analyze capture-mark-recapture data for a stage-structured (i.e., age, size, sex, or any other discrete category) population

In this repositoru, we will simulate data for a capture-mark-recapture study, and we will analyze that data in JAGS.
 
Here is the scenario: We are disease ecologists, and we want to understand the disease dynamics of our system. When we capture an individual, we assign it to a certain ‘state.’ A state in our case is a disease status, but can also be defined as a geographic location, an age, a breeding status, etc.
 
 
Before we begin, I’d like to lay out some information:
 
 
Our parameters of interest: state-specific apparent survival probability, recruitment probability, transition rates (i.e., movement between two ‘states’), abundance for each state, and population growth rate.


The type of data we ‘collect’: capture-mark-recapture- where each individual is captured, uniquely marked (e.g., elastomer tags, pit tag, toe clips, etc.), and released. Then, each season, new individuals are marked, and old individuals are recaptured. When an individual is captured, we record its state.

 
Our modeling approach: Multi-state Jolly-Seber Model

 
Statistical framework: Bayesian

 
For more information on these models, I highly recommend the book: Bayesian population analysis using WinBUGs by Marc Kéry and Michael Schaub.

