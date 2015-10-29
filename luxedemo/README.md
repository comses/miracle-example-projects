## Archive for the luxedemo project

This folder contains
* A small subset of `luxe` output data
* A demo Shiny app that produces an interactive version of Figure 1 in Huang et al. (2013)
* Two demo R scripts that produce parameterized versions of Figure 1 and 8 in Huang et al. (2013)

### apps/shiny-demo

This Shiny app produces an interactive version of Figure 1 in Huang et al. (2013): "Effects of agent heterogeneity in the presence of a land-market: A systematic test in an agent-based laboratory". In the model, agents select a land parcel based on its proximity to the city center and its open space amenity. `sdp` is the standard deviation of buyers' preferences for proximity to city center. `sdbudget` is the standard deviation of buyers' budgets. A higher value of these variables means increasing level of heterogeneity among agents in terms of preferences or budgets, while a zero value means agents are identical to each other in terms of preference (or budget, or both).

### src/Figure1.R

This script produces a parameterized version of Figure 1 in Huang et al. (2013). In the model, agents select a land parcel based on its proximity to the city center and its open space amenity. `sdp2` is the standard deviation of buyers' preferences for proximity to city center. `sdb3` is the standard deviation of buyers' budgets. A higher value of these variables means increasing level of heterogeneity among agents in terms of preferences or budgets, while a zero value means agents are identical to each other in terms of preference (or budget, or both).

### src/Figure8.R

This script produces a parameterized version of Figure 8 in Huang et al. (2013). The output figure is similar to Figure 1 but in 3D to show the bid-rent curve and sequence of development. Height of the bars shows the transaction prices, while colour of the bars shows the development sequence (red ones are developed first, while blue ones last). It takes the same two parameters as `Figure1.R`: `sdp2` is the standard deviation of buyers' preferences for proximity to city center. `sdb3` is the standard deviation of buyers' budgets. A higher value of these variables means increasing level of heterogeneity among agents in terms of preferences or budgets, while a zero value means agents are identical to each other in terms of preference (or budget, or both).