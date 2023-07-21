
All code in this directory is inspired or taken from the following repository: https://github.com/WayScience/phenotypic_profiling_model/tree/main

The goal in this directory is three fold on two different resolutions:
1. To create a logistic regression model that can predict the probability of a treatment that correlates with a specific cell state using the morphological features of the cells.
    * At the well level
    * At the single cell level
2. To create a logistic regression model that can predict the probability of a treatment that correlates with a specific cell state using cytokine and chemokine quantification data.
    * At the well level.
3. To create a logistic regression model that can predict the probability of a treatment that correlates with a specific cell state using the morphological features of the cells and cytokine and chemokine quantification data.
    * At the well level.
    * At the single cell level.
        * This will be harder to accomplish as the cytokine and chemokine quantification data is only available at the well level.
