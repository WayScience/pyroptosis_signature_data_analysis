## Linear Modeling of CellProfiler features to understand differences between treatments and control.

Fit a series of linear models (one per feature) to determine which features contribute the most to the differences between treatments and controls (highest beta coefficients) adjusting for cell count per well.

A couple of models are used to determine feature variance and the condition's contribution to the variance.

files
`1a.fit_linear_model_1beta.ipynb`
`1b.fit_linear_model_2beta.ipynb`
`1c.fit_linear_model_3beta.ipynb`
`1d.fit_linear_model_4beta.ipynb`

fit a different model to the same data

and  files


`2a.Linear_Modeling_Visualization.ipynb`
`2b.Linear_Modeling_Visualization.ipynb`
`2c.Linear_Modeling_Visualization.ipynb`
`2d.Linear_Modeling_Visualization.ipynb`

to visualize each model into one pdf

#### Only for the single beta + cell count models:
I have included a visualization from [CFReT Analysis](https://github.com/gwaybio/CFReT_data/blob/368f0baf7b027f5614de7173bbc7336add81fcb3/4.analyze-data/notebooks/linear_model/1.visualize_linear_model.ipynb) to look at the beta coefficients across morphology feature types and cellular compartments:
* Feature Types:
    * Texture
    * Radial Distribution
    * Intensity
    * Granularity
    * Correlation
* Cellular Compartments:
    * Nucleus
    * Cytoplasm
    * Cell

For more information about the features, please see the [CellProfiler documentation](https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-4.2.6/modules/measurement.html).

## The Models:

for each linear model,
$y$ is each feature
$x$ is the inputed variables
$\epsilon$ is the residual variance not explained by factors in the model

## Linear Model a:
$y = \beta _{0}x+ \beta _{1}x+ \epsilon$ \
where;
$\beta _{0}$ is the beta coefficient attributed to cell count.
$\beta _{1}$ is the beta coefficient attributed to Inducer, Inducer Dose, Inhibitor, and Inhibitor Dose.


## Linear Model b:
$y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \epsilon$ \
where;
$\beta _{0}$ is the beta coefficient attributed to cell count.
$\beta _{1}$ is the beta coefficient attributed to Inducer, Inhibitor, and Inhibitor Dose.
$\beta _{2}$ is the beta coefficient attributed to Inducer dose.

## Linear Model c:
$y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \beta _{3}x+ \epsilon$ \
where;
$\beta _{0}$ is the beta coefficient attributed to cell count.
$\beta _{1}$ is the beta coefficient attributed to Inducer.
$\beta _{2}$ is the beta coefficient attributed to Inducer dose.
$\beta _{3}$ is the beta coefficient attributed to Inhibitor, and Inhibitor Dose.

## Linear Model d:
$y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \beta _{3}x+ \beta _{4}x+ \epsilon$ \
where;
$\beta _{0}$ is the beta coefficient attributed to cell count.
$\beta _{1}$ is the beta coefficient attributed to Inducer.
$\beta _{2}$ is the beta coefficient attributed to Inducer dose.
$\beta _{3}$ is the beta coefficient attributed to Inhibitor.
$\beta _{4}$ is the beta coefficient attributed to Inhibitor Dose.

