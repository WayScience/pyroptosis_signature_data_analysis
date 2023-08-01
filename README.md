# pyroptosis_signature_data_analysis
Analysis Repository for the Interstellar project. For data processing for this project please see [pyroptosis_signature](https://github.com/WayScience/pyroptosis_signature)

### Workflow:

Cells were incubated for 1h in the inhibitors, and then incubated for 6h in each inducer.
For more information on plate layout and experimental design please see [pyroptosis_signature](https://github.com/WayScience/pyroptosis_signature).


# Treatments of Interest

| Inducer           | Number of Doses   | Function                      |
| :----------       | :---------------: | :---------                    |
| Media             | 1                 | Control                       |
| DMSO              | 1                 | Vehicle Control               |
| Thapsigargin      | 2                 | ER Stress Induced Apoptosis   |
| H202              | 2                 | Cell Death                    |
| Flagellin         | 3                 | Induce Pyroptosis             |
| LPS               | 5                 | Induce Pyroptosis             |
| LPS + Nigericin   | 2 * 3 = (6)       | Induce Pyroptosis             |
| Disulfiram        | 3                 | Pyroptosis Inhibitor          |


# Other Treatments and combinations

| Inducers       | # of Doses  | Dose Amount(s)                                                        | Inhibitors | # of Doses | Function                      |
| :----------    | :---------: | :---------                                                            | :--------  | :--------: | :---------                    |
| Media          | 1           | 0                                                                     | Media      | 1          | Control                       |
| DMSO           | 1           | 0.0250%                                                               | DMSO       | 2          | Vehicle Control               |
| DMSO           | 1           | 0.0250%                                                               | Z-VAD-FMK  | 2          | Inhibit Inflamasome           |
| Thapsigargin   | 2           | 1uM <br> 10uM                                                         | DMSO       | 1          | ER Stress Induced Apoptosis   |
| H202           | 2           | 100nM <br> 100uM                                                      | DMSO       | 1          | Cell Death                    |
| H202           | 2           | 100nM <br> 100uM                                                      | Z-VAD-FMK  | 1          | Inhibit Inflamasome           |
| H202           | 2           | 100nM <br> 100uM                                                      | Disulfiram | 1          | Cell Death                    |
| Flagellin      | 2           | 0.1ug/mL <br> 1ug/mL                                                  | DMSO       | 1          | Induce Pyroptosis             |
| Flagellin      | 2           | 0.1ug/mL <br> 1ug/mL                                                  | Disulfiram | 1          | Induce Pyroptosis <br> Inhibit Inflamasome          |
| LPS            | 5           | 0.01ug/mL <br> 0.1ug/mL <br> 1.0ug/mL <br> 10.0ug/mL <br> 100.0 ug/mL | DMSO       | 1          | Induce Pyroptosis             |
| LPS            | 1           | 10.0ug/mL                                                             | Z-VAD-FMK  | 1          | Induce Pyroptosis / Inhibit Inflamasome          |
| LPS            | 1           | 10.0ug/mL                                                             | Disulfiram | 3          | Induce Pyroptosis / Inhibit Inflamasome      |
| LPS + Nigercin | 2 * 3 = (6) | 1ug/mL+1ug/mL <br> 1ug/mL+3ug/mL <br> 1ug/mL+10ug/mL <br> 100ug/mL+1ug/mL <br> 100ug/mL+3ug/mL <br> 100ug/mL+10ug/mL | DMSO | 1 | Induce Pyroptosis |
| LPS + Nigercin | 1 * 1 = (1) | 1ug/mL+10ug/mL         | Z-VAD-FMK  | 1               | Induce Pyroptosis <br> Inhibit Inflamasome          |
| LPS + Nigercin | 1 * 1 = (1) | 1ug/mL+10ug/mL        | Disulfiram | 1               | Induce Pyroptosis <br> Inhibit Inflamasome          |
| Disulfiram     | 3           | 0.1uM <br> 1uM <br> 2.5uM       | DMSO       | 1               | Pyroptosis Inhibitor          |


