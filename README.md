# pyroptosis_signature_data_analysis
Analysis Repository for the Interstellar project. For data processing for this project please see [pyroptosis_signature](https://github.com/WayScience/pyroptosis_signature)

### Workflow:

Cells were incubated for 1h in the inhibitors, and then incubated for 6h in each inducer.
For more information on plate layout and experimental design please see [pyroptosis_signature](https://github.com/WayScience/pyroptosis_signature).


# Treatments of Interest (All in this table are with DMSO 0.025% inhibitor)

| Inducer           | Number of Doses   | Dose Concentration(s) | Function                      |
| :----------       | :---------------: | :-----| :---------                    |
| DMSO              | 1                 | 0.025%| Vehicle Control               |
| Thapsigargin      | 2                 | 1uM <br> 10uM | ER Stress Induced Apoptosis   |
| H202              | 2                 | 100nM <br> 100uM | Cell Death                    |
| Flagellin         | 2                 | 0.1ug/mL <br> 1.0ug/mL | Induce Pyroptosis             |
| LPS               | 3                 | 0.01ug/mL <br> 1.0ug/mL <br> 100ug/mL | Induce Pyroptosis             |
| LPS + Nigericin   | 1 * 3 = (3)       | 1ug/mL + 1uM <br> 1ug/mL + 3uM <br> 1ug/mL + 10uM | Induce Pyroptosis             |


# Other Treatments and combinations

| Inducers       | # of Doses  | Dose Concentration(s)                                                        | Inhibitors | # of Doses | Dose Concentration(s) | Function                      |
| :----------    | :---------: | :---------                                                            | :--------  | :--------: |:---- |:---------                    |
| Media          | 1           | 0                                                                     | Media      | 1          | 0 | Control                       |
| DMSO           | 1           | 0.0250%                                                               | DMSO       | 2          | 0.025% <br> 1.0% | Vehicle Control               |
| DMSO           | 1           | 0.0250%                                                               | Z-VAD-FMK  | 2          | 30uM <br> 100 uM | Inhibit Inflamasome           |
| Thapsigargin   | 2           | 1uM <br> 10uM                                                         | DMSO       | 1          | 0.025% | ER Stress Induced Apoptosis   |
| H202           | 2           | 100nM <br> 100uM                                                      | DMSO       | 1          | 0.025% | Cell Death                    |
| H202           | 2           | 100nM <br> 100uM                                                      | Z-VAD-FMK  | 1          | 100uM | Inhibit Inflamasome           |
| H202           | 2           | 100nM <br> 100uM                                                      | Disulfiram | 1          | 1uM | Cell Death                    |
| Flagellin      | 2           | 0.1ug/mL <br> 1ug/mL                                                  | DMSO       | 1          | 0.025% | Induce Pyroptosis             |
| Flagellin      | 2           | 0.1ug/mL <br> 1ug/mL                                                  | Disulfiram | 1          | 1uM | Induce Pyroptosis <br> Inhibit Inflamasome          |
| LPS            | 5           | 0.01ug/mL <br> 0.1ug/mL <br> 1.0ug/mL <br> 10.0ug/mL <br> 100.0 ug/mL | DMSO       | 1          | 0.025% | Induce Pyroptosis             |
| LPS            | 1           | 10.0ug/mL                                                             | Z-VAD-FMK  | 1          | 100uM | Induce Pyroptosis / Inhibit Inflamasome          |
| LPS            | 1           | 10.0ug/mL                                                             | Disulfiram | 3          | 0.1uM <br> 1uM <br> 2.5uM | Induce Pyroptosis / Inhibit Inflamasome      |
| LPS + Nigercin | 2 * 3 = (6) | 1ug/mL+1ug/mL <br> 1ug/mL+3ug/mL <br> 1ug/mL+10ug/mL <br> 100ug/mL+1ug/mL <br> 100ug/mL+3ug/mL <br> 100ug/mL+10ug/mL | DMSO | 1 | 0.025% | Induce Pyroptosis |
| LPS + Nigercin | 1 * 1 = (1) | 1ug/mL+10ug/mL         | Z-VAD-FMK  | 1               | 100uM | Induce Pyroptosis <br> Inhibit Inflamasome          |
| LPS + Nigercin | 1 * 1 = (1) | 1ug/mL+10ug/mL        | Disulfiram | 1               | 1uM | Induce Pyroptosis <br> Inhibit Inflamasome          |
| Disulfiram     | 3           | 0.1uM <br> 1uM <br> 2.5uM       | DMSO       | 1               | 0.025% | Pyroptosis Inhibitor          |


