# pyroptosis_signature_data_analysis
Analysis Repository for the Interstellar project. For data processing for this project please see [pyroptosis_signature](https://github.com/WayScience/pyroptosis_signature)

### Workflow:

Cells were incubated for 1h in the inhibitors, and then incubated for 6h in each inducer.
For more information on plate layout and experimental design please see [pyroptosis_signature](https://github.com/WayScience/pyroptosis_signature).


# Treatments of Interest

| Treatments        | Number of Doses   | Function                      |
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

| Treatments        | Number of Doses   | Inhibitors | Number of Doses | Function                      |
| :----------       | :---------------: | :--------- | :----------:    | :-----------                  |
| Media             | 1                 | Media      | 1               | Control                       |
| DMSO              | 1                 | DMSO       | 2               | Vehicle Control               |
| DMSO              | 1                 | Z-VAD-FMK  | 2               | Inhibit Inflamasome           |
| Thapsigargin      | 2                 | DMSO       | 1               | ER Stress Induced Apoptosis   |
| H202              | 2                 | DMSO       | 1               | Cell Death                    |
| H202              | 2                 | Z-VAD-FMK  | 1               | Inhibit Inflamasome           |
| H202              | 2                 | Disulfiram | 1               | Cell Death                    |
| Flagellin         | 2                 | DMSO       | 1               | Induce Pyroptosis             |
| Flagellin         | 2                 | Disulfiram | 1               | Induce Pyroptosis / Inhibit Inflamasome          |
| LPS               | 5                 | DMSO       | 1               | Induce Pyroptosis             |
| LPS               | 1                 | Z-VAD-FMK  | 1               | Induce Pyroptosis / Inhibit Inflamasome          |
| LPS               | 1                 | Disulfiram | 3               | Induce Pyroptosis / Inhibit Inflamasome          |
| LPS + Nigercin    | 2 * 3 = (6)       | DMSO       | 1               | Induce Pyroptosis             |
| LPS + Nigercin    | 1 * 1 = (1)       | Z-VAD-FMK  | 1               | Induce Pyroptosis / Inhibit Inflamasome          |
| LPS + Nigercin    | 1 * 1 = (1)       | Disulfiram | 1               | Induce Pyroptosis / Inhibit Inflamasome          |
| Disulfiram        | 3                 | DMSO       | 1               | Pyroptosis Inhibitor          |


