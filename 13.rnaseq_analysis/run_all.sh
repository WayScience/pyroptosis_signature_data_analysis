#!/bin/bash
jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb



conda activate rnaseq_r_env
cd notebooks || exit

papermill 0.SHSY5Y_expression_search.ipynb 0.SHSY5Y_expression_search.ipynb
papermill 1.search_whole_tissue.ipynb 1.search_whole_tissue.ipynb
papermill 2.plot_whole_tissue_search.ipynb 2.plot_whole_tissue_search.ipynb
papermill 3.U118mg_analysis.ipynb 3.U118mg_analysis.ipynb

cd ../ || exit

conda deactivate
