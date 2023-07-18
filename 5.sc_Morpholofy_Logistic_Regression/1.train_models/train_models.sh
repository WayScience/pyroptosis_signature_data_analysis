#!/bin/bash

papermill train_single_class.ipynb train_single_class.ipynb
papermill train_multi_class.ipynb train_multi_class.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

