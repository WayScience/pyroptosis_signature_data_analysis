#!/bin/bash

papermill test_single_class_model.ipynb test_single_class_model.ipynb
papermill test_multi_class_model.ipynb test_multi_class_model.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

