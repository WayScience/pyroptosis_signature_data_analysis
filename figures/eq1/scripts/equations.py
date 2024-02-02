#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
from IPython.display import Math

# In[2]:


def latex2png(latex, filename, fontsize=300):
    """
    Render latex code to png image
    """
    from sympy import preview

    return preview(
        latex,
        viewer="file",
        filename=filename,
        euler=False,
        dvioptions=["-D", f"{str(fontsize)}"],
    )


# $$eq1: A \cap (B \cup C)^c$$
# $$eq2: B \cap (A \cup C)^c$$
# $$eq3: C \cap (A \cup B)^c$$
#
# $$ Where:
# \newline
# A = Apoptosis \; vs. \; Control
# \newline
# B = Apoptosis \; vs. \; Pyroptosis
# \newline
# C = Pyroptosis \; vs. \; Control$$

# In[3]:


# set output path
out_path = pathlib.Path("../figures").resolve()
# make sure output path exists
out_path.mkdir(parents=True, exist_ok=True)


# In[4]:


latex_string = r"""
$$\; eq1: A \cap (B \cup C)^c ; \; \; eq2: B \cap (A \cup C)^c ; \; \; eq3: C \cap (A \cup B)^c\; $$
$$Where: $$
$$A = Apoptosis \; vs. \; Control \; ANOVA \; features$$
$$B = Apoptosis \; vs. \; Pyroptosis \; ANOVA \; features$$
$$C = Pyroptosis \; vs. \; Control \; ANOVA \; features$$
"""

latex2png(latex_string, "../figures/equations.png", fontsize=300)
Math(latex_string)
