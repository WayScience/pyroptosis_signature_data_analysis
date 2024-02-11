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
$$\textbf{eq.1:}\; A \cap (B \cup C)^c ; \; \; \textbf{eq.2:} \; B \cap (A \cup C)^c ; \; \; \textbf{eq.3:}\; C \cap (A \cup B)^c\; $$
$$Where: $$
$$A = Apoptosis \; vs. \; Control \; ANOVA \; features$$
$$B = Apoptosis \; vs. \; Pyroptosis \; ANOVA \; features$$
$$C = Pyroptosis \; vs. \; Control \; ANOVA \; features$$
"""

latex2png(latex_string, "../figures/equations.png", fontsize=300)
Math(latex_string)


# Equations 4-6

# ## $\tilde{A} = corr(Y,U) \; Where; \; Y=Dataset1 \; U$
# ## $\tilde{B} = corr(X,V)$
# ## $u_k = \frac{1}{P} \sum^P_{p=1} \tilde a^2_{pk}$
# ## $v_k = \frac{1}{Q} \sum^Q_{q=1} \tilde b^2_{qk}$
# ## $RI_u = u_k * r^2_k$
# ## $RI_v = v_k * r^2_k$

# In[5]:


latex_string = r"""
$$\textbf{eq.4:}\; X = \mathbb{R}^{N, P}; \; \textbf{eq.5:}\; Y = \mathbb{R}^{N, Q};\;\textbf{eq.6:}\; k = min(P, Q, N)$$
$$\textbf{eq.7:}\; \tilde{A} = corr(X,U); \; Where: \; U=Canonical \; variates \; for \; dataset\;X$$
$$\textbf{eq.8:}\; \tilde{B} = corr(Y,V); \; Where: \; V=Canonical \; variates \; for \; dataset\;Y$$


$$Where: \; \tilde{A} \; is \; \mathbb{R}^{P, k}; \; \tilde{B} \;is \;\mathbb{R}^{Q, k}$$
$$\textbf{eq.9:}\;u_k = \frac{1}{P} \sum^P_{p=1} \tilde a^2_{pk} \; \textbf{eq.10:}\;v_k = \frac{1}{Q} \sum^Q_{q=1} \tilde b^2_{qk}$$

$$Where: \; r^2_k \; is \; the \; Coefficient \; of \; determination \; for \; the \; k^{th} \; canonical \; variate$$

$$\textbf{eq.11:}\; RI_u = u_k * r^2_k; \;\; \textbf{eq.12:}\;\; RI_v = v_k * r^2_k$$
"""

latex2png(latex_string, "../figures/equations4_9.png", fontsize=300)
Math(latex_string)
