# Folder description

Scripts to generate srDFT functionals for Dalton.

* functionals.py; normal DFT functionals.
* functionals_mu.py; range-separated functionals.
* functionals_special.py; special case of some functionals. Might be changed later, to just be directly embedded into the functionals that include these special cases.
* print_functional_to_DALTON.py; Script to write functional and it derivatives to Dalton format.
* print_functional_to_DALTON_TPSSc_special_case.py; Script to write TPSSc to Dalton format, special case because of MAX() function. This should be changed later for shorter code and faster code generation.
* constants.txt; Constants for the functionals.
* constants_libxc.txt; Constants for the functionals, that are identical to those used in libxc.
* write_dalton_file.py; This is where all the ugly stuff is written in!


# Example of implementing code-generation of a functional (spin-PBE)

```python
E = mu_func.PBEc_mu(parameters).subs({rho_s: 0})*rho_c
d1E_rhoc = E.diff(rho_c)
d1E_gammacc = E.diff(gamma_cc)
d2E_rhoc2 = d1E_rhoc.diff(rho_c)
d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)

Kernel = [E]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1]
diff_idx = [[0]]
dalprint.dalton_functional_printer(Kernel, "ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
# # #
Kernel = [E,d1E_rhoc,d1E_gammacc]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1,2]
diff_idx = [[0],[1,3]]
dalprint.dalton_functional_printer(Kernel, "D1ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
# # #
Kernel = [E,d1E_rhoc,d1E_gammacc,d2E_rhoc2,d2E_rhocgammacc,d2E_gammacc2]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1,2,3]
diff_idx = [[0],[1,3],[1,4,6]]
dalprint.dalton_functional_printer(Kernel, "D2ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
```

## Specify the functional and create the derivatives. (having put the functionals into the right python functional files (just open the files to see how that works))

```python
E = mu_func.PBEc_mu(parameters).subs({rho_s: 0})*rho_c
d1E_rhoc = E.diff(rho_c)
d1E_gammacc = E.diff(gamma_cc)
d2E_rhoc2 = d1E_rhoc.diff(rho_c)
d2E_rhocgammacc = d1E_rhoc.diff(gamma_cc)
d2E_gammacc2 = d1E_gammacc.diff(gamma_cc)
```

## Specify all the values for the dalton_printer_function.

```python
Kernel = [E,d1E_rhoc,d1E_gammacc,d2E_rhoc2,d2E_rhocgammacc,d2E_gammacc2]
description = ["Implemented by E.R. Kjellgren.\n"]
diff_order = [1,2,3]
diff_idx = [[0],[1,3],[1,4,6]]
dalprint.dalton_functional_printer(Kernel, "D2ESRC_PBE_GWS_ERF", ["rho_c","gamma_cc"],["E","d1E","d2E"],description=description,shortrange=True,diff_order=diff_order,diff_idx=diff_idx, output_files=[daltonfile, f2pyfile])
```

* Kernel; a list of all the terms that is generated. Have to be order by differentiation order.
* description; some notes about the functional if wanted.
* diff_order; specify how many terms there is at every order of differentiation. [1,2,3] means one term for zeroth order. Two terms for first order. Three terms for second order.
* diff_idx; a list of lists of indexes to tell where to put the derivatives inside the Dalton vectors. [[0],[1,3],[1,4,6]] mean that for the zeroth order derivative it is put into place 0, zero meaning it is a scalar. [1,3] means that the first derivatives are put into index 1 and index 3 in the d1E() inside Dalton. 
* shortrange=True; specify that the functionals should have "mu" in their input. 
* output_files=[daltonfile, f2pyfile]; specify what files to write the output to.

The "diff_order" and "diff_idx" specification is a little tedious, but is needed because it is hard to know Python side what is supposed to be put where inside Dalton.


